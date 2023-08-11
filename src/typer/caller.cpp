#include <cassert>       // assert
#include <limits>        // std::numeric_limits<T>
#include <memory>        // std::unique_ptr
#include <sstream>       // std::ostringstream
#include <string>        // std::string
#include <unordered_map> // std::unordered_map
#include <vector>        // std::vector

#include <parallel_hashmap/phmap.h>

#include <paw/align.hpp>
#include <paw/station.hpp>

#include <seqan/basic.h>
#include <seqan/hts_io.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/absolute_position.hpp> // gyper::absolute_pos
#include <graphtyper/graph/constructor.hpp>
#include <graphtyper/graph/graph_serialization.hpp> // gyper::load_graph
#include <graphtyper/graph/haplotype_calls.hpp>
#include <graphtyper/graph/haplotype_extractor.hpp>
#include <graphtyper/graph/reference_depth.hpp>
#include <graphtyper/index/indexer.hpp> // gyper::index (global)
#include <graphtyper/index/ph_index.hpp>
#include <graphtyper/typer/alignment.hpp>
#include <graphtyper/typer/bucket.hpp>
#include <graphtyper/typer/caller.hpp>
#include <graphtyper/typer/event.hpp>
#include <graphtyper/typer/primers.hpp> // gyper::Primers
#include <graphtyper/typer/read.hpp>
#include <graphtyper/typer/variant_candidate.hpp>
#include <graphtyper/typer/variant_map.hpp>
#include <graphtyper/typer/vcf.hpp>
#include <graphtyper/utilities/hts_parallel_reader.hpp> // gyper::HtsParallelReader
#include <graphtyper/utilities/hts_reader.hpp>          // gyper::HtsReader
#include <graphtyper/utilities/io.hpp>
#include <graphtyper/utilities/logging.hpp> // BOOST_LOG_TRIVIAL
#include <graphtyper/utilities/options.hpp> // gyper::Options

namespace gyper
{
struct HaplotypeInfo
{
  phmap::flat_hash_set<gyper::Event, gyper::EventHash> ever_together;
  phmap::flat_hash_set<gyper::Event, gyper::EventHash> always_together;
};

} // namespace gyper

using Thap = gyper::HaplotypeInfo;

namespace
{
#ifndef NDEBUG
std::string const debug_read_name = "HISEQ1:33:H9YY4ADXX:1:2110:2792:58362/2";
long debug_event_pos{3603668};
char debug_event_type{'X'};
std::size_t debug_event_size{1};
#endif // NDEBUG

void merge_haplotypes2(std::map<gyper::Event, Thap> & into, std::map<gyper::Event, Thap> & from)
{
  // check trivial case
  if (into.size() == 0)
  {
    into = std::move(from);
    from.clear();
    return;
  }

  // Insert any elements from "from" into "into"
  for (auto from_it = from.begin(); from_it != from.end(); ++from_it)
  {
#ifndef NDEBUG
    if (from_it->first.pos == debug_event_pos)
      gyper::print_log(gyper::log_severity::info, __HERE__, " var ", from_it->first.to_string());
#endif // NDEBUG

    auto insert_p = into.insert(*from_it);

    if (insert_p.second)
    {
      // we hadn't seen this event before in "into"
      // nothing needs to be done for "ever_together"
      // go through all events in the newly inserted "always_together" and remove any that have been seen in into

#ifndef NDEBUG
      if (insert_p.first->first.pos == debug_event_pos)
        gyper::print_log(gyper::log_severity::info, __HERE__, " new var ", insert_p.first->first.to_string());
#endif // NDEBUG

      gyper::HaplotypeInfo & into_hap_info = insert_p.first->second;
      auto it = into_hap_info.always_together.begin();

      for (; it != into_hap_info.always_together.end();) // no increment
      {
        if (into.count(*it) > 0)
        {
          // gyper::print_log(gyper::log_severity::info, __HERE__, " Found old event ", it->to_string());
          it = into_hap_info.always_together.erase(it);
        }
        else
        {
          // gyper::print_log(gyper::log_severity::info, __HERE__, " Keeping ", it->to_string());
          ++it;
        }
      }
    }
    else
    {
      // the event was already in "into"
      gyper::HaplotypeInfo & into_hap_info = insert_p.first->second;
      gyper::HaplotypeInfo & from_hap_info = from_it->second;

      // The event has been seen before, add ever_together events
      std::move(from_hap_info.ever_together.begin(),
                from_hap_info.ever_together.end(),
                std::inserter(into_hap_info.ever_together, into_hap_info.ever_together.begin()));

      // take the intersection for the always_together events
      phmap::flat_hash_set<gyper::Event, gyper::EventHash> intersection;

      // go through the smaller table
      if (from_hap_info.always_together.size() <= into_hap_info.always_together.size())
      {
        for (auto const & at : from_hap_info.always_together)
        {
          if (into_hap_info.always_together.count(at) > 0)
            intersection.insert(at);
        }
      }
      else
      {
        for (auto const & at : into_hap_info.always_together)
        {
          if (from_hap_info.always_together.count(at) > 0)
            intersection.insert(at);
        }
      }

#ifndef NDEBUG
      if (insert_p.first->first.pos == debug_event_pos)
      {
        gyper::print_log(gyper::log_severity::info, __HERE__, " intersection of ", insert_p.first->first.to_string());
        gyper::print_log(gyper::log_severity::info,
                         __HERE__,
                         " ",
                         from_hap_info.always_together.size(),
                         " ",
                         into_hap_info.always_together.size(),
                         " ",
                         intersection.size());
      }
#endif // NDEBUG

      insert_p.first->second.always_together = std::move(intersection);
    }
  }

  // Clear memory asap
  from.clear();
}

bool is_clipped(bam1_t const & b, uint32_t const min_count = 1)
{
  if (b.core.n_cigar == 0)
    return false;

  auto it = b.data + b.core.l_qname;

  // Check first
  uint32_t opAndCnt;
  memcpy(&opAndCnt, it, sizeof(uint32_t));
  uint32_t cigar_count = opAndCnt >> 4;

  if ((opAndCnt & 15) == 4 && cigar_count >= min_count)
  {
    return true;
  }

  // Check last
  memcpy(&opAndCnt, it + sizeof(uint32_t) * (b.core.n_cigar - 1), sizeof(uint32_t));
  cigar_count = opAndCnt >> 4;

  if ((opAndCnt & 15) == 4 && cigar_count >= min_count)
  {
    // std::cerr << "MIDNSHP=X*******"[opAndCnt & 15] << " last ";
    return true;
  }

  return false;
}

void _determine_num_jobs_and_num_parts(long const & jobs, long & num_parts, long const NUM_SAMPLES)
{
  using namespace gyper;

  gyper::Options const & copts = *(Options::const_instance());
  num_parts = jobs;
  long const MAX_FILES_OPEN = copts.max_files_open > jobs ? copts.max_files_open : jobs;

  if (jobs >= NUM_SAMPLES)
  {
    // Special case where there are more threads than samples OR more threads than max files open.
    // In this case, each threads gets 1 sample
    num_parts = std::min(NUM_SAMPLES, MAX_FILES_OPEN);
  }
  else if (NUM_SAMPLES > MAX_FILES_OPEN)
  {
    // More samples than maximum allowed files open at the same time. Pick num_parts such that
    //  ceil(NUM_SAMPLES / num_parts) * jobs <= MAX_FILES_OPEN
    long const MAX_FILES_OPEN_PER_THREAD = (MAX_FILES_OPEN + jobs - 1) / jobs;
    num_parts = (NUM_SAMPLES + MAX_FILES_OPEN_PER_THREAD - 1) / MAX_FILES_OPEN_PER_THREAD;
    // jobs unchanged
  }
  // else num_parts and jobs unchanged
}

} // namespace

namespace gyper
{
std::vector<std::string> call(
  std::vector<std::string> const & hts_paths,
  std::vector<double> const & avg_cov_by_readlen,
  std::string const & graph_path,
  PHIndex const & ph_index,
  std::string const & output_dir,
  std::string const & reference_fn,
  std::string const & region,
  Primers const * primers,
  std::map<std::pair<uint16_t, uint16_t>, std::map<std::pair<uint16_t, uint16_t>, int8_t>> & ph,
  bool const is_writing_calls_vcf,
  bool const is_writing_hap,
  std::vector<std::unordered_map<uint32_t, uint32_t>> * allele_hap_gts_ptr)
{
  if (hts_paths.size() == 0)
  {
    print_log(log_severity::error, __HERE__, " No input BAM/CRAM files.");
    std::exit(1);
  }

  if (graph_path.size() > 0)
  {
    load_graph(graph_path); // Loads the graph into the global variable 'graph'
  }

  using TMap = std::map<std::pair<uint16_t, uint16_t>, std::map<std::pair<uint16_t, uint16_t>, int8_t>>;

  // Split hts_paths
  std::vector<std::string> paths;
  std::vector<std::unique_ptr<std::vector<std::string>>> spl_hts_paths;
  std::vector<std::unique_ptr<std::vector<double>>> spl_avg_cov;
  std::vector<TMap> spl_ph;
  assert(Options::const_instance()->max_files_open > 0);

  long jobs{Options::const_instance()->threads}; // Running jobs
  long num_parts{1};                             // Number of parts to split the work into
  long const NUM_SAMPLES{static_cast<long>(hts_paths.size())};

  _determine_num_jobs_and_num_parts(jobs, num_parts, NUM_SAMPLES);

  if (jobs > num_parts)
    jobs = num_parts;

  assert(jobs >= 1);
  spl_ph.resize(jobs); // spl_ph per thread

  auto emplace_paths = [&spl_hts_paths, &spl_avg_cov](std::vector<std::string>::const_iterator & hts_paths_it,
                                                      std::vector<double>::const_iterator & cov_it,
                                                      long const num_parts,
                                                      long const num_parts_start,
                                                      long const num_parts_total,
                                                      long const num_samples,
                                                      long const num_splits)
  {
    assert(num_parts > 0);
    assert(num_samples > 0);
    assert(num_splits > 0);

    for (long i{num_parts_start}; i < (num_parts + num_parts_start); ++i)
    {
      long const advance = num_samples / num_parts_total + (i < (num_samples % num_parts_total));
      assert(advance > 0);

      auto next_it = hts_paths_it;
      auto cov_next_it = cov_it;

      for (long j{0}; j < num_splits; ++j)
      {
        long const split_advance = (advance / num_splits) + (j < (advance % num_splits));

        // ..just in case, no point in having empty work
        if (split_advance == 0)
          continue;

        std::advance(next_it, split_advance);
        std::advance(cov_next_it, split_advance);

        spl_hts_paths.emplace_back(new std::vector<std::string>(hts_paths_it, next_it));
        hts_paths_it = next_it;

        spl_avg_cov.emplace_back(new std::vector<double>(cov_it, cov_next_it));
        cov_it = cov_next_it;
      }
    }
  };

  // jobs=4, NUM_SAMPLES=1000, num_parts=4
  if (jobs <= 2 || NUM_SAMPLES <= 20 || NUM_SAMPLES < (4 * jobs))
  {
    // special case where there are 4 or less samples allocated to each thread (or we have few threads..)
    auto hts_paths_it = hts_paths.begin();
    auto cov_it = avg_cov_by_readlen.begin();

    emplace_paths(hts_paths_it, cov_it, num_parts, 0, num_parts, NUM_SAMPLES, 1);

    // make sure everything is consumed
    assert(hts_paths_it == hts_paths.end());
    assert(cov_it == avg_cov_by_readlen.end());
  }
  else if (num_parts < (4 * jobs))
  {
    // BOOST_LOG_TRIVIAL(info) << __HERE__ << " Case 2: Some threads get less than 4 \"work packages\""
    //                        << " and there are more than 4 samples allocated to each thread";

    // Some threads do not get two "work packages" and there are more than 4 samples allocated to each thread
    long const num_samples_first = NUM_SAMPLES / 2;
    long current_phase_parts = jobs;
    assert(jobs > 1);
    _determine_num_jobs_and_num_parts(jobs - 1, current_phase_parts, num_samples_first); // skip main thread

    auto hts_paths_it = hts_paths.begin();
    auto cov_it = avg_cov_by_readlen.begin();

    //#ifndef NDEBUG
    //    BOOST_LOG_TRIVIAL(info) << __HERE__ << " " <<
    //#endif // NDEBUG

    emplace_paths(hts_paths_it, cov_it, current_phase_parts, 0, current_phase_parts, num_samples_first, 1);

    assert(hts_paths_it == (hts_paths.begin() + num_samples_first));
    assert(cov_it == (avg_cov_by_readlen.begin() + num_samples_first));

    long const num_samples_second = NUM_SAMPLES / 4;

    if (num_samples_second > 0)
    {
      current_phase_parts = jobs;
      assert(jobs > 1);
      _determine_num_jobs_and_num_parts(jobs, current_phase_parts, num_samples_second);

      emplace_paths(hts_paths_it, cov_it, current_phase_parts, 0, current_phase_parts, num_samples_second, 2);

      assert(hts_paths_it == (hts_paths.begin() + num_samples_first + num_samples_second));
      assert(cov_it == (avg_cov_by_readlen.begin() + num_samples_first + num_samples_second));
    }

    // long const num_parts_remaining = num_parts - num_parts_three_quarters;
    long const num_samples_remaining = NUM_SAMPLES - num_samples_first - num_samples_second;
    assert(num_samples_remaining > 0);
    _determine_num_jobs_and_num_parts(jobs, current_phase_parts, num_samples_remaining);

    emplace_paths(hts_paths_it, cov_it, current_phase_parts, 0, current_phase_parts, num_samples_remaining, 4);

    // make sure everything is consumed
    assert(hts_paths_it == hts_paths.end());
    assert(cov_it == avg_cov_by_readlen.end());
  }
  else
  {
    // BOOST_LOG_TRIVIAL(info) << __HERE__ << " Case 3: We reduce the size of the last two packages "
    //                        << "by half when there are threads*2 packages left";
    // There are more than 4 sample per thread and more jobs than threads
    // In this case we reduce the size of the last two "work packages" by half when there are threads*2 packages left
    auto hts_paths_it = hts_paths.begin();
    auto cov_it = avg_cov_by_readlen.begin();
    assert(num_parts >= (2 * jobs));
    long const first_phase_parts = num_parts - 2 * jobs;

    emplace_paths(hts_paths_it, cov_it, first_phase_parts, 0, num_parts, NUM_SAMPLES, 1);
    emplace_paths(hts_paths_it, cov_it, jobs, first_phase_parts, num_parts, NUM_SAMPLES, 2);
    emplace_paths(hts_paths_it, cov_it, jobs, first_phase_parts + jobs, num_parts, NUM_SAMPLES, 4);

    // make sure everything is consumed
    assert(hts_paths_it == hts_paths.end());
    assert(cov_it == avg_cov_by_readlen.end());
  }

  print_log(log_severity::debug, __HERE__, " Number of pools = ", spl_hts_paths.size());
  long const NUM_POOLS = spl_hts_paths.size();
  paths.resize(NUM_POOLS);

  // Run in parallel
  {
    paw::Station call_station(jobs); // last parameter is queue_size

    for (long i{0}; i < (NUM_POOLS - 1l); ++i)
    {
      call_station.add_work_with_thread_id(parallel_reader_genotype_only,
                                           &paths[i],
                                           spl_hts_paths[i].get(),
                                           spl_avg_cov[i].get(),
                                           &output_dir,
                                           &reference_fn,
                                           &region,
                                           &ph_index,
                                           primers,
                                           &spl_ph,
                                           is_writing_calls_vcf,
                                           is_writing_hap,
                                           allele_hap_gts_ptr);
    }

    // Do the last pool on the current thread
    call_station.add_to_thread(jobs - 1,
                               parallel_reader_genotype_only,
                               jobs - 1,
                               &paths[NUM_POOLS - 1],
                               spl_hts_paths[NUM_POOLS - 1].get(),
                               spl_avg_cov[NUM_POOLS - 1].get(),
                               &output_dir,
                               &reference_fn,
                               &region,
                               &ph_index,
                               primers,
                               &spl_ph,
                               is_writing_calls_vcf,
                               is_writing_hap,
                               allele_hap_gts_ptr);

    std::string thread_info = call_station.join();
    print_log(log_severity::info, "Finished calling. Thread work: ", thread_info);
  }

  if (is_writing_hap)
  {
    // gather ph
    ph = std::move(spl_ph[0]);

    for (long i{1}; i < static_cast<long>(spl_ph.size()); ++i)
    {
      // merge
      std::map<std::pair<uint16_t, uint16_t>, std::map<std::pair<uint16_t, uint16_t>, int8_t>> & ph_parts = spl_ph[i];

      for (auto && ph_part : ph_parts)
      {
#ifndef NDEBUG
        std::size_t const old_size = ph_part.second.size();
#endif // NDEBUG
        auto insert1_it = ph.insert(std::move(ph_part));

        if (insert1_it.second)
        {
          continue; // not found
        }
#ifndef NDEBUG
        assert(ph_part.second.size() == old_size);
#endif // NDEBUG

        std::map<std::pair<uint16_t, uint16_t>, int8_t> & ph2_map = ph_part.second;

        // found - so it was not inserted
        for (auto const & ph2 : ph2_map)
        {
          int8_t const new_flags = ph2.second;
          auto insert2_it = insert1_it.first->second.insert(ph2);

          if (insert2_it.second)
            continue; // not found

          int8_t & flags = insert2_it.first->second;
          flags |= new_flags;
        }

        ph2_map.clear();
      }
    }
  }

  // BOOST_LOG_TRIVIAL(info) << "Finished a calling pass for all samples.";
  return paths;
}

void run_first_pass(bam1_t * hts_rec,
                    HtsReader & hts_reader,
                    long file_i,
                    std::vector<BucketFirstPass> & buckets,
                    std::map<Event, Thap> & pool_haplotypes,
                    long const BUCKET_SIZE,
                    long const region_begin,
                    std::vector<char> const & reference_sequence)
{
  long const REF_SIZE{static_cast<long>(reference_sequence.size())};
  int32_t global_max_pos_end{0};
  std::vector<uint32_t> cov_up(REF_SIZE);
  std::vector<uint32_t> cov_down(REF_SIZE);
  std::map<Event, Thap> sample_haplotypes;

  // make sure the buckets are empty
#ifndef NDEBUG
  for (auto & bucket : buckets)
  {
    assert(bucket.events.size() == 0);
  }
#endif // NDEBUG

  while (true)
  {
    assert(hts_rec);

    while (hts_rec->core.n_cigar == 0 || hts_rec->core.pos < region_begin)
    {
      hts_rec = hts_reader.get_next_read(hts_rec);

      if (!hts_rec)
        break;
    }

    if (!hts_rec)
      break;

    std::array<char, 16> static constexpr CIGAR_MAP = {
      {'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', 'B', '*', '*', '*', '*', '*', '*'}};

    auto const & core = hts_rec->core;
    auto it = hts_rec->data;
    auto cigar_it = it + core.l_qname;
    auto seq_it = cigar_it + (core.n_cigar << 2);

    // Index of first aligned read base in read sequence
    long read_offset{0};

    // Index of first aligned read base in reference_sequence
    long ref_offset{static_cast<long>(core.pos) - region_begin};
    assert(ref_offset >= 0);

    long const bucket_index{ref_offset / BUCKET_SIZE};

    if (bucket_index >= static_cast<long>(buckets.size()))
      buckets.resize(bucket_index + 1);

    assert(bucket_index >= 0);
    assert(bucket_index < static_cast<long>(buckets.size()));
    long const N_CIGAR = core.n_cigar;

    if (ref_offset >= REF_SIZE)
    {
      print_log(log_severity::warning, __HERE__, " Unexpected ref_offset=", ref_offset, " >= REF_SIZE=", REF_SIZE);
      break;
    }

    Read read;
    read.name = reinterpret_cast<char *>(it);
    // read.mate_pos = static_cast<int32_t>(core.mpos);

    if ((core.flag & IS_FIRST_IN_PAIR) != 0u)
      read.name.append("/1");
    else
      read.name.append("/2");

    read.flags = core.flag;
    read.sequence.resize(core.l_qseq);

    for (int i = 0; i < core.l_qseq; ++i)
      read.sequence[i] = seq_nt16_str[bam_seqi(seq_it, i)];

    auto const qual_it = bam_get_qual(hts_rec);
    std::vector<std::map<Event, EventSupport>::iterator> cigar_events;
    bool const is_read_clipped = is_clipped(*hts_rec);

    // Check if

    for (long i{0}; i < N_CIGAR; ++i)
    {
      uint32_t opAndCnt;
      memcpy(&opAndCnt, cigar_it, sizeof(uint32_t));
      cigar_it += sizeof(uint32_t);
      uint32_t const cigar_count = opAndCnt >> 4;
      char const cigar_operation = CIGAR_MAP[opAndCnt & 15];
      assert(cigar_count > 0);

      if (ref_offset >= REF_SIZE)
        break;

      switch (cigar_operation)
      {
      case 'M':
      case '=': // '=' and 'X' are typically not used by aligners (at least not by default),
      case 'X': // but we keep it here just in case
      {
        for (long r{0}; r < cigar_count; ++r)
        {
          long const ref_pos = ref_offset + r;

          if (ref_pos >= REF_SIZE)
            break;

          char const ref = reference_sequence[ref_pos];
          long const read_pos = read_offset + r;

          if (read_pos >= static_cast<long>(read.sequence.size()))
            break;

          assert(read_pos < static_cast<long>(read.sequence.size()));
          char const read_base = read.sequence[read_pos];

          if (read_base == ref || (ref != 'A' && ref != 'C' && ref != 'G' && ref != 'T') ||
              (read_base != 'A' && read_base != 'C' && read_base != 'G' && read_base != 'T'))
          {
            continue;
          }

          Event new_snp_event(ref_pos + region_begin, 'X', {read_base});

#ifndef NDEBUG
          if (debug_event_type == 'X' && new_snp_event.pos == debug_event_pos)
          {
            print_log(log_severity::info,
                      __HERE__,
                      " ",
                      new_snp_event.to_string(),
                      " in file_i=",
                      file_i,
                      " core.pos=",
                      core.pos,
                      " read=",
                      read.name);
          }
#endif // NDEBUG

          std::map<Event, EventSupport>::iterator snp_it =
            add_snp_event_to_bucket(buckets, std::move(new_snp_event), region_begin, BUCKET_SIZE);

          EventSupport & event_support = snp_it->second;
          uint8_t const qual = *(qual_it + read_pos);

          if (qual >= 25)
            ++event_support.hq_count;
          else
            ++event_support.lq_count;

          if (core.qual != 255 && core.qual > event_support.max_mapq)
            event_support.max_mapq = core.qual;

          event_support.proper_pairs += ((read.flags & IS_PROPER_PAIR) != 0);
          event_support.first_in_pairs += ((read.flags & IS_FIRST_IN_PAIR) != 0);
          event_support.sequence_reversed += ((read.flags & IS_SEQ_REVERSED) != 0);
          event_support.clipped += is_read_clipped;

          if (event_support.uniq_pos1 == -1)
          {
            event_support.uniq_pos1 = core.pos;
          }
          else if (event_support.uniq_pos2 == -1)
          {
            if (event_support.uniq_pos1 != core.pos)
            {
              // due to something in bamshrink, the input is not always sorted
              // assert(core.pos > event_support.uniq_pos1);
              event_support.uniq_pos2 = core.pos;
            }
          }
          else if (event_support.uniq_pos3 == -1 &&
                   event_support.uniq_pos2 != core.pos /*&& event_support.uniq_pos1 != core.pos*/)
          {
            assert(event_support.uniq_pos1 != core.pos); // should not happend if input is sorted
            event_support.uniq_pos3 = core.pos;
          }

          long const max_distance = std::min(read_pos, core.l_qseq - 1 - read_pos);
          assert(max_distance >= 0);

          if (max_distance > event_support.max_distance)
            event_support.max_distance = max_distance;

          cigar_events.push_back(snp_it);
        }

        read_offset += cigar_count;
        ref_offset += cigar_count;
        break;
      }

      case 'I':
      {
        assert(cigar_count > 0);

        auto begin_it = read_offset < static_cast<long>(read.sequence.size()) ? read.sequence.begin() + read_offset
                                                                              : read.sequence.end();

        auto end_it = (read_offset + cigar_count) < static_cast<long>(read.sequence.size())
                      ? read.sequence.begin() + read_offset + cigar_count
                      : read.sequence.end();

        if (begin_it == end_it)
          break;

        // Make sure all bases are ACGT
        if (std::all_of(begin_it, end_it, [](char c) { return c == 'A' || c == 'C' || c == 'G' || c == 'T'; }))
        {
          std::vector<char> event_sequence(begin_it, end_it);
          Event new_event = make_insertion_event(region_begin + ref_offset, std::move(event_sequence));

          // Add to bucket events
          auto indel_event_it = add_indel_event_to_bucket(buckets,
                                                          std::move(new_event),
                                                          region_begin,
                                                          BUCKET_SIZE,
                                                          reference_sequence,
                                                          ref_offset);

#ifndef NDEBUG
          if (read.name == debug_read_name)
            print_log(log_severity::info, __HERE__, " new ins=", indel_event_it->first.to_string());
#endif // NDEBUG

          auto & event_support = indel_event_it->second;
          ++event_support.hq_count;

          if (core.qual != 255 && core.qual > event_support.max_mapq)
            event_support.max_mapq = core.qual;

          event_support.proper_pairs += ((read.flags & IS_PROPER_PAIR) != 0);
          // event_support.first_in_pairs += ((read.flags & IS_FIRST_IN_PAIR) != 0);
          event_support.sequence_reversed += ((read.flags & IS_SEQ_REVERSED) != 0);
          event_support.clipped += is_read_clipped;
          cigar_events.push_back(indel_event_it);
        }

        read_offset += cigar_count;
        break;
      }

      case 'D':
      {
        assert(cigar_count > 0);

        if (ref_offset + cigar_count >= REF_SIZE)
        {
          ref_offset += cigar_count;
          break;
        }

        Event new_event = make_deletion_event(reference_sequence, ref_offset, region_begin + ref_offset, cigar_count);

        if (std::all_of(new_event.sequence.begin(),
                        new_event.sequence.end(),
                        [](char c) { return c == 'A' || c == 'C' || c == 'G' || c == 'T'; }))
        {
          auto indel_event_it = add_indel_event_to_bucket(buckets,
                                                          std::move(new_event),
                                                          region_begin,
                                                          BUCKET_SIZE,
                                                          reference_sequence,
                                                          ref_offset);
#ifndef NDEBUG
          if (read.name == debug_read_name)
            print_log(log_severity::info, __HERE__, " new del=", indel_event_it->first.to_string());
#endif // NDEBUG

          auto & event_support = indel_event_it->second;
          ++event_support.hq_count;

          if (core.qual != 255 && core.qual > event_support.max_mapq)
            event_support.max_mapq = core.qual;

          event_support.proper_pairs += ((read.flags & IS_PROPER_PAIR) != 0);
          // event_support.first_in_pairs += ((read.flags & IS_FIRST_IN_PAIR) != 0);
          event_support.sequence_reversed += ((read.flags & IS_SEQ_REVERSED) != 0);
          event_support.clipped += is_read_clipped;
          cigar_events.push_back(indel_event_it);
        }

        ref_offset += cigar_count;
        break;
      }

      case 'S':
      {
        read_offset += cigar_count;
        break;
      } // else 'H' or 'P' which move neither reference nor read
      }
    }

    // TODO set as options
    int constexpr HIGH_EVENT_COUNT{12};
    int constexpr VHIGH_EVENT_COUNT{18};

    if (static_cast<int>(cigar_events.size()) >= HIGH_EVENT_COUNT)
    {
      for (auto event_it = cigar_events.begin(); event_it != cigar_events.end(); ++event_it)
      {
        EventSupport & info = (*event_it)->second;

        if (cigar_events.size() >= VHIGH_EVENT_COUNT)
        {
          // very high amount of events
          if (info.hq_count > 0)
            --info.hq_count;
          else if (info.lq_count > 0)
            --info.lq_count;
        }
        else
        {
          if (info.hq_count > 0)
          {
            --info.hq_count;
            ++info.lq_count;
          }
        }
      }
    }

    if (static_cast<int>(cigar_events.size()) < VHIGH_EVENT_COUNT)
    {
      for (long e{1}; e < static_cast<long>(cigar_events.size()); ++e)
      {
        Event const & event = cigar_events[e]->first;

        // Add event to previous snp events
        for (long prev{0}; prev < e; ++prev)
        {
          auto & prev_info = cigar_events[prev]->second;
          auto insert_it = prev_info.phase.insert({event, 0});
          ++insert_it.first->second;
        }
      }
    }

    read.alignment.pos = static_cast<int32_t>(core.pos);
    read.alignment.pos_end = region_begin + std::min(ref_offset, REF_SIZE - 1);

    assert((read.alignment.pos - region_begin) >= 0);
    assert((read.alignment.pos - region_begin) < static_cast<long>(cov_up.size()));
    ++cov_up[read.alignment.pos - region_begin];

    assert((read.alignment.pos_end - region_begin) >= 0);
    assert((read.alignment.pos_end - region_begin) < static_cast<long>(cov_down.size()));
    ++cov_down[read.alignment.pos_end - region_begin];

    auto & bucket = buckets[bucket_index];
    int32_t const end_with_clip = read.alignment.pos_end + read.alignment.num_clipped_end;

    if (end_with_clip > bucket.max_pos_end)
    {
      bucket.max_pos_end = end_with_clip;
      global_max_pos_end = std::max(global_max_pos_end, end_with_clip);
    }

    bucket.global_max_pos_end = global_max_pos_end;

#ifndef NDEBUG
    // test if the input is sorted
    long prev_pos = hts_rec->core.pos;
#endif // NDEBUG
    hts_rec = hts_reader.get_next_read(hts_rec);

    if (!hts_rec)
      break;

#ifndef NDEBUG
    if (hts_rec->core.pos < prev_pos)
    {
      print_log(log_severity::warning,
                __HERE__,
                " file_i=",
                file_i,
                " is not sorted. ",
                hts_rec->core.pos,
                " < ",
                prev_pos);
      assert(hts_rec->core.pos >= prev_pos);
    }
#endif
  }

  // Check if we have too many buckets
  if ((static_cast<long>(buckets.size()) - 1l) * BUCKET_SIZE >= REF_SIZE)
  {
    long const new_size = ((REF_SIZE - 1) / BUCKET_SIZE) + 1;
    assert(new_size < static_cast<long>(buckets.size()));
    buckets.resize(new_size);
  }

  long const NUM_BUCKETS{static_cast<long>(buckets.size())};

  auto update_coverage = [&cov_up, &cov_down](long & cov, long const pos, long const b, long const BUCKET_SIZE) -> void
  {
    long offset{pos + 1};

    if (offset > b * BUCKET_SIZE)
    {
      offset = b * BUCKET_SIZE;

      // Add coverage until naive begin
      while (offset <= pos)
      {
        cov += (static_cast<long>(cov_up[offset]) - static_cast<long>(cov_down[offset]));
        ++offset;
      }
    }
  };

  // Remove SNPs with low support
  {
    long depth{0};

    for (long b{0}; b < NUM_BUCKETS; ++b)
    {
      auto & bucket = buckets[b];

      for (auto snp_it = bucket.events.begin(); snp_it != bucket.events.end();) // no increment
      {
        Event const & snp = snp_it->first;

        if (snp.type != 'X')
        {
          ++snp_it;
          continue; // not a SNP
        }

        EventSupport const & info = snp_it->second;

        long cov{depth};
        long const begin = std::max(0l, snp.pos - region_begin);
        update_coverage(cov, begin, b, BUCKET_SIZE);

        if (info.has_good_support(cov))
        {
#ifndef NDEBUG
          if (debug_event_type == 'X' && snp.pos == debug_event_pos)
          {
            print_log(log_severity::info,
                      __HERE__,
                      " debug snp has good support ",
                      snp.to_string(),
                      " ",
                      info.to_string(),
                      " file_i=",
                      file_i);
          }
#endif // DEBUG
          ++snp_it;
        }
        else
        {
#ifndef NDEBUG
          if (debug_event_type == 'X' && snp.pos == debug_event_pos)
          {
            print_log(log_severity::info,
                      __HERE__,
                      " debug snp has bad support ",
                      snp.to_string(),
                      " ",
                      info.to_string(),
                      " file_i=",
                      file_i);
          }
#endif // DEBUG

          snp_it = bucket.events.erase(snp_it);
        }
      }

      // Add coverage from bucket b
      if ((b * BUCKET_SIZE) >= REF_SIZE)
        break;

      assert((b * BUCKET_SIZE) < REF_SIZE);
      long offset{b * BUCKET_SIZE};
      long const end_offset = std::min(REF_SIZE, (b + 1) * BUCKET_SIZE);

      while (offset < end_offset)
      {
        depth += static_cast<long>(cov_up[offset]) - static_cast<long>(cov_down[offset]);
        ++offset;
      }
    }
  }

  long depth{0};

  // Find indels with good enough support to allow realignments, destroy the others
  for (long b{0}; b < NUM_BUCKETS; ++b)
  {
    auto & bucket = buckets[b];

    // Remove bad indels
    for (Tindel_events::iterator it = bucket.events.begin(); it != bucket.events.end();) // no increment
    {
      Event const & indel = it->first;

      if (indel.type == 'X')
      {
        ++it;
        continue; // do not remove SNPs here
      }

      assert(indel.type == 'I' || indel.type == 'D');
      EventSupport const & indel_info = it->second;
      assert(!indel_info.has_realignment_support);

      long const naive_pad = 4.0 + static_cast<double>(indel.sequence.size()) / 3.0;
      long const naive_begin = std::max(0l, indel.pos - naive_pad - region_begin);
      long const naive_end = std::min(REF_SIZE, indel.pos + indel_info.span + naive_pad - region_begin);

      double const correction = indel.type == 'I' ? static_cast<double>(indel.sequence.size() / 2.0 + 8.0) / 8.0
                                                  : static_cast<double>(indel.sequence.size() / 3.0 + 10.0) / 10.0;

      double const count = correction * (indel_info.hq_count + indel_info.lq_count);

#ifndef NDEBUG
      if (debug_event_type == indel.type && debug_event_pos == indel.pos && debug_event_size == indel.sequence.size())
      {
        print_log(log_severity::info,
                  __HERE__,
                  "Indel [pos,pos+span], size [begin,end]: ",
                  indel.to_string(),
                  " [",
                  region_begin,
                  " ",
                  indel.pos,
                  ",",
                  (indel.pos + indel_info.span),
                  "] ",
                  indel.sequence.size(),
                  " ",
                  correction,
                  " [",
                  naive_begin,
                  ",",
                  naive_end,
                  "]",
                  " count=",
                  count,
                  " depth=",
                  depth);
      }
#endif // NDEBUG

      long cov{depth}; // starting coverage depth
      long offset{naive_begin};

      // Fix coverage starting before here
      if (offset <= b * BUCKET_SIZE)
      {
        while (offset < b * BUCKET_SIZE)
        {
          cov -= (static_cast<long>(cov_up[offset]) - static_cast<long>(cov_down[offset]));
          ++offset;
        }
      }
      else
      {
        offset = b * BUCKET_SIZE;

        // Add coverage until naive begin
        while (offset < naive_begin)
        {
          cov += (static_cast<long>(cov_up[offset]) - static_cast<long>(cov_down[offset]));
          ++offset;
        }
      }

      // reduce coverage by reads who do not overlap the whole naive interval
      while (offset <= naive_end)
      {
        cov -= static_cast<long>(cov_down[offset]);
        ++offset;
      }

      double const corrected_cov{std::max(static_cast<double>(cov), count)};
      double const anti_count_d{corrected_cov - count};
      uint32_t log_qual = get_log_qual_double(count, anti_count_d, 10.0);

      if (indel_info.hq_count >= 6 && count >= 8.0 && log_qual >= 60 && indel_info.sequence_reversed > 0 &&
          indel_info.sequence_reversed < indel_info.hq_count && indel_info.proper_pairs >= 3 &&
          indel_info.max_mapq >= 20 && (indel_info.clipped == 0 || (indel_info.clipped + 3) <= indel_info.hq_count))
      {
#ifndef NDEBUG
        if (debug_event_type == indel.type && debug_event_pos == indel.pos && debug_event_size == indel.sequence.size())
        {
          print_log(log_severity::info,
                    __HERE__,
                    " Indel has good support in file_i,log_qual=",
                    file_i,
                    ",",
                    log_qual,
                    ",",
                    count,
                    ",",
                    anti_count_d,
                    " cov=",
                    cov);
        }
#endif // NDEBUG

        it->second.has_indel_good_support = true;
        it->second.has_realignment_support = true;
        it->second.max_log_qual = log_qual;
        it->second.max_log_qual_file_i = file_i;

        ++it;
      }
      else if (count >= 3.0 && log_qual > 0 &&
               // indel_info.sequence_reversed > 0 &&
               // indel_info.sequence_reversed < indel_info.hq_count &&
               indel_info.proper_pairs >= 1 && (indel_info.hq_count >= 5 || indel_info.max_mapq >= 25) &&
               /*indel_info.hq_count >= 2 &&*/
               indel_info.max_mapq >= 10 && indel_info.clipped < indel_info.hq_count)
      {
        // Realignment support, meaning low support but needs realignment to confirm/deny
#ifndef NDEBUG
        if (debug_event_type == indel.type && debug_event_pos == indel.pos && debug_event_size == indel.sequence.size())
        {
          print_log(log_severity::info,
                    __HERE__,
                    " Indel has realignment support in file_i,log_qual,count,acount=",
                    file_i,
                    ",",
                    log_qual,
                    ",",
                    count,
                    ",",
                    anti_count_d,
                    " cov=",
                    cov);
        }
#endif // NDEBUG

        assert(it->second.has_indel_good_support == false);
        it->second.has_realignment_support = true;
        it->second.max_log_qual = log_qual;
        it->second.max_log_qual_file_i = file_i;
        ++it;
      }
      else
      {
        // erase indel if it is not good enough
        // BOOST_LOG_TRIVIAL(info) << __HERE__ << " erasing indel=" << indel.to_string();
#ifndef NDEBUG
        if (debug_event_type == indel.type && debug_event_pos == indel.pos && debug_event_size == indel.sequence.size())
        {
          print_log(log_severity::info,
                    __HERE__,
                    " erasing indel with bad support in file_i,log_qual,count,acount=",
                    file_i,
                    ",",
                    log_qual,
                    ",",
                    count,
                    ",",
                    anti_count_d,
                    " cov=",
                    cov);
        }
#endif // NDEBUG
        it = bucket.events.erase(it);
      }
    }

    // Add coverage from bucket b
    if ((b * BUCKET_SIZE) >= REF_SIZE)
      break;

    assert(b * BUCKET_SIZE < REF_SIZE);
    long offset{b * BUCKET_SIZE};
    long const end_offset = std::min(REF_SIZE, (b + 1) * BUCKET_SIZE);

    while (offset < end_offset)
    {
      depth += static_cast<long>(cov_up[offset]) - static_cast<long>(cov_down[offset]);
      ++offset;
    }
  }

  depth = 0;

  // Find HQ SNPs and their haplotypes
  for (long b{0}; b < NUM_BUCKETS; ++b)
  {
    auto & bucket = buckets[b];

    // Check event phase
    for (auto event_it = bucket.events.begin(); event_it != bucket.events.end();) // no increment
    {
      Event const & event = event_it->first;
      EventSupport const & info = event_it->second;

      long const begin = std::max(0l, event.pos - region_begin);
      long cov{depth}; // starting coverage depth

      auto it_bool = sample_haplotypes.insert({event, Thap()});
      auto & always_together = it_bool.first->second.always_together;
      auto & ever_together = it_bool.first->second.ever_together;

      assert(begin >= b * BUCKET_SIZE);
      update_coverage(cov, begin, b, BUCKET_SIZE);
      double support_ratio = static_cast<double>(info.get_raw_support()) / static_cast<double>(cov);

      if (support_ratio < 0.3)
        support_ratio = 0.3;

      // 0 low coverage
      // 1 support
      // 2 anti support
      // 3 ambigous
      auto is_good_support =
        [&cov_down, &region_begin, &support_ratio](long local_cov,
                                                   long local_offset,
                                                   std::map<Event, EventSupport>::const_iterator event_it,
                                                   std::map<Event, EventSupport>::const_iterator event_it2,
                                                   std::map<Event, uint16_t>::const_iterator find_it,
                                                   std::map<Event, uint16_t> const & map) -> uint16_t
      {
        bool const is_indel = event_it->first.type != 'X' || event_it2->first.type != 'X';

        if (is_indel)
        {
          // at least one is an indel
          if (find_it == map.end() || find_it->second == 0)
          {
            return IS_ANY_ANTI_HAP_SUPPORT;
          }
          else
          {
            return IS_ANY_HAP_SUPPORT | IS_ANY_ANTI_HAP_SUPPORT;
          }
        }

        long const end = std::max(0l, static_cast<long>(event_it2->first.pos) - region_begin);

        // reduce coverage by reads who do not overlap the whole interval
        while (local_offset <= end)
        {
          local_cov -= static_cast<long>(cov_down[local_offset]);
          ++local_offset;
        }

        // Make sure there is some coverage
        if (local_cov <= 2)
          return 0; // low coverage

        double const support = find_it == map.end() ? 0.0 : find_it->second;

        // both are snps
        if ((support / static_cast<double>(local_cov) / support_ratio) < 0.22)
        {
          return IS_ANY_ANTI_HAP_SUPPORT;
        }
        else if ((support / static_cast<double>(local_cov) / support_ratio) > 0.78)
        {
          return IS_ANY_HAP_SUPPORT;
        }

        return IS_ANY_ANTI_HAP_SUPPORT | IS_ANY_HAP_SUPPORT;
      };

      // check this bucket
      for (auto event_it2 = std::next(event_it); event_it2 != bucket.events.end(); ++event_it2)
      {
        Event const & other_event = event_it2->first;

        if (other_event.pos == event.pos && other_event.type == event.type)
        {
          // If they share a position and type they trivially cant share haplotype
          continue;
        }

        auto find_it = info.phase.find(other_event);
        uint16_t const flags = is_good_support(cov, begin + 1, event_it, event_it2, find_it, info.phase);

        if ((flags & IS_ANY_HAP_SUPPORT) != 0u)
        {
          ever_together.insert(other_event);

          // stricter criteria on position in always_together than for ever_together
          if (other_event.pos <= (event.pos + 10))
            always_together.insert(other_event);
        }
      }

      // check next bucket
      if (b + 1 < NUM_BUCKETS)
      {
        auto const & next_bucket = buckets[b + 1];

        for (auto event_it2 = next_bucket.events.begin(); event_it2 != next_bucket.events.end(); ++event_it2)
        {
          Event const & other_event = event_it2->first;
          assert(other_event.pos > event.pos);
          assert(other_event.pos < (event.pos + 2 * BUCKET_SIZE));
          auto find_it = info.phase.find(other_event);
          uint16_t const flags = is_good_support(cov, begin + 1, event_it, event_it2, find_it, info.phase);

          if ((flags & IS_ANY_HAP_SUPPORT) != 0u)
          {
            ever_together.insert(other_event);

            // stricter criteria on position in always_together than for ever_together
            if (other_event.pos <= (event.pos + 10))
              always_together.insert(other_event);
          }
        }
      }

      // check next bucket
      if (b + 2 < NUM_BUCKETS)
      {
        auto const & next_bucket = buckets[b + 2];

        for (auto event_it2 = next_bucket.events.begin(); event_it2 != next_bucket.events.end(); ++event_it2)
        {
          Event const & other_event = event_it2->first;
          assert(other_event.pos > event.pos);

          if (other_event.pos >= (event.pos + 2 * BUCKET_SIZE))
            break;

          auto find_it = info.phase.find(other_event);
          uint16_t const flags = is_good_support(cov, begin + 1, event_it, event_it2, find_it, info.phase);

          if ((flags & IS_ANY_HAP_SUPPORT) != 0u)
          {
            ever_together.insert(other_event);
            // no need to check always_together, position difference will be too high
          }
        }
      }

      // remove SNP events, they are stored in the ph map from now on
      if (event_it->first.type == 'X')
        event_it = bucket.events.erase(event_it);
      else
        ++event_it;
    }

    // Add coverage from bucket b
    if ((b * BUCKET_SIZE) >= REF_SIZE)
      break;

    assert(b * BUCKET_SIZE < REF_SIZE);
    long offset{b * BUCKET_SIZE};
    long const end_offset = std::min(REF_SIZE, (b + 1) * BUCKET_SIZE);

    while (offset < end_offset)
    {
      depth += static_cast<long>(cov_up[offset]) - static_cast<long>(cov_down[offset]);
      ++offset;
    }
  } // for (long b{0}; b < static_cast<long>(buckets.size()); ++b)

  merge_haplotypes2(pool_haplotypes, sample_haplotypes);
}

void run_first_pass_lr(bam1_t * hts_rec,
                       HtsReader & hts_reader,
                       long file_i,
                       std::vector<BucketLR> & buckets,
                       long const BUCKET_SIZE,
                       long const region_begin,
                       std::vector<char> const & reference_sequence)
{
  long const REF_SIZE{static_cast<long>(reference_sequence.size())};
  // std::vector<uint32_t> cov_up(REF_SIZE);
  // std::vector<uint32_t> cov_down(REF_SIZE);
  // long constexpr MAX_READ_SIZE{1000000};
  Options const & copts = *(Options::const_instance());
  int min_pos{-1}; // Can be set in case of extremely high coverage

  while (true)
  {
    assert(hts_rec);

    while (hts_rec->core.n_cigar == 0 || hts_rec->core.pos < min_pos || hts_rec->core.l_qseq < 150 ||
           hts_rec->core.qual < copts.lr_mapq_filter || (hts_rec->core.flag & copts.sam_flag_filter) != 0u)
    {
      hts_rec = hts_reader.get_next_read(hts_rec);

      if (!hts_rec)
        break;
    }

    if (!hts_rec)
      break;

    std::array<char, 16> static constexpr CIGAR_MAP = {
      {'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', 'B', '*', '*', '*', '*', '*', '*'}};

    auto const & core = hts_rec->core;
    auto it = hts_rec->data;
    auto cigar_it = it + core.l_qname;
    auto seq_it = cigar_it + (core.n_cigar << 2);

    // Index of first aligned read base in read sequence
    long read_offset{0};

    // Index of first aligned read base in reference_sequence
    long ref_offset{static_cast<long>(core.pos) - region_begin};
    long const N_CIGAR = core.n_cigar;

    if (ref_offset >= REF_SIZE)
    {
      print_log(log_severity::error, __HERE__, " Unexpected ref_offset = ", ref_offset);
      std::exit(1);
    }

    Read read;
    read.name = reinterpret_cast<char *>(it);
    assert(static_cast<long>(read.name.size()) <= core.l_qname);

    // BOOST_LOG_TRIVIAL(info) << __HERE__ << " working on read=" << read.name;

    read.flags = core.flag;
    read.sequence.resize(core.l_qseq);

    for (int i{0}; i < core.l_qseq; ++i)
      read.sequence[i] = seq_nt16_str[bam_seqi(seq_it, i)];

    auto const qual_it = bam_get_qual(hts_rec);

    // uint32_t const clip_threshold = core.l_qseq < 5 ? 1 : (core.l_qseq / 5);
    // bool const is_read_clipped = is_clipped(*hts_rec, clip_threshold);

    for (long i{0}; i < N_CIGAR; ++i)
    {
      uint32_t opAndCnt;
      memcpy(&opAndCnt, cigar_it, sizeof(uint32_t));
      cigar_it += sizeof(uint32_t);

      char const cigar_operation = CIGAR_MAP[opAndCnt & 15];
      uint32_t const cigar_count = opAndCnt >> 4;
      assert(cigar_count > 0);

      switch (cigar_operation)
      {
      case 'M':
      case '=': // '=' and 'X' are typically not used by aligners (at least not by default),
      case 'X': // but we keep it here just in case
      {
        // assert((ref_offset + cigar_count - 1l) < REF_SIZE);
        // auto ref_it = reference_sequence.begin() + ref_offset;

        if ((ref_offset + static_cast<long>(cigar_count)) < 0)
        {
          read_offset += cigar_count;
          ref_offset += cigar_count;
          break;
          // no bases in this cigar operation overlap the region of interest
        }

        for (long r{0}; r < cigar_count; ++r)
        {
          long const ref_pos = ref_offset + r;

          if (ref_pos < 0)
            continue;
          else if (ref_pos >= REF_SIZE)
            break;

          char const ref = reference_sequence[ref_pos];
          long const read_pos = read_offset + r;

          if (read_pos >= static_cast<long>(read.sequence.size()))
          {
            print_log(log_severity::warning,
                      __HERE__,
                      " read_pos >= sequence.size() : ",
                      read_pos,
                      " >= ",
                      read.sequence.size());
            break;
          }

          assert(read_pos < static_cast<long>(read.sequence.size()));
          char const read_base = read.sequence[read_pos];
          char qual = *(qual_it + read_pos);
          assert(qual <= 60);

          if (qual == 0 || (ref != 'A' && ref != 'C' && ref != 'G' && ref != 'T') ||
              (read_base != 'A' && read_base != 'C' && read_base != 'G' && read_base != 'T'))
          {
            continue;
          }

          if (qual > 60)
            qual = 60;

          // transform qual to range 15-27
          long const tr_qual = 15l + std::lround((static_cast<double>(qual) * 12.0) / 60.0);

          assert(tr_qual >= 15l);
          assert(tr_qual <= 27l);
          qual = static_cast<char>(tr_qual);
          bool high_cov = add_base_to_bucket(buckets,                // buckets
                                             ref_pos + region_begin, // pos
                                             read_base,
                                             qual,
                                             region_begin,
                                             BUCKET_SIZE);

          if (high_cov)
            min_pos = ref_pos + region_begin;
        }

        read_offset += cigar_count;
        ref_offset += cigar_count;
        break;
      }

      case 'I':
      {
        if (ref_offset < 0 || (ref_offset + cigar_count >= REF_SIZE))
        {
          read_offset += cigar_count;
          break;
        }

        assert(cigar_count > 0);

        auto begin_it = read_offset < static_cast<long>(read.sequence.size()) ? read.sequence.begin() + read_offset
                                                                              : read.sequence.end();

        auto end_it = (read_offset + cigar_count) < static_cast<long>(read.sequence.size())
                      ? read.sequence.begin() + read_offset + cigar_count
                      : read.sequence.end();

        if (begin_it == end_it)
        {
          read_offset += cigar_count;
          break;
        }

        /*
        // Make sure all bases are ACGT
        if (std::all_of(begin_it,
                        end_it,
                        [](char c){
              return c == 'A' || c == 'C' || c == 'G' || c == 'T';
            }))
        {
          std::vector<char> event_sequence(begin_it, end_it);
          Event new_event = make_insertion_event(region_begin + ref_offset,
                                                 std::vector<char>(event_sequence));

          // Add to bucket events
          auto indel_event_it = add_indel_event_to_bucket(buckets,
                                                          std::move(new_event),
                                                          region_begin,
                                                          BUCKET_SIZE,
                                                          reference_sequence,
                                                          ref_offset);

#ifndef NDEBUG
          if (read.name == debug_read_name)
            print_log(log_severity::info, __HERE__, " new ins=", indel_event_it->first.to_string());
#endif // NDEBUG

          auto & event_support = indel_event_it->second;
          ++event_support.hq_count;

          if (core.qual != 255 && core.qual > event_support.max_mapq)
            event_support.max_mapq = core.qual;

          event_support.proper_pairs += ((read.flags & IS_PROPER_PAIR) != 0);
          //event_support.first_in_pairs += ((read.flags & IS_FIRST_IN_PAIR) != 0);
          event_support.sequence_reversed += ((read.flags & IS_SEQ_REVERSED) != 0);
          event_support.clipped += is_read_clipped;

          //cigar_events.push_back(indel_event_it);
        }
        */

        read_offset += cigar_count;
        break;
      }

      case 'D':
      {
        assert(cigar_count > 0);

        if (ref_offset < 0 || (ref_offset + cigar_count >= REF_SIZE))
        {
          ref_offset += cigar_count;
          break;
        }

        /*
        Event new_event =
          make_deletion_event(reference_sequence, ref_offset, region_begin + ref_offset, cigar_count);

        auto indel_event_it = add_indel_event_to_bucket(buckets,
                                                        std::move(new_event),
                                                        region_begin,
                                                        BUCKET_SIZE,
                                                        reference_sequence,
                                                        ref_offset);
#ifndef NDEBUG
        if (read.name == debug_read_name)
          print_log(log_severity::info, __HERE__, " new del=", indel_event_it->first.to_string());
#endif // NDEBUG

        auto & event_support = indel_event_it->second;
        ++event_support.hq_count;

        if (core.qual != 255 && core.qual > event_support.max_mapq)
          event_support.max_mapq = core.qual;

        event_support.proper_pairs += ((read.flags & IS_PROPER_PAIR) != 0);
        //event_support.first_in_pairs += ((read.flags & IS_FIRST_IN_PAIR) != 0);
        event_support.sequence_reversed += ((read.flags & IS_SEQ_REVERSED) != 0);
        event_support.clipped += is_read_clipped;
        */

        ref_offset += cigar_count;
        break;
      }

      case 'S':
      {
        read_offset += cigar_count;
        break;
      } // else 'H' or 'P' which move neither reference nor read
      }
    }

    // read.alignment.pos = static_cast<int32_t>(core.pos);
    // read.alignment.pos_end = region_begin + std::min(ref_offset, REF_SIZE - 1);
    //
    // assert((read.alignment.pos - region_begin) >= 0);
    // assert((read.alignment.pos - region_begin) < static_cast<long>(cov_up.size()));
    //++cov_up[read.alignment.pos - region_begin];
    //
    // assert((read.alignment.pos_end - region_begin) >= 0);
    // assert((read.alignment.pos_end - region_begin) < static_cast<long>(cov_down.size()));
    //++cov_down[read.alignment.pos_end - region_begin];

    // test if the input is sorted
    long prev_pos = hts_rec->core.pos;

    hts_rec = hts_reader.get_next_read(hts_rec);

    if (!hts_rec)
      break;

    if (hts_rec->core.pos < prev_pos)
    {
      print_log(log_severity::warning,
                __HERE__,
                " file_i=",
                file_i,
                " is not sorted. ",
                hts_rec->core.pos,
                " < ",
                prev_pos);
      assert(hts_rec->core.pos >= prev_pos);

      while (hts_rec && hts_rec->core.pos < prev_pos)
        hts_rec = hts_reader.get_next_read(hts_rec);
    }
  }

  // Check if we have too many buckets
  if ((static_cast<long>(buckets.size()) - 1l) * BUCKET_SIZE >= REF_SIZE)
  {
    long const new_size = ((REF_SIZE - 1) / BUCKET_SIZE) + 1;
    assert(new_size < static_cast<long>(buckets.size()));
    buckets.resize(new_size);
  }

  /*
  // Check if we have too many buckets
  if ((static_cast<long>(buckets.size()) - 1l) * BUCKET_SIZE >= REF_SIZE)
  {
    long const new_size = ((REF_SIZE - 1) / BUCKET_SIZE) + 1;
    assert(new_size < static_cast<long>(buckets.size()));
    buckets.resize(new_size);
  }

  long const NUM_BUCKETS{static_cast<long>(buckets.size())};
  long depth{0};

  // Find indels with good enough support to allow realignments, destroy the others
  for (long b{0}; b < NUM_BUCKETS; ++b)
  {
    auto & bucket = buckets[b];

    // Remove bad indels
    for (Tindel_events::iterator it = bucket.events.begin(); it != bucket.events.end();)  // no increment
    {
      Event const & indel = it->first;

      if (indel.type == 'X')
      {
        ++it;
        continue; // do not remove SNPs here
      }

      assert(indel.type == 'I' || indel.type == 'D');
      EventSupport const & indel_info = it->second;
      assert(!indel_info.has_realignment_support);

      long const naive_pad = 4.0 + static_cast<double>(indel.sequence.size()) / 3.0;
      long const naive_begin = std::max(0l, indel.pos - naive_pad - region_begin);
      long const naive_end = std::min(REF_SIZE, indel.pos + indel_info.span + naive_pad - region_begin);

      double const correction = indel.type == 'I' ?
                                static_cast<double>(indel.sequence.size() / 2.0 + 40.0) / 40.0 :
                                static_cast<double>(indel.sequence.size() / 3.0 + 50.0) / 50.0;

      double const count = correction * (indel_info.hq_count + indel_info.lq_count);

#ifndef NDEBUG
      if (debug_event_type == indel.type && debug_event_pos == indel.pos && debug_event_size == indel.sequence.size())
      {
        print_log(log_severity::info, __HERE__, "Indel [pos,pos+span], size [begin,end]: "
                                << indel.to_string() << " ["
                                << region_begin << " "
                                <<  indel.pos << "," << (indel.pos + indel_info.span)
                                << "] " << indel.sequence.size() << " " << correction
                                << " [" << naive_begin << "," << naive_end << "]"
                               , " count=", count, " depth=", depth);
      }
#endif // NDEBUG


      long cov{depth}; // starting coverage depth
      long offset{naive_begin};

      // Fix coverage starting before here
      if (offset <= b * BUCKET_SIZE)
      {
        while (offset < b * BUCKET_SIZE)
        {
          cov -= (static_cast<long>(cov_up[offset]) - static_cast<long>(cov_down[offset]));
          ++offset;
        }
      }
      else
      {
        offset = b * BUCKET_SIZE;

        // Add coverage until naive begin
        while (offset < naive_begin)
        {
          cov += (static_cast<long>(cov_up[offset]) - static_cast<long>(cov_down[offset]));
          ++offset;
        }
      }

      // reduce coverage by reads who do not overlap the whole naive interval
      while (offset <= naive_end)
      {
        cov -= static_cast<long>(cov_down[offset]);
        ++offset;
      }

      double const corrected_cov{std::max(static_cast<double>(cov), count)};
      double const anti_count_d{corrected_cov - count};
      uint32_t log_qual = get_log_qual_double(count, anti_count_d, 10.0);

      if (indel_info.hq_count >= 6 &&
          count >= 8.0 &&
          log_qual >= 60 &&
          indel_info.sequence_reversed > 0 &&
          indel_info.sequence_reversed < indel_info.hq_count &&
          indel_info.proper_pairs >= 3 &&
          indel_info.max_mapq >= 20 &&
          (indel_info.clipped == 0 || (indel_info.clipped + 3) <= indel_info.hq_count))
      {
#ifndef NDEBUG
        if (debug_event_type == indel.type && debug_event_pos == indel.pos && debug_event_size == indel.sequence.size())
        {
          print_log(log_severity::info, __HERE__, " Indel has good support in file_i,log_qual="
                                  << file_i << "," << log_qual << "," << count << "," << anti_count_d
                                 , " cov=", cov);
        }
#endif // NDEBUG

        it->second.has_indel_good_support = true;
        it->second.has_realignment_support = true;
        it->second.max_log_qual = log_qual;
        it->second.max_log_qual_file_i = file_i;

        ++it;
      }
      else if (count >= 3.0 &&
               log_qual > 0 &&
               //indel_info.sequence_reversed > 0 &&
               //indel_info.sequence_reversed < indel_info.hq_count &&
               indel_info.proper_pairs >= 1 &&
               indel_info.max_mapq >= 25 &&
               indel_info.clipped < indel_info.hq_count)
      {
        // Realignment support, meaning low support but needs realignment to confirm/deny
#ifndef NDEBUG
        if (debug_event_type == indel.type && debug_event_pos == indel.pos && debug_event_size == indel.sequence.size())
        {
          print_log(log_severity::info, __HERE__, " Indel has realignment support in file_i,log_qual,count,acount="
                                  << file_i << "," << log_qual << "," << count << "," << anti_count_d
                                 , " cov=", cov);
        }
#endif // NDEBUG

        assert(it->second.has_indel_good_support == false);
        it->second.has_realignment_support = true;
        it->second.max_log_qual = log_qual;
        it->second.max_log_qual_file_i = file_i;
        ++it;
      }
      else
      {
        // erase indel if it is not good enough
        //BOOST_LOG_TRIVIAL(info) << __HERE__ << " erasing indel=" << indel.to_string();
#ifndef NDEBUG
        if (debug_event_type == indel.type && debug_event_pos == indel.pos && debug_event_size == indel.sequence.size())
        {
          print_log(log_severity::info, __HERE__, " erasing indel with bad support in file_i,log_qual,count,acount="
                                  << file_i << "," << log_qual << "," << count << "," << anti_count_d
                                 , " cov=", cov);
        }
#endif // NDEBUG
        it = bucket.events.erase(it);
      }
    }

    // Add coverage from bucket b
    if ((b * BUCKET_SIZE) >= REF_SIZE)
      break;

    assert(b * BUCKET_SIZE < REF_SIZE);
    long offset{b * BUCKET_SIZE};
    long const end_offset = std::min(REF_SIZE, (b + 1) * BUCKET_SIZE);

    while (offset < end_offset)
    {
      depth += static_cast<long>(cov_up[offset]) - static_cast<long>(cov_down[offset]);
      ++offset;
    }
  }
  */
}

void realign_to_indels(std::vector<Tindel_events::iterator> const & realignment_indels,
                       std::vector<Bucket> & buckets,
                       long const max_read_size,
                       long const BUCKET_SIZE,
                       long const region_begin,
                       std::vector<char> const & reference_sequence)
{
  long const REF_SIZE = reference_sequence.size();
  long const PAD{50}; // TODO make an option
  using Tuint = uint16_t;
  paw::AlignmentOptions<Tuint> opts;
  opts.set_match(SCORE_MATCH).set_mismatch(SCORE_MISMATCH);
  opts.set_gap_open(SCORE_GAP_OPEN).set_gap_extend(SCORE_GAP_EXTEND);
  opts.left_column_free = true;
  opts.right_column_free = true;
  opts.is_clip = true;
#ifndef NDEBUG
  long _alignment_counter{0};
#endif // NDEBUG

  // Realign reads and add reference support
  for (Tindel_events::iterator indel_it : realignment_indels)
  {
    Event const & indel = indel_it->first;
    EventSupport const & indel_info = indel_it->second;
    long const indel_span = indel.pos + indel_info.span;

#ifndef NDEBUG
    if (debug_event_type == indel.type && debug_event_pos == indel.pos && debug_event_size == indel.sequence.size())
    {
      print_log(log_severity::info, __HERE__, " Realignment to indel=", indel.to_string(), " span=", indel_info.span);
    }
#endif // NDEBUG

    // create reference sequence with the indel
    long const begin_padded = std::max(0l, indel.pos - max_read_size - 2 * PAD - region_begin);
    assert(begin_padded < REF_SIZE);
    long const end_padded = indel.pos + max_read_size + 2 * PAD - region_begin;

    auto end_it = end_padded >= REF_SIZE ? reference_sequence.end() : reference_sequence.begin() + end_padded;

    std::vector<char> new_ref(reference_sequence.begin() + begin_padded, end_it);
    std::vector<int32_t> ref_pos(new_ref.size());
    std::iota(ref_pos.begin(), ref_pos.end(), 0);
    // BOOST_LOG_TRIVIAL(info) << "New ref before:\n" << std::string(new_ref.begin(), new_ref.end());

    {
      bool const is_applied = apply_indel_event(new_ref, ref_pos, indel, begin_padded + region_begin);
      assert(is_applied);

      if (!is_applied)
        continue;
    }

    // BOOST_LOG_TRIVIAL(info) << "New ref after:\n" << std::string(new_ref.begin(), new_ref.end());
    assert(buckets.size() > 0);
    long b = begin_padded / BUCKET_SIZE;
    long b_end = std::min(static_cast<long>(buckets.size()) - 1, end_padded / BUCKET_SIZE);
    assert(b < static_cast<long>(buckets.size()));

    std::vector<char> const new_ref_stable(new_ref);
    std::vector<int32_t> const ref_pos_stable(ref_pos);

    // Find the smallest bucket that might contain reads that are overlapping the indel
    while (b > 0 && buckets[b].global_max_pos_end > (indel.pos - PAD))
      --b;

    // Check any reads that might overlap the indel
    for (; b <= b_end; ++b)
    {
      assert(b >= 0);
      assert(b < static_cast<long>(buckets.size()));
      auto & bucket = buckets[b];

      if (bucket.max_pos_end <= (indel.pos - PAD))
        continue;

      for (auto & read : bucket.reads)
      {
#ifndef NDEBUG
        if (read.name == debug_read_name)
        {
          print_log(log_severity::info, __HERE__, " Found read=", debug_read_name);
        }
#endif // NDEBUG

        if (read.alignment.pos < 0 || read.sequence.size() == 0)
          continue;

        // If the variant already has the indel then we can just update the score and skip the realignment
        if (read.alignment.has_indel_event(indel_it))
          continue;

        if ((read.alignment.num_clipped_end == 0 && read.alignment.pos_end < static_cast<long>(indel.pos)) ||
            (read.alignment.pos_end + read.alignment.num_clipped_end +
               std::min(static_cast<long>(read.alignment.num_clipped_end), PAD) <
             indel.pos) ||
            (read.alignment.num_clipped_begin == 0 && read.alignment.pos > indel_span) ||
            (read.alignment.pos - read.alignment.num_clipped_begin -
               std::min(static_cast<long>(read.alignment.num_clipped_begin), PAD) >
             indel_span))
        {
          continue;
        }

        // if (read.alignment.score == static_cast<long>(read.sequence.size())) // Cannot improve score
        //  continue; // I cannot make this optimization unless I also add anti-support
        assert(new_ref == new_ref_stable);
        assert(ref_pos == ref_pos_stable);
        bool reset_new_ref{false};
        std::vector<ReadIndelEvent> applied_events{};
        applied_events.push_back({0, indel_it});

        for (auto it = read.alignment.indel_events.cbegin(); it != read.alignment.indel_events.cend(); ++it)
        {
          Event const & indel_event = it->event_it->first;
          EventSupport const & event_info = it->event_it->second;

          if (event_info.has_realignment_support)
          {
            bool ret = apply_indel_event(new_ref, ref_pos, indel_event, begin_padded + region_begin, false);

            if (ret)
            {
#ifndef NDEBUG
              if (read.name == debug_read_name)
              {
                print_log(log_severity::info, __HERE__, " Applied event ", it->event_it->first.to_string());
              }
#endif // NDEBUG

              applied_events.push_back({0, it->event_it});
              reset_new_ref = true;
            }
            else
            {
#ifndef NDEBUG
              if (read.name == debug_read_name)
              {
                apply_indel_event(new_ref, ref_pos, indel_event, begin_padded + region_begin, true);
                print_log(log_severity::info, __HERE__, " Could not apply event ", it->event_it->first.to_string());
              }
#endif // NDEBUG

              applied_events.push_back({READ_ANTI_SUPPORT, it->event_it});
            }
          }
        }

        long const old_score = read.alignment.score;

        // TODO write greedy alignment
        paw::pairwise_alignment(read.sequence, new_ref, opts);
#ifndef NDEBUG
        ++_alignment_counter;
#endif // NDEBUG
        assert(opts.get_alignment_results());
        // paw::AlignmentResults<Tuint> & ac = opts.get_alignment_cache();
        // auto const clip = opts.get_alignment_results()->apply_clipping(ac,
        //                                                               read.sequence,
        //                                                               new_ref,
        //                                                               opts.get_match(),
        //                                                               opts.get_mismatch(),
        //                                                               opts.get_gap_open(),
        //                                                               opts.get_gap_extend(),
        //                                                               opts.get_clip());

        paw::AlignmentResults const & ar = *opts.get_alignment_results();
        // auto p = ar.get_database_begin_end(read.sequence, new_ref);

        if (ar.database_begin == 0 || ar.database_end == static_cast<long>(new_ref.size()))
        {
          // auto aligned_strings = ar.get_aligned_strings(read.sequence, new_ref);
          //
          // BOOST_LOG_TRIVIAL(warning) << __HERE__ << " Insufficient padding in read alignment, skipping read. "
          //                           << read.to_string() << "\n"
          //                           << aligned_strings.first << '\n'
          //                           << aligned_strings.second;

          if (reset_new_ref)
          {
            new_ref = new_ref_stable;
            ref_pos = ref_pos_stable;
          }

          continue;
        }

#ifndef NDEBUG
        if (read.name == debug_read_name)
        {
          // if (clip.first > 0 || clip.second < static_cast<long>(read.sequence.size()))
          {
            print_log(log_severity::debug,
                      __HERE__,
                      " name=",
                      read.name,
                      " pos=",
                      read.alignment.pos,
                      " begin,end_clip=",
                      ar.clip_begin,
                      ",",
                      ar.clip_end,
                      " old_score=",
                      old_score,
                      " new_score=",
                      ar.score);
          }
        }
#endif // NDEBUG

        if (ar.score <= old_score)
        {
          // If the score is worse and the sequence is not clipped it is anti-support
          if (ar.score < old_score)
          {
            //#ifndef NDEBUG
            //            if (read.name == debug_read_name)
            //            {
            //              auto aligned_strings = ar.get_aligned_strings(read.sequence, new_ref);
            //
            //              BOOST_LOG_TRIVIAL(debug) << __HERE__ << " ANTI old_score, new_score, read = " << old_score
            //              << ", "
            //                                       << ar.score << ", "
            //                                       << read.name << '\n' << aligned_strings.first << '\n'
            //                                       << aligned_strings.second;
            //            }
            //#endif // NDEBUG

            read.alignment.add_indel_event(READ_ANTI_SUPPORT, read.flags, read.mapq, indel_it);
          }
          else if (indel_it->first.pos >= (ref_pos[ar.database_begin] + begin_padded) &&
                   indel_it->first.pos <= (ref_pos[ar.database_end] + begin_padded))
          {
            assert(ar.score == old_score);
#ifndef NDEBUG
            if (read.name == debug_read_name)
            {
              print_log(log_severity::debug, __HERE__, " SAME SCORE ", old_score, " read=", read.to_string());
            }
#endif // NDEBUG

            read.alignment.add_indel_event(READ_MULTI_SUPPORT, read.flags, read.mapq, indel_it);
          }

          if (reset_new_ref)
          {
            new_ref = new_ref_stable;
            ref_pos = ref_pos_stable;
          }

          continue;
        }

        // We got a better score
        read.alignment.replace_indel_events(read.flags, read.mapq, std::move(applied_events));

        // Update alignment
#ifndef NDEBUG
        {
          long const _next_alignment_pos = ref_pos[ar.database_begin] + region_begin + begin_padded;

          if (read.alignment.pos > (_next_alignment_pos + 300) || _next_alignment_pos > (read.alignment.pos + 300))
          {
            print_log(log_severity::warning,
                      __HERE__,
                      " change pos from=",
                      read.alignment.pos,
                      " to=",
                      _next_alignment_pos,
                      " (",
                      ref_pos[ar.database_begin],
                      " + ",
                      region_begin,
                      " + ",
                      begin_padded,
                      ")");
          }
        }
#endif // NDEBUG

        read.alignment.pos = ref_pos[ar.database_begin] + region_begin + begin_padded;
        read.alignment.pos_end = ref_pos[ar.database_end] + region_begin + begin_padded;
        read.alignment.score = ar.score;
        read.alignment.num_clipped_begin = ar.clip_begin;
        read.alignment.num_clipped_end = static_cast<long>(read.sequence.size()) - ar.clip_end;

        // count number of inserted bases at the beginning of the read
        {
          long num_ins{0};

          while (ar.database_begin + num_ins + 1 < static_cast<long>(ref_pos.size()) &&
                 ref_pos[ar.database_begin + num_ins] == ref_pos[ar.database_begin + num_ins + 1])
          {
            ++num_ins;
          }

          read.alignment.num_ins_begin = num_ins;
        }

        //#ifndef NDEBUG
        //        if (read.name == debug_read_name)
        //        {
        //          auto aligned_strings = ar.get_aligned_strings(read.sequence, new_ref);
        //          BOOST_LOG_TRIVIAL(info) << __HERE__ << " score=" << ar.score << "\n"
        //                                  << aligned_strings.first << "\n" << aligned_strings.second;
        //        }
        //#endif // NDEBUG

        if (reset_new_ref)
        {
          new_ref = new_ref_stable;
          ref_pos = ref_pos_stable;
        }
      }
    }
  }

  /*
#ifndef NDEBUG
  print_log(log_severity::info, __HERE__, " alignment counter = ", _alignment_counter);
#endif // NDEBUG
  */
  // Print counts etc
  for (Tindel_events::iterator indel_it : realignment_indels)
  {
    Event const & indel = indel_it->first;
    EventSupport & indel_info = indel_it->second;

    if (indel_info.has_indel_good_support)
      continue;

    double const correction = indel.type == 'I' ? static_cast<double>(indel.sequence.size() / 2.0 + 8.0) / 8.0
                                                : static_cast<double>(indel.sequence.size() / 3.0 + 10.0) / 10.0;

    double const count = correction * (indel_info.hq_count + indel_info.lq_count);

    bool const is_good_count = (indel_info.hq_count >= 5 && count >= 5.5) ||
                               (indel_info.span >= 5 && indel_info.hq_count >= 4 && count >= 5.0) ||
                               (indel_info.span >= 15 && indel_info.hq_count >= 3 && count >= 4.5);

#ifndef NDEBUG
    if (debug_event_type == indel.type && debug_event_pos == indel.pos && debug_event_size == indel.sequence.size())
    {
      print_log(log_severity::info,
                __HERE__,
                " Realignment results for indel=",
                indel.to_string(),
                " ",
                indel_info.hq_count,
                ",",
                indel_info.anti_count,
                " log_qual=",
                indel_info.log_qual(10),
                " info=",
                indel_info.to_string(),
                " is_good_count=",
                is_good_count,
                " is_good_info=",
                indel_info.is_good_indel());
    }
#endif // NDEBUG

    // Check if it is good after realignment
    if (is_good_count && indel_info.is_good_indel())
    {
      indel_info.has_indel_good_support = true;
    }
#ifndef NDEBUG
    else
    {
      assert(!indel_info.has_indel_good_support);
    }
#endif // NDEBUG
  }
}

void read_hts_and_return_realignment_indels(bam1_t * hts_rec,
                                            HtsReader & hts_reader,
                                            std::vector<Bucket> & buckets,
                                            long & max_read_size,
                                            long const BUCKET_SIZE,
                                            long const region_begin,
                                            std::vector<char> const & reference_sequence)
{
  long const REF_SIZE{static_cast<long>(reference_sequence.size())};
  int32_t global_max_pos_end{0};

  while (true)
  {
    assert(hts_rec);

    while (hts_rec->core.n_cigar == 0 || hts_rec->core.pos < region_begin)
    {
      hts_rec = hts_reader.get_next_read(hts_rec);

      if (!hts_rec)
        break;
    }

    if (!hts_rec)
      break;

    std::array<char, 16> static constexpr CIGAR_MAP = {
      {'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', 'B', '*', '*', '*', '*', '*', '*'}};

    auto const & core = hts_rec->core;
    auto it = hts_rec->data;
    auto cigar_it = it + core.l_qname;
    auto seq_it = cigar_it + (core.n_cigar << 2);

    // Index of first aligned read base in read sequence
    long read_offset{0};

    // Index of first aligned read base in reference_sequence
    long ref_offset{static_cast<long>(core.pos) - region_begin};

    // index in buckets vector
    long const bucket_index{ref_offset / BUCKET_SIZE};

    // take note of the maximum read size
    if (core.l_qseq > max_read_size)
      max_read_size = core.l_qseq;

    assert(bucket_index < static_cast<long>(buckets.size()));

    long const N_CIGAR = core.n_cigar;

    if (ref_offset < 0 || ref_offset >= REF_SIZE)
    {
      print_log(log_severity::error, __HERE__, " Unexpected ref_offset = ", ref_offset);
      std::exit(1);
    }

    Read read;
    read.name = reinterpret_cast<char *>(it);

    if ((core.flag & IS_FIRST_IN_PAIR) != 0u)
      read.name.append("/1");
    else
      read.name.append("/2");

#ifndef NDEBUG
    if (read.name == debug_read_name)
    {
      print_log(log_severity::info, __HERE__, " Found debug read=", debug_read_name);
    }
#endif // NDEBUG

    read.mate_pos = static_cast<int32_t>(core.mpos);
    read.flags = core.flag;
    read.mapq = core.qual;

    read.sequence.resize(core.l_qseq);

    for (int i = 0; i < core.l_qseq; ++i)
      read.sequence[i] = seq_nt16_str[bam_seqi(seq_it, i)];

    read.qual.resize(core.l_qseq);

    {
      auto qual_it = bam_get_qual(hts_rec);

      for (int i = 0; i < core.l_qseq; ++qual_it, ++i)
        read.qual[i] = static_cast<char>(*qual_it);
    }

    read.alignment.score = 0;

    for (long i{0}; i < N_CIGAR; ++i)
    {
      if (ref_offset >= REF_SIZE)
      {
#ifndef NDEBUG
        print_log(log_severity::warning,
                  __HERE__,
                  " While processing read=",
                  reinterpret_cast<char *>(hts_rec->data),
                  " went beyond the region if interest");
#endif // NDEBUG
        break;
      }

      uint32_t opAndCnt;
      memcpy(&opAndCnt, cigar_it, sizeof(uint32_t));
      cigar_it += sizeof(uint32_t);

      char const cigar_operation = CIGAR_MAP[opAndCnt & 15];
      uint32_t const count = opAndCnt >> 4;
      assert(count > 0);

      switch (cigar_operation)
      {
      case 'M':
      case '=': // '=' and 'X' are typically not used by aligners (at least not by default),
      case 'X': // but we keep it here just in case
      {
        assert(ref_offset < REF_SIZE);
        auto ref_it = reference_sequence.begin() + ref_offset;

        for (long r{0}; r < count; ++r, ++ref_it)
        {
          char const alt = read.sequence[r + read_offset];

          if (alt != *ref_it && alt != 'N' && *ref_it != 'N')
          {
            read.alignment.score -= SCORE_MISMATCH;
          }
          else
          {
            read.alignment.score += SCORE_MATCH;
          }
        }

        read_offset += count;
        ref_offset += count;
        break;
      }

      case 'I':
      {
        assert(count > 0);

        auto begin_it = read_offset < static_cast<long>(read.sequence.size()) ? read.sequence.begin() + read_offset
                                                                              : read.sequence.end();

        auto end_it = (read_offset + count) < static_cast<long>(read.sequence.size())
                      ? read.sequence.begin() + read_offset + count
                      : read.sequence.end();

        if (begin_it == end_it)
        {
          break;
        }

        std::vector<char> event_sequence(begin_it, end_it);

        Event new_event = make_insertion_event(region_begin + ref_offset, std::move(event_sequence));

        auto indel_event_it = add_indel_event_to_bucket(buckets,
                                                        std::move(new_event),
                                                        region_begin,
                                                        BUCKET_SIZE,
                                                        reference_sequence,
                                                        ref_offset);

        if (!indel_event_it->second.has_realignment_support)
          read.alignment.score -= (SCORE_GAP_OPEN + (count - 1) * SCORE_GAP_EXTEND);
        else
          read.alignment.score += SCORE_MATCH * count;

        read.alignment.add_indel_event(read_offset, read.flags, read.mapq, indel_event_it);
        read_offset += count;
        break;
      }

      case 'D':
      {
        assert(count > 0);

        if (ref_offset + count >= REF_SIZE)
        {
#ifndef NDEBUG
          print_log(log_severity::warning,
                    __HERE__,
                    " While processing read=",
                    reinterpret_cast<char *>(hts_rec->data),
                    " went beyond the region if interest");
#endif // NDEBUG
          break;
        }

        Event new_event = make_deletion_event(reference_sequence, ref_offset, region_begin + ref_offset, count);

        auto indel_event_it = add_indel_event_to_bucket(buckets,
                                                        std::move(new_event),
                                                        region_begin,
                                                        BUCKET_SIZE,
                                                        reference_sequence,
                                                        ref_offset);

        if (!indel_event_it->second.has_realignment_support)
          read.alignment.score -= (SCORE_GAP_OPEN + (count - 1) * SCORE_GAP_EXTEND);

        read.alignment.add_indel_event(read_offset, read.flags, read.mapq, indel_event_it);
        ref_offset += count;
        break;
      }

      case 'S':
      {
        read_offset += count;
        set_bit(read.flags, IS_CLIPPED);
        read.alignment.score -= SCORE_CLIP;

        if (i == 0)
        {
          read.alignment.num_clipped_begin = count;
        }
        else
        {
          // if ((i + 1 != N_CIGAR) && (i + 2 != N_CIGAR))
          //{
          //  // write cigar
          //  std::ostringstream ss;
          //  auto cigar_it2 = hts_rec->data + hts_rec->core.l_qname;
          //
          //  for (long j{0}; j < static_cast<long>(hts_rec->core.n_cigar); ++j)
          //  {
          //    uint32_t opAndCnt;
          //    memcpy(&opAndCnt, cigar_it2, sizeof(uint32_t));
          //    cigar_it2 += sizeof(uint32_t);
          //
          //    char const cigar_operation = CIGAR_MAP[opAndCnt & 15];
          //    uint32_t const count = opAndCnt >> 4;
          //    ss << count << cigar_operation;
          //  }
          //
          //  BOOST_LOG_TRIVIAL(warning) << __HERE__ << " Unexpected cigar: " << ss.str() << " in read "
          //                             << reinterpret_cast<char *>(hts_rec->data) << " "
          //                             << hts_rec->core.n_cigar << " "
          //                             << N_CIGAR;
          //}

          assert((i + 1 == N_CIGAR) || (i + 2 == N_CIGAR));
          read.alignment.num_clipped_end = count;
        }

        break;
      } // else 'H' or 'P' which move neither reference nor read
      }
    }

    read.alignment.pos = static_cast<int32_t>(core.pos);
    read.alignment.pos_end = region_begin + ref_offset;
    auto & bucket = buckets[bucket_index];
    int32_t const end_with_clip = read.alignment.pos_end + read.alignment.num_clipped_end;

    if (end_with_clip > bucket.max_pos_end)
    {
      bucket.max_pos_end = end_with_clip;
      global_max_pos_end = std::max(global_max_pos_end, end_with_clip);
    }

    bucket.global_max_pos_end = global_max_pos_end;
    bucket.reads.push_back(std::move(read));

    hts_rec = hts_reader.get_next_read(hts_rec);

    if (!hts_rec)
      break;
  }

  // realign
}

void parallel_first_pass(std::vector<std::string> * hts_paths_ptr,
                         std::map<Event, Thap> * pool_haplotypes_ptr,
                         std::vector<std::vector<BucketFirstPass>> * file_buckets_first_pass_ptr,
                         std::vector<std::string> * sample_names_ptr,
                         std::vector<char> * reference_sequence_ptr,
                         long const BUCKET_SIZE,
                         long const region_begin,
                         long const lowest_file_i)
{
  assert(hts_paths_ptr);
  assert(pool_haplotypes_ptr);

  std::vector<std::string> const & hts_paths = *hts_paths_ptr;
  std::map<Event, Thap> & pool_haplotypes = *pool_haplotypes_ptr;
  std::vector<char> const & reference_sequence = *reference_sequence_ptr;

  for (long i{0}; i < static_cast<long>(hts_paths.size()); ++i)
  {
    // BAM/CRAM record storage
    HtsStore store;

    // Open BAM/CRAM file
    HtsReader hts_reader(store);
    hts_reader.open(hts_paths[i], std::string("."), "");

    // Add sample
    if (hts_reader.samples.size() > 1)
    {
      print_log(log_severity::error,
                __HERE__,
                " We found file with multiple samples, sorry, ",
                "this is currently not supported.");
      std::exit(1);
    }

    assert(hts_reader.samples.size() == 1);
    (*sample_names_ptr)[lowest_file_i + i] = hts_reader.samples[0];

    // Get the first read
    bam1_t * hts_rec = hts_reader.get_next_read();

    // Check if there were any reads
    if (!hts_rec)
    {
      hts_reader.close();
      continue;
    }

    run_first_pass(hts_rec,
                   hts_reader,
                   lowest_file_i + i,
                   (*file_buckets_first_pass_ptr)[lowest_file_i + i],
                   pool_haplotypes,
                   BUCKET_SIZE,
                   region_begin,
                   reference_sequence);

    // Close BAM/CRAM file
    hts_reader.close();
  } // for (long file_i{0}; file_i < static_cast<long>(hts_paths.size()); ++file_i)
}

void parallel_first_pass_lr(std::vector<std::string> * hts_paths_ptr,
                            std::vector<std::vector<BucketLR>> * file_buckets_first_pass_ptr,
                            std::vector<std::string> * sample_names_ptr,
                            std::vector<char> * reference_sequence_ptr,
                            std::string const * region_str_ptr,
                            long const BUCKET_SIZE,
                            long const region_begin,
                            long const lowest_file_i)
{
  assert(hts_paths_ptr);

  std::vector<std::string> const & hts_paths = *hts_paths_ptr;
  std::vector<char> const & reference_sequence = *reference_sequence_ptr;
  std::string const & region_str = *region_str_ptr;

  for (long i{0}; i < static_cast<long>(hts_paths.size()); ++i)
  {
    // BAM/CRAM record storage
    HtsStore store;

    // Open BAM/CRAM file
    HtsReader hts_reader(store);
    hts_reader.open(hts_paths[i], region_str, "");

    // Add sample
    if (hts_reader.samples.size() > 1)
    {
      print_log(log_severity::error,
                __HERE__,
                " We found file with multiple samples, sorry, ",
                "this is currently not supported.");
      std::exit(1);
    }

    assert(hts_reader.samples.size() == 1);
    (*sample_names_ptr)[lowest_file_i + i] = hts_reader.samples[0];

    // Get the first read
    bam1_t * hts_rec = hts_reader.get_next_read();

    // Check if there were any reads
    if (!hts_rec)
    {
      hts_reader.close();
      continue;
    }

    run_first_pass_lr(hts_rec,
                      hts_reader,
                      lowest_file_i + i,
                      (*file_buckets_first_pass_ptr)[lowest_file_i + i],
                      BUCKET_SIZE,
                      region_begin,
                      reference_sequence);

    // Close BAM/CRAM file
    hts_reader.close();
  } // for (long file_i{0}; file_i < static_cast<long>(hts_paths.size()); ++file_i)
}

void parallel_second_pass(std::string const * hts_path_ptr,
                          std::vector<char> const * reference_sequence_ptr,
                          Tindel_events * indel_events_ptr,
                          std::vector<Tindel_events::iterator> * indel_to_realign_ptr,
                          long const BUCKET_SIZE,
                          long const NUM_BUCKETS,
                          long const region_begin)
{
  if (indel_to_realign_ptr->size() == 0)
    return;

  std::string const & hts_path = *hts_path_ptr;
  std::vector<char> const & reference_sequence = *reference_sequence_ptr;
  Tindel_events const & indel_events = *indel_events_ptr;
  std::vector<Tindel_events::iterator> & indel_to_realign = *indel_to_realign_ptr;

  std::vector<Bucket> buckets;
  buckets.resize(NUM_BUCKETS);
  long max_read_size{100};

  // I/O
  {
    // BAM/CRAM record storage
    HtsStore store;

    // Open BAM/CRAM file
    HtsReader hts_reader(store);
    hts_reader.open(hts_path, std::string("."), "");

    // Get the first read
    bam1_t * hts_rec = hts_reader.get_next_read();

    // Check if there were any reads
    if (!hts_rec)
    {
      hts_reader.close();
      return;
    }

    read_hts_and_return_realignment_indels(hts_rec,
                                           hts_reader,
                                           buckets,
                                           max_read_size,
                                           BUCKET_SIZE,
                                           region_begin,
                                           reference_sequence);

    // Close BAM/CRAM file
    hts_reader.close();
  } // I/O ends

  std::vector<Tindel_events::iterator> nearby_good_events;
  long constexpr NEARBY_BP{60};

  for (Tindel_events::iterator indel : indel_to_realign)
  {
    Tindel_events::iterator indel_it = indel;

    // Check before the indel
    while (indel_it != indel_events.begin())
    {
      --indel_it;

      if (indel_it->second.has_indel_good_support && indel_it->first.pos + NEARBY_BP >= indel->first.pos)
      {
        // It also needs to be found in the file
        if (is_indel_in_bucket(buckets, indel_it->first, region_begin, BUCKET_SIZE))
          nearby_good_events.push_back(indel_it);
      }
    }

    // Check after the indel
    indel_it = std::next(indel);

    while (indel_it != indel_events.end())
    {
      if (indel_it->second.has_indel_good_support && indel_it->first.pos - NEARBY_BP <= indel->first.pos)
      {
        // It also needs to be found in the file
        if (is_indel_in_bucket(buckets, indel_it->first, region_begin, BUCKET_SIZE))
          nearby_good_events.push_back(indel_it);
      }

      ++indel_it;
    }
  }

  std::sort(nearby_good_events.begin(),
            nearby_good_events.end(),
            [](Tindel_events::iterator const & a, Tindel_events::iterator const & b) -> bool
            { return a->first.pos < b->first.pos; });

  auto end_unique_it = std::unique(nearby_good_events.begin(), nearby_good_events.end());
  std::move(nearby_good_events.begin(), end_unique_it, std::back_inserter(indel_to_realign));

  std::sort(indel_to_realign.begin(),
            indel_to_realign.end(),
            [](Tindel_events::iterator const & a, Tindel_events::iterator const & b) -> bool
            {
              return a->second.has_indel_good_support > b->second.has_indel_good_support ||
                     (a->second.has_indel_good_support == b->second.has_indel_good_support &&
                      a->first.pos < b->first.pos);
            });

#ifndef NDEBUG
  print_log(log_severity::debug, __HERE__, " Nearby realignment indels=", nearby_good_events.size());
#endif // NDEBUG

  realign_to_indels(indel_to_realign, // realignment_indels,
                    buckets,
                    max_read_size,
                    BUCKET_SIZE,
                    region_begin,
                    reference_sequence);

#ifndef NDEBUG
  print_log(log_severity::debug, __HERE__, " Realignment DONE.");
#endif // NDEBUG
}

void streamlined_discovery(std::vector<std::string> const & hts_paths,
                           std::string const & ref_path,
                           std::string const & region_str,
                           gyper::Vcf & vcf)
{
  long const NUM_FILES = hts_paths.size();
  GenomicRegion genomic_region(region_str);
  uint64_t const chromosome_offset = genomic_region.get_absolute_position(1);
  long const region_begin{genomic_region.begin};
  std::vector<char> reference_sequence;
  open_and_read_reference_genome(reference_sequence, ref_path, genomic_region);
  long constexpr BUCKET_SIZE{50}; // TODO make an option
  std::vector<std::string> sample_names;
  sample_names.resize(NUM_FILES);

  std::vector<std::unique_ptr<std::vector<std::string>>> spl_hts_paths{};
  std::vector<std::unique_ptr<std::map<Event, Thap>>> spl_snp_hq_haplotypes{};
  long jobs{1};

  {
    long num_parts{1};
    _determine_num_jobs_and_num_parts(jobs, num_parts, NUM_FILES);

    {
      auto it = hts_paths.cbegin();
      spl_hts_paths.reserve(num_parts);
      spl_snp_hq_haplotypes.reserve(num_parts);

      for (long i = 0; i < num_parts; ++i)
      {
        auto end_it = it + NUM_FILES / num_parts + (i < (NUM_FILES % num_parts));
        assert(std::distance(hts_paths.cbegin(), end_it) <= NUM_FILES);
        spl_hts_paths.emplace_back(new std::vector<std::string>(it, end_it));
        spl_snp_hq_haplotypes.emplace_back(new std::map<Event, Thap>());
        it = end_it;
      }
    }
  }

#ifndef NDEBUG
  print_log(log_severity::info, __HERE__, " First discovery pass starting.");
#endif // NDEBUG

  long const NUM_POOLS{static_cast<long>(spl_hts_paths.size())};

#ifndef NDEBUG
  print_log(log_severity::debug, __HERE__, " Number of pools = ", NUM_POOLS);
#endif // NDEBUG

  Tindel_events indel_events;
  long NUM_BUCKETS{0};
  std::map<Event, Thap> haplotypes;

  // FIRST PASS
  {
    std::vector<std::vector<BucketFirstPass>> file_buckets_first_pass(NUM_FILES);

    // Parallel first pass
    {
      paw::Station first_station(jobs);
      long file_i{0};

      for (long i{0}; i < NUM_POOLS - 1; ++i)
      {
        first_station.add_work(parallel_first_pass,
                               &(*spl_hts_paths[i]),
                               &(*spl_snp_hq_haplotypes[i]),
                               &file_buckets_first_pass,
                               &sample_names,
                               &reference_sequence,
                               BUCKET_SIZE,
                               region_begin,
                               file_i);

        file_i += spl_hts_paths[i]->size();
      }

      // Do the last pool on the current thread
      first_station.add_to_thread(jobs - 1,
                                  parallel_first_pass,
                                  &(*spl_hts_paths[NUM_POOLS - 1]),
                                  &(*spl_snp_hq_haplotypes[NUM_POOLS - 1]),
                                  &file_buckets_first_pass,
                                  &sample_names,
                                  &reference_sequence,
                                  BUCKET_SIZE,
                                  region_begin,
                                  file_i);

      first_station.join();
    }

    // Find SNPs and their SNP haplotypes
    haplotypes = std::move(*spl_snp_hq_haplotypes[0]);

    for (long i{1}; i < NUM_POOLS; ++i)
    {
      std::map<Event, Thap> & pool_haplotypes = *spl_snp_hq_haplotypes[i];
      merge_haplotypes2(haplotypes, pool_haplotypes);
    }

    // Create buckets for second pass
    for (auto && buckets_first : file_buckets_first_pass)
    {
      if (static_cast<long>(buckets_first.size()) > NUM_BUCKETS)
        NUM_BUCKETS = buckets_first.size();

      for (long b{0}; b < static_cast<long>(buckets_first.size()); ++b)
      {
        auto && bucket_first = buckets_first[b];

        for (auto && indel_event : bucket_first.events)
        {
          // Only add indels
          assert(indel_event.first.type != 'X');

          if (indel_event.first.type == 'X')
            continue;

#ifndef NDEBUG
          if (debug_event_type == indel_event.first.type && debug_event_pos == indel_event.first.pos &&
              debug_event_size == indel_event.first.sequence.size())
          {
            print_log(log_severity::info,
                      __HERE__,
                      " Indel=",
                      indel_event.first.to_string(),
                      " good_support=",
                      indel_event.second.has_indel_good_support);
          }
#endif

          assert(indel_event.second.has_realignment_support);
          assert(indel_event.first.type == 'D' || indel_event.first.type == 'I');
          auto insert_it = indel_events.insert(indel_event);

          if (!insert_it.second)
          {
            // nothing was inserted
            EventSupport & old_info = insert_it.first->second;
            old_info.has_indel_good_support |= indel_event.second.has_indel_good_support;

            // Found - Check if support is better
            if (indel_event.second.max_log_qual > old_info.max_log_qual)
            {
              old_info.max_log_qual = indel_event.second.max_log_qual;
              old_info.max_log_qual_file_i = indel_event.second.max_log_qual_file_i;
            }
          }
        }
      }
    }
  }

  // FIRST PASS ENDS
#ifndef NDEBUG
  print_log(log_severity::info, __HERE__, " First pass DONE. Starting second pass.");
#endif // NDEBUG

  // SECOND PASS BEGINS
  {
    std::vector<std::vector<Tindel_events::iterator>> indel_to_realign;
    indel_to_realign.resize(NUM_FILES);

    // Add indels to events
    for (auto it = indel_events.begin(); it != indel_events.end(); ++it)
    {
      it->second.clear();

      assert(it->second.anti_count == 0);
      assert(it->second.multi_count == 0);
      assert(it->second.has_realignment_support);

      if (!it->second.has_indel_good_support)
      {
#ifndef NDEBUG
        if (debug_event_type == it->first.type && debug_event_pos == it->first.pos &&
            debug_event_size == it->first.sequence.size())
        {
          print_log(log_severity::info,
                    __HERE__,
                    " realignment indel=",
                    it->first.to_string(),
                    " qual,file_i=",
                    it->second.max_log_qual,
                    " ",
                    it->second.max_log_qual_file_i);
        }
#endif // NDEBUG

        assert(it->second.max_log_qual_file_i < static_cast<long>(indel_to_realign.size()));
        indel_to_realign[it->second.max_log_qual_file_i].push_back(it);
      }
    }

    {
      paw::Station second_station(jobs);

      for (long file_i{0}; file_i < (NUM_FILES - 1); ++file_i)
      {
        if (indel_to_realign[file_i].size() > 0)
        {
          second_station.add_work(parallel_second_pass,
                                  &(hts_paths[file_i]),
                                  &reference_sequence,
                                  &indel_events,
                                  &(indel_to_realign[file_i]),
                                  BUCKET_SIZE,
                                  NUM_BUCKETS,
                                  region_begin);
        }
      }

      // Last file is run on current thread
      long const file_i{(NUM_FILES - 1)};

      if (indel_to_realign[file_i].size() > 0)
      {
        second_station.add_to_thread(jobs - 1,
                                     parallel_second_pass,
                                     &(hts_paths[file_i]),
                                     &reference_sequence,
                                     &indel_events,
                                     &(indel_to_realign[file_i]),
                                     BUCKET_SIZE,
                                     NUM_BUCKETS,
                                     region_begin);
      }

      second_station.join();
    }

    long event_index{1}; // tracks the index of the event in the haplotypes map

    for (auto event_it = haplotypes.begin(); event_it != haplotypes.end(); ++event_it, ++event_index)
    {
      Event const & event = event_it->first;

      // if the event is an indel, check if it was removed in the second iteration
      if (event.type != 'X')
      {
        assert(event.type == 'I' || event.type == 'D');
        auto find_it = indel_events.find(event);
        assert(find_it != indel_events.end());

        if (find_it == indel_events.end() || !find_it->second.has_indel_good_support)
          continue; // skip
      }

      uint32_t const abs_pos = event.pos + chromosome_offset;

      Variant variant{};
      variant.abs_pos = abs_pos;
      assert(event.pos >= region_begin);
      assert((event.pos - region_begin) < static_cast<long>(reference_sequence.size()));

      if (event.type == 'X')
      {
        variant.seqs.push_back({reference_sequence[event.pos - region_begin]});
        variant.seqs.push_back(event.sequence);
        variant.type = 'X';
      }
      else if (event.type == 'I')
      {
        variant.seqs.push_back({});
        variant.seqs.push_back({event.sequence});
        variant.add_base_in_front(true);
        variant.type = 'I';
      }
      else if (event.type == 'D')
      {
        variant.seqs.push_back({event.sequence});
        variant.seqs.push_back({});
        variant.add_base_in_front(true);
        variant.type = 'D';
      }
      else
      {
        assert(false);
      }

      {
        std::ostringstream ss_hap;
        std::ostringstream ss_anti;
        bool is_hap_stream_empty{true};
        bool is_anti_stream_empty{true};
        long next_event_index = event_index + 1;

        // Check if the next event within BUCKET_SIZE bp are in phased with this one ever
        for (auto event_it2 = std::next(event_it); event_it2 != haplotypes.end(); ++event_it2, ++next_event_index)
        {
          Event const & next_event = event_it2->first;

          if (next_event.pos >= event.pos + 2 * BUCKET_SIZE)
            break;

          // skip if the event is a bad indel
          if (next_event.type != 'X')
          {
            assert(next_event.type == 'I' || next_event.type == 'D');
            auto find_it = indel_events.find(next_event);
            assert(find_it != indel_events.end());

            if (find_it == indel_events.end() || !find_it->second.has_indel_good_support)
              continue; // skip
          }

          if (event_it->second.always_together.count(next_event) > 0)
          {
            // these two were always seen together
            if (is_hap_stream_empty)
              is_hap_stream_empty = false;
            else
              ss_hap << ',';

            ss_hap << next_event_index;
          }
          else if (event_it->second.ever_together.count(next_event) == 0)
          {
            // these two events were never seen together - anti haplotype support
            if (is_anti_stream_empty)
              is_anti_stream_empty = false;
            else
              ss_anti << ',';

            ss_anti << next_event_index;
          }
        }

        variant.infos["GT_ID"] = std::to_string(event_index);

        if (!is_hap_stream_empty)
          variant.infos["GT_HAPLOTYPE"] = ss_hap.str();

        if (!is_anti_stream_empty)
          variant.infos["GT_ANTI_HAPLOTYPE"] = ss_anti.str();
      }

      vcf.variants.push_back(std::move(variant));
    }
  }
}

void streamlined_lr_genotyping(std::vector<std::string> const & hts_paths,
                               std::string const & ref_path,
                               std::string const & region_str,
                               gyper::Vcf & vcf)
{
  print_log(log_severity::info, "Start lr WIP on region ", region_str);

  long const NUM_FILES = hts_paths.size();
  GenomicRegion genomic_region(region_str);
  uint64_t const chromosome_offset = genomic_region.get_absolute_position(1);
  long const region_begin{genomic_region.begin};
  std::vector<char> reference_sequence;
  open_and_read_reference_genome(reference_sequence, ref_path, genomic_region);
  long const REF_SIZE = reference_sequence.size();
  long constexpr BUCKET_SIZE{50}; // TODO make an option
  std::vector<std::string> sample_names;
  sample_names.resize(NUM_FILES);

  std::vector<std::unique_ptr<std::vector<std::string>>> spl_hts_paths{};
  long jobs{1};

  {
    long num_parts{1};
    _determine_num_jobs_and_num_parts(jobs, num_parts, NUM_FILES);

    {
      auto it = hts_paths.cbegin();
      spl_hts_paths.reserve(num_parts);

      for (long i = 0; i < num_parts; ++i)
      {
        auto end_it = it + NUM_FILES / num_parts + (i < (NUM_FILES % num_parts));
        assert(std::distance(hts_paths.cbegin(), end_it) <= NUM_FILES);
        spl_hts_paths.emplace_back(new std::vector<std::string>(it, end_it));
        it = end_it;
      }
    }
  }

#ifndef NDEBUG
  print_log(log_severity::info, __HERE__, " First pass starting.");
#endif // NDEBUG

  long const NUM_POOLS{static_cast<long>(spl_hts_paths.size())};

#ifndef NDEBUG
  print_log(log_severity::debug, __HERE__, " Number of pools = ", NUM_POOLS);
#endif // NDEBUG

  // Tindel_events indel_events;
  // Tindel_events events;
  long num_buckets{0};
  std::set<SnpEvent> snp_events;

  // FIRST PASS
  std::vector<std::vector<BucketLR>> file_buckets_first_pass(NUM_FILES);

  // Parallel first pass
  {
    paw::Station first_station(jobs);
    long file_i{0};

    for (long i{0}; i < NUM_POOLS - 1; ++i)
    {
      first_station.add_work(parallel_first_pass_lr,
                             &(*spl_hts_paths[i]),
                             &file_buckets_first_pass,
                             &sample_names,
                             &reference_sequence,
                             &region_str,
                             BUCKET_SIZE,
                             region_begin,
                             file_i);

      file_i += spl_hts_paths[i]->size();
    }

    // Do the last pool on the current thread
    first_station.add_to_thread(jobs - 1,
                                parallel_first_pass_lr,
                                &(*spl_hts_paths[NUM_POOLS - 1]),
                                &file_buckets_first_pass,
                                &sample_names,
                                &reference_sequence,
                                &region_str,
                                BUCKET_SIZE,
                                region_begin,
                                file_i);

    first_station.join();
  }

  // Merge files that have the same sample
  std::vector<std::string> new_sample_names;
  std::vector<std::vector<BucketLR>> new_file_buckets;
  std::unordered_map<std::string, long> seen_sample_names;

  for (long s{0}; s < NUM_FILES; ++s)
  {
    assert(s < static_cast<long>(sample_names.size()));

    std::string const & sample_name = sample_names[s];
    assert(new_sample_names.size() == new_file_buckets.size());
    auto insert_it = seen_sample_names.insert({sample_name, (new_sample_names.size() - 1)});

    if (insert_it.second)
    {
      // new sample name
      new_file_buckets.push_back(std::move(file_buckets_first_pass[s]));
      new_sample_names.push_back(sample_name);
    }
    else
    {
      // merge with old bucket
      long const old_s = insert_it.first->second;

      assert(s < static_cast<long>(sample_names.size()));
      assert(s < static_cast<long>(file_buckets_first_pass.size()));
      assert(old_s < static_cast<long>(new_sample_names.size()));
      assert(old_s < static_cast<long>(new_file_buckets.size()));
      assert(old_s < s);

      std::vector<BucketLR> & old_buckets = file_buckets_first_pass[old_s];
      std::vector<BucketLR> const & new_buckets = file_buckets_first_pass[s];

      merge_bucket_lr(old_buckets, new_buckets);
    }
  }

  assert(new_sample_names.size() == new_file_buckets.size());
  file_buckets_first_pass.clear();
  sample_names.clear();

  // Get SNP events
  for (auto && buckets : new_file_buckets)
  {
    if (static_cast<long>(buckets.size()) > num_buckets)
      num_buckets = buckets.size();

    // Find candidate SNPs in the pileups
    for (long b{0}; b < static_cast<long>(buckets.size()); ++b)
    {
      auto & bucket = buckets[b];

      if (bucket.pileup.size() == 0)
        continue;

      long const pos = b * BUCKET_SIZE;
      std::vector<BaseCount> const & pileup = bucket.pileup;

      for (long p{0}; p < static_cast<long>(pileup.size()); ++p)
      {
        assert(pos >= 0);

        if (pos + p >= REF_SIZE)
          break;

        char const ref_base = reference_sequence[pos + p];
        long const ref_index = base2index(ref_base);

        if (ref_index < 0)
          continue;

        BaseCount const & base_count = pileup[p];

        std::array<int32_t, 4> const & bc = base_count.acgt;
        std::array<int64_t, 4> const & qs = base_count.acgt_qualsum;
        std::array<long, 4> idx{{0, 1, 2, 3}};

        std::sort(idx.begin(), idx.end(), [&qs](int64_t i1, int64_t i2) { return qs[i1] < qs[i2]; });

        long const first = idx[3];
        long const second = idx[2];
        long const third = idx[1];

        if (first != ref_index && bc[first] >= 3 &&
            (((qs[first] - qs[second]) >= 30) || ((qs[first] - qs[third]) >= 50)))
        {
          SnpEvent snp_event(region_begin + pos + p, index2base[first]);
          snp_events.insert(snp_event);
        }

        if (second != ref_index && bc[second] >= 4 && (qs[second] - qs[third]) >= 50 &&
            ((static_cast<double>(qs[second]) / static_cast<double>(qs[0] + qs[1] + qs[2] + qs[3])) > 0.3))
        {
          SnpEvent snp_event(region_begin + pos + p, index2base[second]);
          snp_events.insert(snp_event);
        }
      }
    }
  }

  print_log(log_severity::info, __HERE__, " Number of events: ", snp_events.size());

  for (auto snp_event_it = snp_events.begin(); snp_event_it != snp_events.end(); ++snp_event_it)
  {
    SnpEvent const & snp_event = *snp_event_it;
    uint32_t const abs_pos = snp_event.pos + chromosome_offset;

    Variant variant{};
    variant.abs_pos = abs_pos;

    assert(snp_event.pos >= region_begin);
    assert((snp_event.pos - region_begin) < static_cast<long>(reference_sequence.size()));
    char const ref_base = reference_sequence[snp_event.pos - region_begin];
    variant.seqs.push_back({ref_base});
    variant.seqs.push_back({snp_event.base});

    // also check next SNPs
    {
      auto next_snp_event_it = std::next(snp_event_it);

      while (next_snp_event_it != snp_events.end())
      {
        if (next_snp_event_it->pos != snp_event_it->pos)
          break;

        variant.seqs.push_back({next_snp_event_it->base});
        snp_event_it = next_snp_event_it;
        std::advance(next_snp_event_it, 1);
      }
    }

    assert(reference_sequence[snp_event.pos - region_begin] != snp_event.base);
    variant.type = 'X';

    long const b = (snp_event.pos - region_begin) / BUCKET_SIZE;
    long const p = (snp_event.pos - region_begin) % BUCKET_SIZE;
    long const cnum = variant.seqs.size();

    // gather genotype calls
    for (auto const & buckets : new_file_buckets)
    {
      SampleCall new_call;
      long const pl_len = (cnum * (cnum + 1)) / 2;
      new_call.phred.resize(pl_len);
      new_call.coverage.resize(cnum);

      if (b >= static_cast<long>(buckets.size()))
      {
        variant.calls.push_back(std::move(new_call));
        continue;
      }

      assert(b < static_cast<long>(buckets.size()));
      auto const & pileup = buckets[b].pileup;
      std::vector<long> new_phred(pl_len);

      if (pileup.size() == BUCKET_SIZE)
      {
        assert(p < static_cast<long>(pileup.size()));
        auto const & base_count = pileup[p];

        std::vector<long> seq_base2index;
        seq_base2index.reserve(cnum);

        for (long y{0}; y < cnum; ++y)
        {
          assert(y < static_cast<long>(variant.seqs.size()));
          assert(variant.seqs[y].size() == 1);

          seq_base2index.push_back(base2index(variant.seqs[y][0]));
        }

        // coverage
        for (long y{0}; y < 4; ++y)
        {
          auto find_it = std::find(seq_base2index.begin(), seq_base2index.end(), y);

          if (find_it == seq_base2index.end())
          {
            // not found, add ambigous coverage
            new_call.ambiguous_depth += base_count.acgt[y];
          }
          else
          {
            // found, add unique depth
            long const i = std::distance(seq_base2index.begin(), find_it);

            assert(i < static_cast<long>(new_call.coverage.size()));
            new_call.coverage[i] += base_count.acgt[y];
          }
        }

        long const total_qualsum = base_count.get_total_qualsum();

        // PHRED
        // auto max_it = std::max_element(base_count.acgt_qualsum.begin(), base_count.acgt_qualsum.end());
        // assert(max_it != base_count.acgt_qualsum.end());
        long i{0};
        // long constexpr ERROR_PHRED{25}; // Penalty of non-support

        // TODO use base_count.acgt_qualsum when calculating PL
        for (long y{0}; y < cnum; ++y)
        {
          for (long x{0}; x <= y; ++x, ++i)
          {
            assert(i < static_cast<long>(new_phred.size()));
            assert(i < static_cast<long>(new_call.phred.size()));
            assert(new_phred[i] == 0); // not previously set

            if (x == y)
            {
              long const x_idx = seq_base2index[y];
              new_phred[i] = total_qualsum - base_count.acgt_qualsum[x_idx]; // Total qualsum of all bases except x_idx
            }
            else
            {
              long const x_idx = seq_base2index[x];
              long const y_idx = seq_base2index[y];
              assert(x_idx != y_idx);

              new_phred[i] = total_qualsum - base_count.acgt_qualsum[x_idx] - base_count.acgt_qualsum[y_idx] +
                             3 * (base_count.acgt[x_idx] + base_count.acgt[y_idx]);
            }
          }
        }

        // Normalize phred score such that there is at least PHRED likelihood score equal to 0
        auto min_it = std::min_element(new_phred.begin(), new_phred.end());
        assert(min_it != new_phred.end());
        long const min_score = *min_it; // find the smallest score

        for (i = 0; i < pl_len; ++i)
        {
          long const pl = new_phred[i] - min_score; // substract the lowest PL to get a 0
          new_call.phred[i] = pl < 255 ? pl : 255;
        }
      }

      variant.calls.push_back(std::move(new_call));
    }

    variant.generate_infos();
    variant.infos.erase("MQ");

    // double QD, top avoid the vcf filter thingie
    /*
    if (variant.infos.count("QD") == 1)
    {
      double const old_qd = std::stod(variant.infos.at("QD"));
      variant.infos["QD"] = std::to_string(2.0 * old_qd);
    }*/

    vcf.variants.push_back(std::move(variant));
  }

  vcf.sample_names = std::move(new_sample_names);

#ifndef NDEBUG
  // print_info(__HERE__, " WIP done");
#endif // NDEBUG
}

} // namespace gyper
