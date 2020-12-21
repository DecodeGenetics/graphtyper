#include <cassert> // assert
#include <limits> // std::numeric_limits<T>
#include <memory> // std::unique_ptr
#include <sstream> // std::ostringstream
#include <string> // std::string
#include <unordered_map> // std::unordered_map
#include <vector> // std::vector

#include <paw/align.hpp>
#include <paw/station.hpp>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/hts_io.h>

#include <boost/log/trivial.hpp> // BOOST_LOG_TRIVIAL

#include <parallel_hashmap/phmap.h>

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
#include <graphtyper/utilities/hts_reader.hpp> // gyper::HtsReader
#include <graphtyper/utilities/io.hpp>
#include <graphtyper/utilities/options.hpp> // gyper::Options


namespace
{

#ifndef NDEBUG
std::string const debug_read_name = "HISEQ1:33:H9YY4ADXX:1:2110:2792:58362/2";
long debug_event_pos{699214};
char debug_event_type{'D'};
std::size_t debug_event_size{3};
#endif // NDEBUG


#ifndef NDEBUG
bool
check_haplotypes(std::map<gyper::Event, std::map<gyper::Event, int8_t> > const & haplotypes)
{
  bool is_ok{true};

  for (auto event_it = haplotypes.begin(); event_it != haplotypes.end(); ++event_it)
  {
    gyper::Event const & event = event_it->first;
    auto const & other_haps = event_it->second;
    long count{0};

    for (auto event_it2 = std::next(event_it);
         event_it2 != haplotypes.end() && event_it2->first.pos < (event.pos + 100l);
         ++event_it2, ++count)
    {
      gyper::Event const & other_event = event_it2->first;

      if (other_haps.count(other_event) == 0)
      {
        BOOST_LOG_TRIVIAL(warning) << __HERE__ << " Missing event link "
                                   << event.to_string() << " "
                                   << other_event.to_string();
        is_ok = false;
      }
    }

    if (count != static_cast<long>(other_haps.size()))
    {
      BOOST_LOG_TRIVIAL(warning) << __HERE__ << " count mismatch " << count << " != " << other_haps.size();
      BOOST_LOG_TRIVIAL(warning) << __HERE__ << " " << event.to_string();

      for (auto event_it2 = std::next(event_it);
           event_it2 != haplotypes.end() && event_it2->first.pos < (event.pos + 100l);
           ++event_it2)
      {
        gyper::Event const & other_event = event_it2->first;
        BOOST_LOG_TRIVIAL(info) << __HERE__ << " " << other_event.to_string();
      }

      for (auto const & o : other_haps)
      {
        BOOST_LOG_TRIVIAL(info) << __HERE__ << " " << o.first.to_string() << " "
                                << static_cast<long>(o.second);
      }

      is_ok = false;
    }
  }

  return is_ok;
}


#endif // NDEBUG


void
merge_haplotypes(std::map<gyper::Event, std::map<gyper::Event, int8_t> > & into,
                 std::map<gyper::Event, std::map<gyper::Event, int8_t> > & from)
{
  using namespace gyper;

  // go through into and find variants not in from, and add anti hap support accordingly
  for (auto it = into.begin(); it != into.end(); ++it)
  {
    Event const & event_into = it->first;
    std::map<gyper::Event, int8_t> & haplotype_into = it->second;

    auto find_it = from.find(event_into);

    if (find_it == from.end())
    {
      // BOOST_LOG_TRIVIAL(info) << __HERE__ << " Did not find " << event_into.to_string();

      for (auto & into2 : haplotype_into)
      {
        find_it = from.find(into2.first);

        if (find_it != from.end())
        {
          into2.second |= gyper::IS_ANY_ANTI_HAP_SUPPORT;
        } // else do nothing, the variant is phase completely irrelevant to "from"
      }
    }
    else
    {
      // found
      for (auto & into2 : haplotype_into)
      {
        auto find_it2 = haplotype_into.find(into2.first);

        if (find_it2 != haplotype_into.end())
        {
          // found
          into2.second |= find_it2->second;
        }
        else
        {
          into2.second |= gyper::IS_ANY_ANTI_HAP_SUPPORT;
        }
      }
    }
  }

  for (auto it = from.begin(); it != from.end(); ++it)
  {
    //Event const & first_from_event = it->first;
    auto insert_it = into.insert(*it);
    bool const is_inserted = insert_it.second;

    if (is_inserted)
    {
      //BOOST_LOG_TRIVIAL(info) << __HERE__ << " new event I hadnt seen before: " << it->first.to_string();
      // add anti-support for the event in previous events
      Event const & inserted_event = insert_it.first->first;
      std::map<Event, std::map<Event, int8_t> >::iterator prev_it(insert_it.first);

      // Go back to previous events and make sure they have this event
      while (prev_it != into.begin())
      {
        --prev_it;

        if ((prev_it->first.pos + 100l) <= inserted_event.pos)
          break;

        // insert it if it isn't already there
        if (prev_it->second.count(inserted_event) == 0)
          prev_it->second[inserted_event] = gyper::IS_ANY_ANTI_HAP_SUPPORT;
      }

      // also, go forward and add anti support to all (nearby) events that this one is missing
      for (auto forw_it = std::next(insert_it.first);
           forw_it != into.end() && forw_it->first.pos < (inserted_event.pos + 100l);
           ++forw_it)
      {
        Event const & next_event = forw_it->first;
        std::map<Event, int8_t> & inserted_hap = insert_it.first->second;

        // insert it if it isn't already there
        if (inserted_hap.count(next_event) == 0)
          inserted_hap[next_event] = gyper::IS_ANY_ANTI_HAP_SUPPORT;
      }
    }
    else
    {
      // We have seen this variant before
      std::map<Event, int8_t> & into_haplotypes = insert_it.first->second;
      std::map<Event, int8_t> const & from_haplotypes = it->second;

      // Add anti hap support to any variants not found in "from"
      for (auto into_hap_it = into_haplotypes.begin();
           into_hap_it != into_haplotypes.end();
           ++into_hap_it)
      {
        Event const & into_next_event = into_hap_it->first;
        int8_t & flags = into_hap_it->second;

        auto find_it = from_haplotypes.find(into_next_event);

        if (find_it == from_haplotypes.end())
          flags |= IS_ANY_ANTI_HAP_SUPPORT;
        else
          flags |= find_it->second;
      }

      for (auto from_hap_it = from_haplotypes.begin(); from_hap_it != from_haplotypes.end(); ++from_hap_it)
      {
        Event const & second_from_event = from_hap_it->first;
        auto find_snp_hap_it = into_haplotypes.find(second_from_event);

        if (find_snp_hap_it != into_haplotypes.end())
        {
          find_snp_hap_it->second |= from_hap_it->second; // Found
        }
        else
        {
          //// Not found
          // this means we saw the first event in a previous pool but not the second event,
          // therefore it is safe to say they are not on the same haplotype from that previous pool
          std::pair<Event, int8_t> new_value = *from_hap_it;
          new_value.second |= gyper::IS_ANY_ANTI_HAP_SUPPORT;
          into_haplotypes.insert(std::move(new_value));
        }
      }
    }
  }

  from.clear();
}


bool
is_clipped(bam1_t const & b)
{
  if (b.core.n_cigar == 0)
    return false;

  auto it = b.data + b.core.l_qname;

  // Check first
  uint32_t opAndCnt;
  memcpy(&opAndCnt, it, sizeof(uint32_t));

  if ((opAndCnt & 15) == 4)
  {
    //std::cerr << "MIDNSHP=X*******"[opAndCnt & 15] << " first ";
    return true;
  }

  // Check last
  memcpy(&opAndCnt, it + sizeof(uint32_t) * (b.core.n_cigar - 1), sizeof(uint32_t));

  if ((opAndCnt & 15) == 4)
  {
    //std::cerr << "MIDNSHP=X*******"[opAndCnt & 15] << " last ";
    return true;
  }

  return false;
}


void
_read_rg_and_samples(std::vector<std::string> & samples,
                     std::vector<std::unordered_map<std::string, int> > & vec_rg2sample_i,
                     std::vector<std::string> const & hts_paths)
{
  using namespace gyper;
  vec_rg2sample_i.resize(hts_paths.size());

  for (long i = 0; i < static_cast<long>(hts_paths.size()); ++i)
  {
    auto const & hts_path = hts_paths[i];
    auto & rg2sample_i = vec_rg2sample_i[i];
    std::size_t const old_size = samples.size();

    if (!Options::const_instance()->get_sample_names_from_filename)
    {
      gyper::get_sample_name_from_bam_header(hts_path, samples, rg2sample_i);
    }

    if (samples.size() == old_size)
    {
      std::string sample_name = hts_path.substr(hts_path.rfind('/') + 1, hts_path.rfind('.'));

      if (std::count(sample_name.begin(), sample_name.end(), '.') > 0)
        sample_name = sample_name.substr(0, sample_name.find('.'));

      samples.push_back(std::move(sample_name));
    }
  }
}


void
_determine_num_jobs_and_num_parts(long & jobs,
                                  long & num_parts,
                                  long const NUM_SAMPLES)
{
  using namespace gyper;

  gyper::Options const & copts = *(Options::const_instance());
  jobs = copts.threads;
  num_parts = jobs;
  long const MAX_FILES_OPEN = copts.max_files_open > jobs ? copts.max_files_open : jobs;

  if (jobs >= NUM_SAMPLES)
  {
    // Special case where there are more threads than samples OR more threads than max files open. POOL_SIZE is 1
    num_parts = std::min(NUM_SAMPLES, MAX_FILES_OPEN);
    jobs = num_parts;
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


} // anon namespace


namespace gyper
{

std::vector<std::string>
call(std::vector<std::string> const & hts_paths,
     std::vector<double> const & avg_cov_by_readlen,
     std::string const & graph_path,
     PHIndex const & ph_index,
     std::string const & output_dir,
     std::string const & reference_fn,
     std::string const & region,
     Primers const * primers,
     std::map<std::pair<uint16_t, uint16_t>, std::map<std::pair<uint16_t, uint16_t>, int8_t> > & ph,
     long const minimum_variant_support,
     double const minimum_variant_support_ratio,
     bool const is_writing_calls_vcf,
     bool const is_discovery,
     bool const is_writing_hap)
{
  std::vector<std::string> paths;

  if (hts_paths.size() == 0)
  {
    BOOST_LOG_TRIVIAL(error) << __HERE__ << " No input BAM files.";
    std::exit(1);
  }

  if (graph_path.size() > 0)
    load_graph(graph_path); // Loads the graph into the global variable 'graph'

  // Split hts_paths
  std::vector<std::unique_ptr<std::vector<std::string> > > spl_hts_paths;
  std::vector<std::unique_ptr<std::vector<double> > > spl_avg_cov;
  std::vector<std::unique_ptr<std::map<std::pair<uint16_t, uint16_t>,
                                       std::map<std::pair<uint16_t, uint16_t>, int8_t> > > > spl_ph;
  assert(Options::const_instance()->max_files_open > 0);

  long jobs{1}; // Running jobs
  long num_parts{1}; // Number of parts to split the work into
  long const NUM_SAMPLES{static_cast<long>(hts_paths.size())};

  _determine_num_jobs_and_num_parts(jobs, num_parts, NUM_SAMPLES);

  {
    auto it = hts_paths.begin();
    auto cov_it = avg_cov_by_readlen.begin();

    for (long i = 0; i < num_parts; ++i)
    {
      long const advance = NUM_SAMPLES / num_parts + (i < (NUM_SAMPLES % num_parts));

      auto end_it = it + advance;
      assert(std::distance(hts_paths.begin(), end_it) <= NUM_SAMPLES);
      spl_hts_paths.emplace_back(new std::vector<std::string>(it, end_it));
      it = end_it;

      auto cov_end_it = cov_it + advance;
      spl_avg_cov.emplace_back(new std::vector<double>(cov_it, cov_end_it));
      cov_it = cov_end_it;

      spl_ph.emplace_back(new std::map<std::pair<uint16_t, uint16_t>, std::map<std::pair<uint16_t, uint16_t>,
                                                                               int8_t> >());
    }
  }

  BOOST_LOG_TRIVIAL(debug) << __HERE__ << " Number of pools = " << spl_hts_paths.size();
  long const NUM_POOLS = spl_hts_paths.size();
  paths.resize(NUM_POOLS);

  // Run in parallel
  {
    paw::Station call_station(jobs); // last parameter is queue_size

    // Push all but the last pool to the thread pool
    if (!is_discovery)
    {
      for (long i{0}; i < (NUM_POOLS - 1l); ++i)
      {
        call_station.add_work(parallel_reader_genotype_only,
                              &paths[i],
                              spl_hts_paths[i].get(),
                              spl_avg_cov[i].get(),
                              &output_dir,
                              &reference_fn,
                              &region,
                              &ph_index,
                              primers,
                              spl_ph[i].get(),
                              is_writing_calls_vcf,
                              is_writing_hap);
      }

      // Do the last pool on the current thread
      call_station.add_to_thread(jobs - 1,
                                 parallel_reader_genotype_only,
                                 &paths[NUM_POOLS - 1],
                                 spl_hts_paths[NUM_POOLS - 1].get(),
                                 spl_avg_cov[NUM_POOLS - 1].get(),
                                 &output_dir,
                                 &reference_fn,
                                 &region,
                                 &ph_index,
                                 primers,
                                 spl_ph[NUM_POOLS - 1].get(),
                                 is_writing_calls_vcf,
                                 is_writing_hap);
    }
    else
    {
      for (long i{0}; i < (NUM_POOLS - 1l); ++i)
      {
        call_station.add_work(parallel_reader_with_discovery,
                              &paths[i],
                              spl_hts_paths[i].get(),
                              &output_dir,
                              &reference_fn,
                              &region,
                              &ph_index,
                              primers,
                              minimum_variant_support,
                              minimum_variant_support_ratio,
                              is_writing_calls_vcf,
                              is_writing_hap);
      }

      // Do the last pool on the current thread
      call_station.add_to_thread(jobs - 1,
                                 parallel_reader_with_discovery,
                                 &paths[NUM_POOLS - 1],
                                 spl_hts_paths[NUM_POOLS - 1].get(),
                                 &output_dir,
                                 &reference_fn,
                                 &region,
                                 &ph_index,
                                 primers,
                                 minimum_variant_support,
                                 minimum_variant_support_ratio,
                                 is_writing_calls_vcf,
                                 is_writing_hap);
    }

    std::string thread_info = call_station.join();
    BOOST_LOG_TRIVIAL(info) << "Finished calling. Thread work: " << thread_info;
  }

  if (!is_discovery && is_writing_hap)
  {
    assert(spl_ph.size() > 0);
    assert(spl_ph[0]);

    // gather ph
    ph = std::move(*(spl_ph[0].get()));

    for (long i{1}; i < static_cast<long>(spl_ph.size()); ++i)
    {
      // merge
      std::map<std::pair<uint16_t, uint16_t>, std::map<std::pair<uint16_t, uint16_t>, int8_t> > const & ph_parts
        = (*spl_ph[i].get());

      for (auto const & ph_part : ph_parts)
      {
        //std::pair<uint16_t, uint16_t> const & ph1 = ph_part.first;
        std::map<std::pair<uint16_t, uint16_t>, int8_t> const & ph2_map = ph_part.second;

        auto insert1_it = ph.insert(ph_part);

        if (insert1_it.second)
          continue; // not found

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
      }
    }
  }

  BOOST_LOG_TRIVIAL(debug) << __HERE__ << " Finished calling all samples.";
  return paths;
}


std::vector<VariantCandidate>
find_variants_in_cigar(seqan::BamAlignmentRecord const & record,
                       GenomicRegion const & region,
                       std::string const & ref)
{
  std::vector<VariantCandidate> new_var_candidates;
  assert(record.beginPos != -1);
  long ref_abs_pos = absolute_pos.get_absolute_position(region.chr, record.beginPos + 1);
  long const reference_offset = graph.ref_nodes.size() > 0 ?
                                graph.ref_nodes[0].get_label().order +
                                graph.genomic_region.get_absolute_position(1) -
                                1 :
                                0;

  long read_pos{0};
  std::string read(seqan::begin(record.seq), seqan::end(record.seq));
  long begin_clipping_size{0};
  bool is_clipped{false};
  std::string reference_seq;
  std::string read_seq;

  for (long c = 0; c < static_cast<long>(seqan::length(record.cigar)); ++c)
  {
    auto const & cigar = record.cigar[c];

    if (cigar.operation == 'M')
    {
      /// Add both read and reference
      // Get reference sequence overlapping this cigar
      reference_seq += ref.substr(ref_abs_pos - reference_offset, cigar.count);

      // Get read sequence overlapping this cigar
      read_seq += read.substr(read_pos, cigar.count);

      // Matches/mismatches advance both the read and the reference position
      read_pos    += cigar.count;
      ref_abs_pos += cigar.count;
    }
    else if (cigar.operation == 'I')
    {
      // Ignore if we haven't moved in the read
      if (read_pos > 0)
      {
        reference_seq.append(cigar.count, '-');
        read_seq += read.substr(read_pos, cigar.count);
      }

      // Insertions advance only the read position
      read_pos += cigar.count;
    }
    else if (cigar.operation == 'D')
    {
      // Ignore if we haven't moved in the read
      if (read_pos > 0)
      {
        reference_seq += ref.substr(ref_abs_pos - reference_offset, cigar.count);
        read_seq.append(cigar.count, '-');
      }

      // Deletions advance only the reference position
      ref_abs_pos += cigar.count;
    }
    else if (cigar.operation == 'S')
    {
      // soft-clips advance only the read position
      read_pos += cigar.count;
      is_clipped = true;

      if (c == 0)
        begin_clipping_size = cigar.count;
    }
    // else hard-clips ('H') advance neither the read nor the reference position
  }

  if (reference_seq == read_seq)
  {
    return new_var_candidates;
  }

  long ref_to_seq_offset = 0;

  // Reset to original ref_abs_pos
  ref_abs_pos = absolute_pos.get_absolute_position(region.chr, record.beginPos + 1);

  Variant new_var =
    make_variant_of_gapped_strings(reference_seq, read_seq, ref_abs_pos, ref_to_seq_offset);

  if (new_var.abs_pos == 0)
    return new_var_candidates;

  std::vector<Variant> new_vars =
    extract_sequences_from_aligned_variant(std::move(new_var), SPLIT_VAR_THRESHOLD - 1);

  new_var_candidates.resize(new_vars.size());

  for (unsigned i = 0; i < new_vars.size(); ++i)
  {
    assert(i < new_var_candidates.size());
    auto & new_var = new_vars[i];
    auto & new_var_candidate = new_var_candidates[i];
    assert(new_var.seqs.size() >= 2);
    new_var.normalize();

    //// Check if high or low quality
    long r = new_var.abs_pos - ref_to_seq_offset;

    // Fix r if the read was clipped
    if (length(record.cigar) > 0 && record.cigar[0].operation == 'S')
      r += record.cigar[0].count;

    long r_end = r + new_var.seqs[1].size();
    ref_to_seq_offset += static_cast<long>(new_var.seqs[0].size()) - static_cast<long>(new_var.seqs[1].size());

    new_var_candidate.abs_pos = new_var.abs_pos;
    new_var_candidate.original_pos = static_cast<uint32_t>(record.beginPos + 1 - begin_clipping_size);
    new_var_candidate.seqs = std::move(new_var.seqs);
    new_var_candidate.flags = record.flag;

    //// In case no qual is available
    if (static_cast<long>(seqan::length(record.qual)) > r)
    {
      int MAX_QUAL;

      if (r_end < static_cast<long>(seqan::length(record.qual)))
      {
        MAX_QUAL = *std::max_element(seqan::begin(record.qual) + r,
                                     seqan::begin(record.qual) + r_end) - 33;
      }
      else
      {
        MAX_QUAL = *std::max_element(seqan::begin(record.qual) + r,
                                     seqan::end(record.qual)) - 33;
      }

      new_var_candidate.flags |= static_cast<uint16_t>(static_cast<bool>(MAX_QUAL < 25)) << IS_LOW_BASE_QUAL_SHIFT;
    }

    //new_var_candidate.flags |= is_in_proper_pair = seqan::hasFlagAllProper(record);
    new_var_candidate.flags |= static_cast<uint16_t>(record.mapQ < 25) << IS_MAPQ_BAD_SHIFT;
    new_var_candidate.flags |= static_cast<uint16_t>(is_clipped) << IS_CLIPPED_SHIFT;
  }

  return new_var_candidates;
}


void
parallel_discover_from_cigar(std::string * output_ptr,
                             std::vector<std::string> const * hts_paths_ptr,
                             GenomicRegion const & region,
                             std::string const & output_dir,
                             std::string const & ref_str,
                             long minimum_variant_support,
                             double minimum_variant_support_ratio)
{
  assert(output_ptr);
  assert(hts_paths_ptr);
  auto const & hts_paths = *hts_paths_ptr;

  if (ref_str.size() == 0)
  {
    BOOST_LOG_TRIVIAL(error) << "Trying to discover variants with no reference string";
    std::exit(1);
  }

  // Determine the size of the region we are discovery variants on
  std::size_t const REGION_SIZE = region.end - region.begin;

  // Extract sample names from SAM
  std::vector<std::string> samples;
  std::vector<std::unordered_map<std::string, int> > vec_rg2sample_i; // Read group to sample index

  // Gather all the sample names
  _read_rg_and_samples(samples, vec_rg2sample_i, hts_paths);
  assert(samples.size() > 0);

  // Set up reference depth tracks
  ReferenceDepth reference_depth;
  reference_depth.set_depth_sizes(samples.size(), REGION_SIZE);

  // Set up variant map
  VariantMap varmap;
  varmap.set_samples(samples);
  varmap.minimum_variant_support = minimum_variant_support;
  varmap.minimum_variant_support_ratio = minimum_variant_support_ratio;

  for (long file_i = 0; file_i < static_cast<long>(hts_paths.size()); ++file_i)
  {
    assert(vec_rg2sample_i.size() == hts_paths.size());

    auto const & sam = hts_paths[file_i];
    auto & rg2sample_i = vec_rg2sample_i[file_i];

    seqan::HtsFile hts(sam.c_str(), "r");
    seqan::BamAlignmentRecord record;

    while (seqan::readRecord(record, hts))
    {
      // Skip supplementary, secondary and QC fail, duplicated and unmapped reads
      if (((record.flag & Options::const_instance()->sam_flag_filter) != 0) ||
          seqan::hasFlagUnmapped(record))
      {
        continue;
      }

      // Skip reads that are clipped on both ends. It is also fine to skip reads with no cigar since they cannot be useful here
      if (length(record.cigar) == 0 ||
          (record.cigar[0].operation == 'S' &&
           record.cigar[length(record.cigar) - 1].operation == 'S'))
      {
        continue;
      }

      if (Options::const_instance()->is_discovery_only_for_paired_reads &&
          (((record.flag & IS_PAIRED) == 0) ||
           ((record.flag & IS_PROPER_PAIR) == 0)))
      {
        continue;
      }

      assert(record.beginPos >= 0);

      // 0-based positions
      int64_t begin_pos = absolute_pos.get_absolute_position(region.chr, record.beginPos + 1);
      int64_t end_pos = begin_pos + seqan::getAlignmentLengthInRef(record);

      // Check if read is within region
      if (begin_pos < region.get_absolute_begin_position() || end_pos > region.get_absolute_end_position())
        continue;

      assert(begin_pos >= 0);
      assert(end_pos >= 0);

      // Determine the sample index
      assert(hts.hts_record);

      long sample_i;

      if (rg2sample_i.size() > 0)
      {
        uint8_t * rg_tag = bam_aux_get(hts.hts_record, "RG");

        if (!rg_tag)
        {
          BOOST_LOG_TRIVIAL(error) << __HERE__ << " Unable to find RG tag in read.";
          std::exit(1);
        }

        std::string read_group(reinterpret_cast<char *>(rg_tag + 1)); // Skip 'Z'

        auto find_rg_it = rg2sample_i.find(read_group);

        if (find_rg_it == rg2sample_i.end())
        {
          BOOST_LOG_TRIVIAL(error) << __HERE__ << " Unable to find read group. " << read_group;
          std::exit(1);
        }

        sample_i = find_rg_it->second;
      }
      else
      {
        sample_i = file_i;
      }

      // Update reference depth (very arbitrarily)
      if (end_pos - begin_pos < 50)
        reference_depth.add_depth(begin_pos, end_pos, sample_i);
      else
        reference_depth.add_depth(begin_pos + 4, end_pos - 4, sample_i);

      // Add variant candidates
      std::vector<VariantCandidate> var_candidates = find_variants_in_cigar(record, region, ref_str);

      if (var_candidates.size() > 0)
      {
        varmap.add_variants(std::move(var_candidates), sample_i);
      }
    }
  }

  // Write variant map
  std::ostringstream variant_map_path;
  variant_map_path << output_dir << "/" << samples[0] << "_variant_map";

  BOOST_LOG_TRIVIAL(debug) << __HERE__ << " Writing variant map to '" << variant_map_path.str();
  varmap.create_varmap_for_all(reference_depth);

#ifndef NDEBUG
  if (Options::const_instance()->stats.size() > 0)
    varmap.write_stats("1");
#endif // NDEBUG

  save_variant_map(variant_map_path.str(), varmap);
  *output_ptr = variant_map_path.str();
}


std::vector<std::string>
discover_directly_from_bam(std::string const & graph_path,
                           std::vector<std::string> const & hts_paths,
                           std::string const & region_str,
                           std::string const & output_dir,
                           long minimum_variant_support,
                           double minimum_variant_support_ratio)
{
  std::vector<std::string> output_file_paths;

  if (graph_path.size() > 0)
    load_graph(graph_path);

  //graph.generate_reference_genome();
  std::string ref_str(graph.reference.begin(), graph.reference.end());

  // parse region
  GenomicRegion const region(region_str);

  std::vector<std::unique_ptr<std::vector<std::string> > > spl_hts_paths;
  long jobs = 1;

  {
    long num_parts = 1;
    long const NUM_FILES = hts_paths.size();
    _determine_num_jobs_and_num_parts(jobs, num_parts, NUM_FILES);

    {
      auto it = hts_paths.cbegin();

      for (long i = 0; i < num_parts; ++i)
      {
        auto end_it = it + NUM_FILES / num_parts + (i < (NUM_FILES % num_parts));
        assert(std::distance(hts_paths.cbegin(), end_it) <= NUM_FILES);
        spl_hts_paths.emplace_back(new std::vector<std::string>(it, end_it));
        it = end_it;
      }
    }
  }

  long const NUM_POOLS = spl_hts_paths.size();
  BOOST_LOG_TRIVIAL(debug) << __HERE__ << " Number of pools = " << NUM_POOLS;
  output_file_paths.resize(NUM_POOLS);

  // Run in parallel
  {
    paw::Station call_station(jobs); // last parameter is queue_size
    //call_station.options.verbosity = 2; // Print messages

    for (long i = 0; i < NUM_POOLS - 1; ++i)
    {
      call_station.add_work(parallel_discover_from_cigar,
                            &(output_file_paths[i]),
                            spl_hts_paths[i].get(),
                            region,
                            output_dir,
                            ref_str,
                            minimum_variant_support,
                            minimum_variant_support_ratio);
    }

    // Do the last pool on the current thread
    call_station.add_to_thread(jobs - 1,
                               parallel_discover_from_cigar,
                               &(output_file_paths[NUM_POOLS - 1]),
                               spl_hts_paths[NUM_POOLS - 1].get(),
                               region,
                               output_dir,
                               ref_str,
                               minimum_variant_support,
                               minimum_variant_support_ratio);

    std::string thread_work_info = call_station.join();
    BOOST_LOG_TRIVIAL(info) << "Finished initial variant discovery step. Thread work info: " << thread_work_info;
  }

  return output_file_paths;
}


void
run_first_pass(bam1_t * hts_rec,
               HtsReader & hts_reader,
               long file_i,
               bool const is_first_in_pool,
               std::vector<BucketFirstPass> & buckets,
               std::map<Event, std::map<Event, int8_t> > & pool_haplotypes,
               long const BUCKET_SIZE,
               long const region_begin,
               std::vector<char> const & reference_sequence)
{
  long const REF_SIZE{static_cast<long>(reference_sequence.size())};
  int32_t global_max_pos_end{0};
  std::vector<uint32_t> cov_up(REF_SIZE);
  std::vector<uint32_t> cov_down(REF_SIZE);
  std::map<Event, std::map<Event, int8_t> > sample_haplotypes;

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

    std::array<char, 16> static constexpr CIGAR_MAP = {{
      'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', 'B', '*', '*', '*', '*', '*', '*'
    }};

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
      BOOST_LOG_TRIVIAL(error) << __HERE__ << " Unexpected ref_offset = " << ref_offset;
      std::exit(1);
    }

    Read read;
    read.name = reinterpret_cast<char *>(it);
    //read.mate_pos = static_cast<int32_t>(core.mpos);

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
      case '=':       // '=' and 'X' are typically not used by aligners (at least not by default),
      case 'X':       // but we keep it here just in case
      {
        //assert((ref_offset + cigar_count - 1l) < REF_SIZE);
        //auto ref_it = reference_sequence.begin() + ref_offset;

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

          if (read_base == ref ||
              (ref != 'A' && ref != 'C' && ref != 'G' && ref != 'T') ||
              (read_base != 'A' && read_base != 'C' && read_base != 'G' && read_base != 'T'))
          {
            continue;
          }

          Event new_snp_event(ref_pos + region_begin, 'X', {read_base});

#ifndef NDEBUG
          if (debug_event_type == 'X' && new_snp_event.pos == debug_event_pos)
          {
            BOOST_LOG_TRIVIAL(info) << __HERE__ << " " << new_snp_event.to_string()
                                    << " in file_i=" << file_i
                                    << " core.pos=" << core.pos
                                    << " read=" << read.name;
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
              //assert(core.pos > event_support.uniq_pos1);
              event_support.uniq_pos2 = core.pos;
            }
          }
          else if (event_support.uniq_pos3 == -1
                  && event_support.uniq_pos2 != core.pos /*&& event_support.uniq_pos1 != core.pos*/)
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

        auto begin_it = read_offset < static_cast<long>(read.sequence.size()) ?
                        read.sequence.begin() + read_offset :
                        read.sequence.end();

        auto end_it = (read_offset + cigar_count) < static_cast<long>(read.sequence.size()) ?
                      read.sequence.begin() + read_offset + cigar_count :
                      read.sequence.end();

        if (begin_it == end_it)
          break;

        // Make sure all bases are ACGT
        if (std::all_of(read.sequence.begin() + read_offset,
                        end_it,
                        [](char c){
              return c == 'A' || c == 'C' || c == 'G' || c == 'T';
            }))
        {
          std::vector<char> event_sequence(read.sequence.begin() + read_offset, end_it);
          Event new_event = make_insertion_event(region_begin + ref_offset,
                                                 std::move(event_sequence));

          // Add to bucket events
          auto indel_event_it = add_indel_event_to_bucket(buckets,
                                                          std::move(new_event),
                                                          region_begin,
                                                          BUCKET_SIZE,
                                                          reference_sequence,
                                                          ref_offset);

#ifndef NDEBUG
          if (read.name == debug_read_name)
            BOOST_LOG_TRIVIAL(info) << __HERE__ << " new ins=" << indel_event_it->first.to_string();
#endif // NDEBUG

          auto & event_support = indel_event_it->second;
          ++event_support.hq_count;

          if (core.qual != 255 && core.qual > event_support.max_mapq)
            event_support.max_mapq = core.qual;

          event_support.proper_pairs += ((read.flags & IS_PROPER_PAIR) != 0);
          //event_support.first_in_pairs += ((read.flags & IS_FIRST_IN_PAIR) != 0);
          event_support.sequence_reversed += ((read.flags & IS_SEQ_REVERSED) != 0);
          event_support.clipped += is_read_clipped;

          read_offset += cigar_count;
          cigar_events.push_back(indel_event_it);
        }

        break;
      }

      case 'D':
      {
        assert(cigar_count > 0);

        if (ref_offset + cigar_count >= REF_SIZE)
          break;

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
          BOOST_LOG_TRIVIAL(info) << __HERE__ << " new del=" << indel_event_it->first.to_string();
#endif // NDEBUG

        auto & event_support = indel_event_it->second;
        ++event_support.hq_count;

        if (core.qual != 255 && core.qual > event_support.max_mapq)
          event_support.max_mapq = core.qual;

        event_support.proper_pairs += ((read.flags & IS_PROPER_PAIR) != 0);
        //event_support.first_in_pairs += ((read.flags & IS_FIRST_IN_PAIR) != 0);
        event_support.sequence_reversed += ((read.flags & IS_SEQ_REVERSED) != 0);
        event_support.clipped += is_read_clipped;
        ref_offset += cigar_count;
        cigar_events.push_back(indel_event_it);
        break;
      }

      case 'S':
      {
        read_offset += cigar_count;
        break;
      } // else 'H' or 'P' which move neither reference nor read
      }
    }

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

    read.alignment.pos = static_cast<int32_t>(core.pos);
    read.alignment.pos_end = region_begin + std::min(ref_offset, REF_SIZE - 1);

    assert((read.alignment.pos - region_begin) < static_cast<long>(cov_up.size()));
    ++cov_up[read.alignment.pos - region_begin];

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
      BOOST_LOG_TRIVIAL(warning) << __HERE__ << " file_i=" << file_i << " is not sorted. "
                                 << hts_rec->core.pos << " < " << prev_pos;
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

  auto update_coverage =
    [&cov_up, &cov_down](long & cov, long const pos, long const b, long const BUCKET_SIZE) -> void
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

      for (auto snp_it = bucket.events.begin(); snp_it != bucket.events.end();)  // no increment
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
            BOOST_LOG_TRIVIAL(info) << __HERE__ << " debug snp has good support " << snp.to_string() << " "
                                    << info.to_string() << " file_i=" << file_i;
          }
#endif // DEBUG
          ++snp_it;
        }
        else
        {
//#ifndef NDEBUG
//          if (debug_event_type == 'X' && snp.pos == debug_event_pos)
//          {
//            BOOST_LOG_TRIVIAL(info) << __HERE__ << " debug snp has bad support " << snp.to_string() << " "
//                                    << info.to_string() << " file_i=" << file_i;
//          }
//#endif // DEBUG

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
                                static_cast<double>(indel.sequence.size() / 2.0 + 8.0) / 8.0 :
                                static_cast<double>(indel.sequence.size() / 3.0 + 10.0) / 10.0;

      double const count = correction * (indel_info.hq_count + indel_info.lq_count);

#ifndef NDEBUG
      if (debug_event_type == indel.type && debug_event_pos == indel.pos && debug_event_size == indel.sequence.size())
      {
        BOOST_LOG_TRIVIAL(info) << __HERE__ << "Indel [pos,pos+span], size [begin,end]: "
                                << indel.to_string() << " ["
                                << region_begin << " "
                                <<  indel.pos << "," << (indel.pos + indel_info.span)
                                << "] " << indel.sequence.size() << " " << correction
                                << " [" << naive_begin << "," << naive_end << "]"
                                << " count=" << count << " depth=" << depth;
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
          BOOST_LOG_TRIVIAL(info) << __HERE__ << " Indel has good support in file_i,log_qual="
                                  << file_i << "," << log_qual << "," << count << "," << anti_count_d
                                  << " cov=" << cov;
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
          BOOST_LOG_TRIVIAL(info) << __HERE__ << " Indel has realignment support in file_i,log_qual,count,acount="
                                  << file_i << "," << log_qual << "," << count << "," << anti_count_d
                                  << " cov=" << cov;
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
          BOOST_LOG_TRIVIAL(info) << __HERE__ << " erasing indel with bad support in file_i,log_qual,count,acount="
                                  << file_i << "," << log_qual << "," << count << "," << anti_count_d
                                  << " cov=" << cov;
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
    auto const & bucket = buckets[b];

    // Check event phase
    for (auto event_it = bucket.events.begin(); event_it != bucket.events.end(); ++event_it)
    {
      Event const & event = event_it->first;
      EventSupport const & info = event_it->second;

      long const begin = std::max(0l, event.pos - region_begin);
      long cov{depth}; // starting coverage depth

      auto it_bool = sample_haplotypes.insert({event, std::map<Event, int8_t>()});
      std::map<Event, int8_t> & haplotypes = it_bool.first->second;

      assert(begin >= b * BUCKET_SIZE);
      update_coverage(cov, begin, b, BUCKET_SIZE);
      double support_ratio = static_cast<double>(info.get_raw_support()) / static_cast<double>(cov);

      if (support_ratio < 0.3)
        support_ratio = 0.3;

      // 0 low coverage or ambigous
      // 1 support
      // 2 anti support
      auto is_good_support =
        [&cov_down, &region_begin, &support_ratio]
          (long local_cov,
          long local_offset,
          std::map<Event, EventSupport>::const_iterator event_it,
          std::map<Event, EventSupport>::const_iterator event_it2,
          std::map<Event, uint16_t>::const_iterator find_it,
          std::map<Event, uint16_t> const & map) -> uint16_t
        {
          bool is_indel = event_it->first.type != 'X' || event_it2->first.type != 'X';

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

          return 0;
        };

      // check this bucket
      for (auto event_it2 = std::next(event_it); event_it2 != bucket.events.end(); ++event_it2)
      {
        Event const & other_event = event_it2->first;

        if (other_event.pos == event.pos && other_event.type == event.type)
        {
          // If they share a position and type they trivially cant share haplotype
          haplotypes[other_event] |= IS_ANY_ANTI_HAP_SUPPORT;
          continue;
        }

        auto find_it = info.phase.find(other_event);
        //assert(other_event.pos > event.pos);
        uint16_t const flags = is_good_support(cov, begin + 1, event_it, event_it2, find_it, info.phase);
        haplotypes[other_event] |= flags;
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
          //if (other_event.pos >= (event.pos + 2 * BUCKET_SIZE))
          //  break;

          auto find_it = info.phase.find(other_event);
          uint16_t const flags = is_good_support(cov, begin + 1, event_it, event_it2, find_it, info.phase);
          haplotypes[other_event] |= flags;
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
          haplotypes[other_event] |= flags;
        }
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
  } // for (long b{0}; b < static_cast<long>(buckets.size()); ++b)

  // Add sample haplotypes to pool haplotypes
#ifndef NDEBUG
  // check for errors in sample_haplotypes
  assert(check_haplotypes(sample_haplotypes));
#endif // NDEBUG

  if (is_first_in_pool)
  {
    pool_haplotypes = std::move(sample_haplotypes);
  }
  else
  {
    merge_haplotypes(pool_haplotypes, sample_haplotypes);
    assert(check_haplotypes(pool_haplotypes));
  }
}


void
realign_to_indels(std::vector<Tindel_events::iterator> const & realignment_indels,
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
    if (debug_event_type == indel.type &&
        debug_event_pos == indel.pos &&
        debug_event_size == indel.sequence.size())
    {
      BOOST_LOG_TRIVIAL(info) << __HERE__ << " Realignment to indel=" << indel.to_string() << " span="
                              << indel_info.span;
    }
#endif // NDEBUG

    //create reference sequence with the indel
    long const begin_padded = std::max(0l, indel.pos - max_read_size - 2 * PAD - region_begin);
    assert(begin_padded < REF_SIZE);
    long const end_padded = indel.pos + max_read_size + 2 * PAD - region_begin;

    auto end_it = end_padded >= REF_SIZE ?
                  reference_sequence.end() :
                  reference_sequence.begin() + end_padded;

    std::vector<char> new_ref(reference_sequence.begin() + begin_padded, end_it);
    std::vector<int32_t> ref_pos(new_ref.size());
    std::iota(ref_pos.begin(), ref_pos.end(), 0);
    //BOOST_LOG_TRIVIAL(info) << "New ref before:\n" << std::string(new_ref.begin(), new_ref.end());

    {
      bool const is_applied = apply_indel_event(new_ref, ref_pos, indel, begin_padded + region_begin);
      assert(is_applied);

      if (!is_applied)
        continue;
    }

    //BOOST_LOG_TRIVIAL(info) << "New ref after:\n" << std::string(new_ref.begin(), new_ref.end());
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
    for ( ; b <= b_end; ++b)
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
          BOOST_LOG_TRIVIAL(info) << __HERE__ << " Found read=" << debug_read_name;
        }
#endif // NDEBUG

        if (read.alignment.pos < 0 || read.sequence.size() == 0)
          continue;

        // If the variant already has the indel then we can just update the score and skip the realignment
        if (read.alignment.has_indel_event(indel_it))
          continue;

        if ((read.alignment.num_clipped_end == 0 && read.alignment.pos_end < static_cast<long>(indel.pos)) ||
            (read.alignment.pos_end + read.alignment.num_clipped_end +
             std::min(static_cast<long>(read.alignment.num_clipped_end), PAD) < indel.pos) ||
            (read.alignment.num_clipped_begin == 0 && read.alignment.pos > indel_span) ||
            (read.alignment.pos - read.alignment.num_clipped_begin -
             std::min(static_cast<long>(read.alignment.num_clipped_begin), PAD) > indel_span))
        {
          continue;
        }

        //if (read.alignment.score == static_cast<long>(read.sequence.size())) // Cannot improve score
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
                BOOST_LOG_TRIVIAL(info) << __HERE__ << " Applied event " << it->event_it->first.to_string();
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
                BOOST_LOG_TRIVIAL(info) << __HERE__ << " Could not apply event " << it->event_it->first.to_string();
              }
#endif // NDEBUG

              applied_events.push_back({READ_ANTI_SUPPORT, it->event_it});
            }
          }
        }

        long const old_score = read.alignment.score;

        // TODO write greedy alignment
        paw::global_alignment(read.sequence, new_ref, opts);
#ifndef NDEBUG
        ++_alignment_counter;
#endif // NDEBUG
        assert(opts.get_alignment_results());
        auto const clip = opts.get_alignment_results()->apply_clipping(read.sequence,
                                                                       new_ref,
                                                                       opts.get_match(),
                                                                       opts.get_mismatch(),
                                                                       opts.get_gap_open(),
                                                                       opts.get_gap_extend(),
                                                                       opts.get_clip());

        paw::AlignmentResults<Tuint> const & ar = *opts.get_alignment_results();
        auto p = ar.get_database_begin_end(read.sequence, new_ref);

        if (p.first == 0 || p.second == ar.database_end)
        {
          //auto aligned_strings = ar.get_aligned_strings(read.sequence, new_ref);
          //
          //BOOST_LOG_TRIVIAL(warning) << __HERE__ << " Insufficient padding in read alignment, skipping read. "
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
          //if (clip.first > 0 || clip.second < static_cast<long>(read.sequence.size()))
          {
            BOOST_LOG_TRIVIAL(debug) << __HERE__ << " name=" << read.name << " pos=" << read.alignment.pos
                                     << " begin,end_clip="  << clip.first << "," << clip.second
                                     << " old_score=" << old_score << " new_score=" << ar.score;
          }
        }
#endif // NDEBUG

        if (ar.score <= old_score)
        {
          // If the score is worse and the sequence is not clipped it is anti-support
          if (ar.score < old_score)
          {
#ifndef NDEBUG
            if (read.name == debug_read_name)
            {
              auto aligned_strings = ar.get_aligned_strings(read.sequence, new_ref);

              BOOST_LOG_TRIVIAL(debug) << __HERE__ << " ANTI old_score, new_score, read = " << old_score << ", "
                                       << ar.score << ", "
                                       << read.name << '\n' << aligned_strings.first << '\n'
                                       << aligned_strings.second;
            }
#endif // NDEBUG

            read.alignment.add_indel_event(READ_ANTI_SUPPORT, read.flags, read.mapq, indel_it);
          }
          else if (indel_it->first.pos >= (ref_pos[p.first] + begin_padded) &&
                   indel_it->first.pos <= (ref_pos[p.second] + begin_padded))
          {
            assert(ar.score == old_score);
#ifndef NDEBUG
            if (read.name == debug_read_name)
            {
              BOOST_LOG_TRIVIAL(debug) << __HERE__ << " SAME SCORE " << old_score
                                       << " read=" << read.to_string();
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
          long const _next_alignment_pos = ref_pos[p.first] + region_begin + begin_padded;

          if (read.alignment.pos > (_next_alignment_pos + 300) ||
              _next_alignment_pos > (read.alignment.pos + 300))
          {
            BOOST_LOG_TRIVIAL(warning) << __HERE__ << " change pos from="
                                       << read.alignment.pos << " to="
                                       << _next_alignment_pos
                                       << " ("
                                       << ref_pos[p.first] << " + " << region_begin << " + " << begin_padded
                                       << ")";
          }
        }
#endif // NDEBUG

        read.alignment.pos = ref_pos[p.first] + region_begin + begin_padded;
        read.alignment.pos_end = ref_pos[p.second] + region_begin + begin_padded;
        read.alignment.score = ar.score;
        read.alignment.num_clipped_begin = clip.first;
        read.alignment.num_clipped_end = static_cast<long>(read.sequence.size()) - clip.second;

        // count number of inserted bases at the beginning of the read
        {
          long num_ins{0};

          while (p.first + num_ins + 1 < static_cast<long>(ref_pos.size()) &&
                 ref_pos[p.first + num_ins] == ref_pos[p.first + num_ins + 1])
          {
            ++num_ins;
          }

          read.alignment.num_ins_begin = num_ins;
        }

#ifndef NDEBUG
        if (read.name == debug_read_name)
        {
          auto aligned_strings = ar.get_aligned_strings(read.sequence, new_ref);
          BOOST_LOG_TRIVIAL(info) << __HERE__ << " score=" << ar.score << "\n"
                                  << aligned_strings.first << "\n" << aligned_strings.second;
        }
#endif // NDEBUG

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
  BOOST_LOG_TRIVIAL(info) << __HERE__ << " alignment counter = " << _alignment_counter;
#endif // NDEBUG
  */
  // Print counts etc
  for (Tindel_events::iterator indel_it : realignment_indels)
  {
    Event const & indel = indel_it->first;
    EventSupport & indel_info = indel_it->second;

    if (indel_info.has_indel_good_support)
      continue;

    double const correction = indel.type == 'I' ?
                              static_cast<double>(indel.sequence.size() / 2.0 + 8.0) / 8.0 :
                              static_cast<double>(indel.sequence.size() / 3.0 + 10.0) / 10.0;

    double const count = correction * (indel_info.hq_count + indel_info.lq_count);

    bool const is_good_count = (indel_info.hq_count >= 5 && count >= 8.0) ||
                               (indel_info.span >= 9 && indel_info.hq_count >= 4 && count >= 6.0) ||
                               (indel_info.span >= 18 && indel_info.hq_count >= 3 && count >= 4.0);

#ifndef NDEBUG
    if (debug_event_type == indel.type &&
        debug_event_pos == indel.pos &&
        debug_event_size == indel.sequence.size())
    {
      BOOST_LOG_TRIVIAL(info) << __HERE__ << " Realignment results for indel=" << indel.to_string() << " "
                              << indel_info.hq_count << "," << indel_info.anti_count
                              << " log_qual=" << indel_info.log_qual(10)
                              << " info=" << indel_info.to_string()
                              << " is_good_count=" << is_good_count
                              << " is_good_info=" << indel_info.is_good_indel();
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


void
read_hts_and_return_realignment_indels(bam1_t * hts_rec,
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

    std::array<char, 16> static constexpr CIGAR_MAP = {{
      'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', 'B', '*', '*', '*', '*', '*', '*'
    }};

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
      BOOST_LOG_TRIVIAL(error) << __HERE__ << " Unexpected ref_offset = " << ref_offset;
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
      BOOST_LOG_TRIVIAL(info) << __HERE__ << " Found debug read=" << debug_read_name;
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
        BOOST_LOG_TRIVIAL(warning) << __HERE__ << " While processing read="
                                   << reinterpret_cast<char *>(hts_rec->data)
                                   << " went beyond the region if interest";
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
      case '=':       // '=' and 'X' are typically not used by aligners (at least not by default),
      case 'X':       // but we keep it here just in case
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

        auto begin_it = read_offset < static_cast<long>(read.sequence.size()) ?
                        read.sequence.begin() + read_offset :
                        read.sequence.end();

        auto end_it = (read_offset + count) < static_cast<long>(read.sequence.size()) ?
                      read.sequence.begin() + read_offset + count :
                      read.sequence.end();

        if (begin_it == end_it)
        {
          break;
        }

        std::vector<char> event_sequence(begin_it, end_it);

        Event new_event = make_insertion_event(region_begin + ref_offset,
                                               std::move(event_sequence));

        auto indel_event_it =
          add_indel_event_to_bucket(buckets,
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
          BOOST_LOG_TRIVIAL(warning) << __HERE__ << " While processing read="
                                     << reinterpret_cast<char *>(hts_rec->data)
                                     << " went beyond the region if interest";
#endif // NDEBUG
          break;
        }

        Event new_event =
          make_deletion_event(reference_sequence, ref_offset, region_begin + ref_offset, count);

        auto indel_event_it =
          add_indel_event_to_bucket(buckets,
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
          //if ((i + 1 != N_CIGAR) && (i + 2 != N_CIGAR))
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


void
parallel_first_pass(std::vector<std::string> * hts_paths_ptr,
                    std::map<Event, std::map<Event, int8_t> > * pool_haplotypes_ptr,
                    std::vector<std::vector<BucketFirstPass> > * file_buckets_first_pass_ptr,
                    std::vector<std::string> * sample_names_ptr,
                    std::vector<char> * reference_sequence_ptr,
                    long const BUCKET_SIZE,
                    long const region_begin,
                    long const lowest_file_i)
{
  assert(hts_paths_ptr);
  assert(pool_haplotypes_ptr);

  std::vector<std::string> const & hts_paths = *hts_paths_ptr;
  std::map<Event, std::map<Event, int8_t> > & pool_haplotypes = *pool_haplotypes_ptr;
  std::vector<char> const & reference_sequence = *reference_sequence_ptr;

  for (long i{0}; i < static_cast<long>(hts_paths.size()); ++i)
  {
    // BAM/CRAM record storage
    HtsStore store;

    // Open BAM/CRAM file
    HtsReader hts_reader(store);
    hts_reader.open(hts_paths[i], std::string("."), "");

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
                   i == 0, // is_first_in_pool
                   (*file_buckets_first_pass_ptr)[lowest_file_i + i],
                   pool_haplotypes,
                   BUCKET_SIZE,
                   region_begin,
                   reference_sequence);

    // Add sample
    if (hts_reader.samples.size() > 1)
    {
      BOOST_LOG_TRIVIAL(error) << __HERE__ << " We found file with multiple samples, sorry, "
                               << "this is currently not supported.";
      std::exit(1);
    }

    assert(hts_reader.samples.size() == 1);
    (*sample_names_ptr)[lowest_file_i + i] = hts_reader.samples[0];

    // Close BAM/CRAM file
    hts_reader.close();
  } // for (long file_i{0}; file_i < static_cast<long>(hts_paths.size()); ++file_i)
}


void
parallel_second_pass(std::string const * hts_path_ptr,
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

  std::sort(nearby_good_events.begin(), nearby_good_events.end(),
            [](Tindel_events::iterator const & a,
               Tindel_events::iterator const & b) -> bool
    {
      return a->first.pos < b->first.pos;
    });

  auto end_unique_it = std::unique(nearby_good_events.begin(), nearby_good_events.end());
  std::move(nearby_good_events.begin(), end_unique_it, std::back_inserter(indel_to_realign));

  std::sort(indel_to_realign.begin(), indel_to_realign.end(),
            [](Tindel_events::iterator const & a,
               Tindel_events::iterator const & b) -> bool
    {
      return a->second.has_indel_good_support > b->second.has_indel_good_support ||
      (a->second.has_indel_good_support == b->second.has_indel_good_support && a->first.pos < b->first.pos);
    });

#ifndef NDEBUG
  BOOST_LOG_TRIVIAL(debug) << __HERE__ << " Nearby realignment indels="
                           << nearby_good_events.size();
#endif // NDEBUG

  realign_to_indels(indel_to_realign, //realignment_indels,
                    buckets,
                    max_read_size,
                    BUCKET_SIZE,
                    region_begin,
                    reference_sequence);

#ifndef NDEBUG
  BOOST_LOG_TRIVIAL(debug) << __HERE__ << " Realignment DONE.";
#endif // NDEBUG
}


void
streamlined_discovery(std::vector<std::string> const & hts_paths,
                      std::string const & ref_path,
                      std::string const & region_str,
                      gyper::Vcf & vcf)
{
  BOOST_LOG_TRIVIAL(info) << "Start WIP on region " << region_str;

  long const NUM_FILES = hts_paths.size();
  GenomicRegion genomic_region(region_str);
  uint64_t const chromosome_offset = genomic_region.get_absolute_position(1);
  long const region_begin{genomic_region.begin};
  std::vector<char> reference_sequence;
  open_and_read_reference_genome(reference_sequence, ref_path, genomic_region);
  long constexpr BUCKET_SIZE{50}; // TODO make an option
  std::vector<std::string> sample_names;
  sample_names.resize(NUM_FILES);

  std::vector<std::unique_ptr<std::vector<std::string> > > spl_hts_paths{};
  std::vector<std::unique_ptr<std::map<Event, std::map<Event, int8_t> > > > spl_snp_hq_haplotypes{};
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
        spl_snp_hq_haplotypes.emplace_back(new std::map<Event, std::map<Event, int8_t> >());
        it = end_it;
      }
    }
  }

#ifndef NDEBUG
  BOOST_LOG_TRIVIAL(info) << __HERE__ << " First pass starting.";
#endif // NDEBUG

  long const NUM_POOLS{static_cast<long>(spl_hts_paths.size())};

#ifndef NDEBUG
  BOOST_LOG_TRIVIAL(debug) << __HERE__ << " Number of pools = " << NUM_POOLS;
#endif // NDEBUG

  Tindel_events indel_events;
  long NUM_BUCKETS{0};
  std::map<Event, std::map<Event, int8_t> > haplotypes;

  // FIRST PASS
  {
    std::vector<std::vector<BucketFirstPass> > file_buckets_first_pass(NUM_FILES);

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
      std::map<Event, std::map<Event, int8_t> > & pool_haplotypes = *spl_snp_hq_haplotypes[i];
      assert(check_haplotypes(pool_haplotypes));
      merge_haplotypes(haplotypes, pool_haplotypes);
      assert(check_haplotypes(haplotypes));
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
          if (indel_event.first.type == 'X')
            continue;

#ifndef NDEBUG
          if (debug_event_type == indel_event.first.type &&
              debug_event_pos == indel_event.first.pos &&
              debug_event_size == indel_event.first.sequence.size())
          {
            BOOST_LOG_TRIVIAL(info) << __HERE__ << " Indel=" << indel_event.first.to_string()
                                    << " good_support=" << indel_event.second.has_indel_good_support;
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
  BOOST_LOG_TRIVIAL(info) << __HERE__ << " First pass DONE. Starting second pass.";
#endif // NDEBUG

  // SECOND PASS BEGINS
  {
    std::vector<std::vector<Tindel_events::iterator> > indel_to_realign;
    indel_to_realign.resize(NUM_FILES);

    // Add indels to events
    for (auto it = indel_events.begin(); it != indel_events.end(); ++it)
    {
      it->second.clear();

      //it->second.hq_count = 0;
      //it->second.lq_count = 0;
      //it->second.sequence_reversed = 0;
      //it->second.proper_pairs = 0;
      //it->second.max_mapq = 0;

      assert(it->second.anti_count == 0);
      assert(it->second.multi_count == 0);
      assert(it->second.has_realignment_support);

      if (!it->second.has_indel_good_support)
      {
#ifndef NDEBUG
        if (debug_event_type == it->first.type &&
            debug_event_pos == it->first.pos &&
            debug_event_size == it->first.sequence.size())
        {
          BOOST_LOG_TRIVIAL(info) << __HERE__ << " realignment indel=" << it->first.to_string()
                                  << " qual,file_i="
                                  << it->second.max_log_qual << " " << it->second.max_log_qual_file_i;
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
      //if (event_it->first.type != 'X')
      //  continue; // only SNPs here

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

      //BOOST_LOG_TRIVIAL(info) << __HERE__ << " SNP=" << snp_it->first.to_string();
      //for (auto ph_it = snp_it->second.begin(); ph_it != snp_it->second.end(); ++ph_it)
      //{
      //  BOOST_LOG_TRIVIAL(info) << __HERE__ << " phased with=" << ph_it->first.to_string()
      //                          << " any_support=" << ((ph_it->second & IS_ANY_HAP_SUPPORT) != 0)
      //                          << " any_anti_support=" << ((ph_it->second & IS_ANY_ANTI_HAP_SUPPORT) != 0);
      //}

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

          auto find_event_it = event_it->second.find(next_event);

          if (find_event_it != event_it->second.end())
          {
            // found
            if (find_event_it->second == IS_ANY_HAP_SUPPORT)
            {
              // ONLY haplotype support
              if (is_hap_stream_empty)
                is_hap_stream_empty = false;
              else
                ss_hap << ',';

              ss_hap << next_event_index;
            }
            else if (find_event_it->second == IS_ANY_ANTI_HAP_SUPPORT)
            {
              // ONLY anti haplotype support
              if (is_anti_stream_empty)
                is_anti_stream_empty = false;
              else
                ss_anti << ',';

              ss_anti << next_event_index;
            }
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

#ifndef NDEBUG
  BOOST_LOG_TRIVIAL(info) << __HERE__ << " WIP done";
#endif // NDEBUG
}


} // namespace gyper
