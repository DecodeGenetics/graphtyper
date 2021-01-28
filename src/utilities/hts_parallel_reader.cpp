#include <string> // std::string
#include <unordered_map> // std::unordered_map
#include <vector> // std::vector

#include <boost/log/trivial.hpp>

#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <htslib/sam.h>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/hts_io.h>

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/haplotype_calls.hpp>
#include <graphtyper/graph/reference_depth.hpp>
#include <graphtyper/typer/alignment.hpp>
#include <graphtyper/typer/genotype_paths.hpp>
#include <graphtyper/typer/primers.hpp>
#include <graphtyper/typer/variant_map.hpp>
#include <graphtyper/typer/vcf.hpp>
#include <graphtyper/typer/vcf_writer.hpp>
#include <graphtyper/utilities/hash_seqan.hpp>
#include <graphtyper/utilities/hts_parallel_reader.hpp>
#include <graphtyper/utilities/hts_store.hpp>
#include <graphtyper/utilities/hts_writer.hpp>
#include <graphtyper/utilities/io.hpp>
#include <graphtyper/utilities/options.hpp>

namespace
{

using TMapGPaths = std::unordered_map<std::string, std::pair<gyper::GenotypePaths, gyper::GenotypePaths> >;

#ifndef NDEBUG
void
check_if_maps_are_empty(std::vector<TMapGPaths> const & maps)
{
  for (long i = 0; i < static_cast<long>(maps.size()); ++i)
  {
    auto const & map = maps[i];

    if (map.size() > 0)
    {
      BOOST_LOG_TRIVIAL(info)
        << __HERE__ << " Found a non-empty read map for with size "
        << map.size()
        << "! This likely means these reads have the BAM_FPAIRED flag set but have no mate read.";

      //for (auto it = map.begin(); it != map.end(); ++it)
      //  BOOST_LOG_TRIVIAL(debug) << "Leftover read name: " << it->first;
    }
  }
}


#endif // NDEBUG

} // anon namespace


namespace gyper
{

HtsParallelReader::~HtsParallelReader()
{
  close();
}


void
HtsParallelReader::open(std::vector<std::string> const & hts_file_paths,
                        std::string const & reference,
                        std::string const & region)
{
  for (auto const & bam : hts_file_paths)
  {
    HtsReader f(store);
    f.open(bam, region, reference);
    f.set_sample_index_offset(samples.size());
    std::copy(f.samples.begin(), f.samples.end(), std::back_inserter(samples));
    f.set_rg_index_offset(num_rg);
    num_rg += f.get_num_rg();
    hts_files.push_back(std::move(f));
  }

  // Output
  long const n = hts_files.size();

  // Read the first record of each hts file
  for (long i{0}; i < n; ++i)
  {
    bam1_t * hts_rec = hts_files[i].get_next_read_in_order();

    if (hts_rec)
    {
      heap.push_back(HtsRecord(hts_rec, i));
      std::push_heap(heap.begin(), heap.end(), Cmp_gt_pair_bam1_t_fun);
    }
  }
}


void
HtsParallelReader::close()
{
  for (auto & hts_f : hts_files)
    hts_f.close();
}


bool
HtsParallelReader::read_record(HtsRecord & hts_record)
{
  if (heap.empty())
    return false;

  auto const i = heap[0].file_index;
  bam1_t * old_record = hts_record.record;
  hts_record.record = heap[0].record;
  hts_record.file_index = i;
  bam1_t * hts_rec;

  // reuse an old bam1_t record if it is available
  if (old_record)
    hts_rec = hts_files[i].get_next_read_in_order(old_record);
  else
    hts_rec = hts_files[i].get_next_read_in_order();

  std::pop_heap(heap.begin(), heap.end(), Cmp_gt_pair_bam1_t_fun);
  assert(hts_record.record == heap.back().record);

  if (hts_rec)
  {
    (heap.end() - 1)->record = hts_rec;
    assert((heap.end() - 1)->file_index == i);
    std::push_heap(heap.begin(), heap.end(), Cmp_gt_pair_bam1_t_fun);
  }
  else
  {
    (heap.end() - 1)->record = nullptr;
    heap.pop_back();
  }

  return true;
}


void
HtsParallelReader::move_record(HtsRecord & to, HtsRecord & from)
{
  if (to.record)
    store.push(to.record);

  to.record = from.record;
  to.file_index = from.file_index;
  from.record = nullptr;
}


std::vector<std::string> const &
HtsParallelReader::get_samples() const
{
  return samples;
}


void
HtsParallelReader::get_sample_and_rg_index(long & sample_i, long & rg_i, HtsRecord const & hts_rec) const
{
  assert(hts_rec.file_index < static_cast<long>(hts_files.size()));
  auto const & hts_f = hts_files[hts_rec.file_index];
  hts_f.get_sample_and_rg_index(sample_i, rg_i, hts_rec.record);
}


long
HtsParallelReader::get_num_rg() const
{
  return num_rg;
}


long
HtsParallelReader::get_num_samples() const
{
  return samples.size();
}


bam_hdr_t *
HtsParallelReader::get_header() const
{
  if (hts_files.size() == 0)
  {
    BOOST_LOG_TRIVIAL(error) << __HERE__ << " Cannot get header with no opened files.";
    std::exit(1);
  }

  bam_hdr_t * hdr = nullptr;
  std::ostringstream ss;

  {
    bam_hdr_t * curr_hdr = hts_files[0].fp->bam_header;
    // Copy all the text
    ss << std::string(curr_hdr->text, curr_hdr->l_text);
  }

  // Add the @RG lines from other files
  for (long i = 1; i < static_cast<long>(hts_files.size()); ++i)
  {
    bam_hdr_t * curr_hdr = hts_files[i].fp->bam_header;
    std::string const rg = "@RG\t";
    const char * t = curr_hdr->text;
    const char * t_end = curr_hdr->text + (curr_hdr->l_text - 4); // 4 is the size of the header tags

    while (t < t_end)
    {
      long line_size = std::distance(t, std::find(t, t + curr_hdr->l_text, '\n'));

      if (std::equal(rg.begin(), rg.begin() + 4, t))
        ss << std::string(t, line_size) << '\n';

      t += line_size + 1;
    }
  }

  std::string new_header = ss.str();
  hdr = sam_hdr_parse(new_header.size(), new_header.c_str());
  hdr->l_text = new_header.size();
  hdr->text = static_cast<char *>(realloc(hdr->text, sizeof(char) * new_header.size()));
  strncpy(hdr->text, new_header.c_str(), new_header.size());
  return hdr;
}


//***
// Free functions

seqan::BamAlignmentRecord
parse(HtsRecord const & hts_rec)
{
  assert(hts_rec.record);

  seqan::BamAlignmentRecord rec;
  seqan::parse(rec, hts_rec.record);
  return rec;
}


void
get_sequence(seqan::IupacString & seq, seqan::IupacString & rseq, bam1_t const * rec)
{
  auto & core = rec->core;
  assert(core.l_qseq > 0);
  seqan::resize(seq, core.l_qseq);
  seqan::resize(rseq, core.l_qseq);

  auto i = bam_get_seq(rec);

  // Get the sequences
  for (int j = 0; j < core.l_qseq; ++j)
  {
    seq[j] = seq_nt16_str[bam_seqi(i, j)];
    rseq[j] = seq[j];
  }

  seqan::reverseComplement(rseq);
}


void
genotype_only(HtsParallelReader const & hts_preader,
              VcfWriter & writer,
              ReferenceDepth & reference_depth,
              std::vector<TMapGPaths> & maps,
              std::pair<GenotypePaths, GenotypePaths> & prev_paths,
              PHIndex const & ph_index,
              Primers const * primers,
              HtsRecord const & hts_rec,
              seqan::IupacString & seq,
              seqan::IupacString & rseq,
              bool update_prev_paths,
              bool const IS_SV_CALLING)
{
  assert(hts_rec.record);
  assert(hts_rec.file_index >= 0);
  assert(hts_rec.file_index < static_cast<long>(maps.size()));

  long rg_i{0};
  long sample_i{0};
  hts_preader.get_sample_and_rg_index(sample_i, rg_i, hts_rec);

  assert(maps.size() > 0);
  assert(sample_i < static_cast<long>(maps.size()));
  assert(rg_i < static_cast<long>(maps.size()));

  TMapGPaths & map_gpaths = maps[rg_i];
  std::string read_name(reinterpret_cast<char *>(hts_rec.record->data));

  if (update_prev_paths)
  {
    get_sequence(seq, rseq, hts_rec.record);
    prev_paths = align_read(hts_rec.record, seq, rseq, ph_index);
  }

  std::pair<GenotypePaths, GenotypePaths> geno_paths(prev_paths);

  auto find_it = map_gpaths.find(read_name);

  if (find_it == map_gpaths.end())
  {
    if ((hts_rec.record->core.flag & IS_PAIRED) != 0u)
    {
      // Read is in a pair and we need to wait for its mate
      update_paths(geno_paths, hts_rec.record);
      map_gpaths[read_name] = std::move(geno_paths);   // Did not find the read name
    }
    else
    {
      GenotypePaths * selected = update_unpaired_read_paths(geno_paths, hts_rec.record);

      if (selected)
      {
        writer.update_haplotype_scores_geno(*selected, sample_i, primers);
      }
    }
  }
  else
  {
    // Do stuff with paired reads....
    update_paths(geno_paths, hts_rec.record);

    if ((geno_paths.first.flags & IS_FIRST_IN_PAIR) == (find_it->second.first.flags & IS_FIRST_IN_PAIR))
    {
      BOOST_LOG_TRIVIAL(error) << __HERE__ << " Reads with name=" << read_name
                               << " both have IS_FIRST_IN_PAIR="
                               << ((geno_paths.first.flags & IS_FIRST_IN_PAIR) != 0);
      std::exit(1);
    }

    // Find the better pair of geno paths
    std::pair<GenotypePaths *, GenotypePaths *> better_paths =
      get_better_paths(find_it->second, geno_paths);

    if (better_paths.first)
    {
      assert(better_paths.second);

      if (IS_SV_CALLING)
      {
        // Add reference depth
        reference_depth.add_genotype_paths(*better_paths.first, sample_i);
        reference_depth.add_genotype_paths(*better_paths.second, sample_i);
      }

      writer.update_haplotype_scores_geno(better_paths, sample_i, primers);
      //writer.update_haplotype_scores_geno(*better_paths.first, sample_i, primers);
      //writer.update_haplotype_scores_geno(*better_paths.second, sample_i, primers);
    }

    map_gpaths.erase(find_it); // Remove the read name afterwards
  }
}


void
genotype_and_discover(HtsParallelReader const & hts_preader,
                      VcfWriter & writer,
                      ReferenceDepth & reference_depth,
                      VariantMap & varmap,
                      std::vector<TMapGPaths> & maps,
                      std::pair<GenotypePaths, GenotypePaths> & prev_paths,
                      PHIndex const & ph_index,
                      Primers const * primers,
                      HtsRecord const & hts_rec,
                      seqan::IupacString & seq,
                      seqan::IupacString & rseq,
                      bool update_prev_paths)
{
  assert(hts_rec.record);
  assert(hts_rec.file_index >= 0);

  long rg_i = 0;
  long sample_i = 0;
  hts_preader.get_sample_and_rg_index(sample_i, rg_i, hts_rec);

  assert(maps.size() > 0u);
  assert(sample_i < static_cast<long>(maps.size()));
  assert(rg_i < static_cast<long>(maps.size()));

  TMapGPaths & map_gpaths = maps[rg_i];
  std::string read_name(reinterpret_cast<char *>(hts_rec.record->data));

  if (update_prev_paths)
  {
    get_sequence(seq, rseq, hts_rec.record); // Updates seq and rseq
    prev_paths = align_read(hts_rec.record, seq, rseq, ph_index);
  }

  std::pair<GenotypePaths, GenotypePaths> geno_paths(prev_paths);
  auto find_it = map_gpaths.find(read_name);

  if (find_it == map_gpaths.end())
  {
    if (hts_rec.record->core.flag & IS_PAIRED)
    {
      // Read is in a pair and we need to wait for its mate
      update_paths(geno_paths, hts_rec.record);

      further_update_paths_for_discovery(geno_paths, hts_rec.record);
      map_gpaths[read_name] = std::move(geno_paths); // Did not find the read name
    }
    else
    {
      GenotypePaths * selected = update_unpaired_read_paths(geno_paths, hts_rec.record);

      if (selected)
      {
        further_update_unpaired_read_paths_for_discovery(*selected, hts_rec.record);

        // Discover new variants
        {
          std::vector<VariantCandidate> new_vars = selected->find_new_variants();

          if (new_vars.size() > 0)
            varmap.add_variants(std::move(new_vars), sample_i);
        }

        // Add reference depth
        reference_depth.add_genotype_paths(*selected, sample_i);

        // Update haplotype likelihood scores
        writer.update_haplotype_scores_geno(*selected, sample_i, primers);
      }
    }
  }
  else
  {
    // Do stuff with paired reads....
    update_paths(geno_paths, hts_rec.record);

    if ((geno_paths.first.flags & IS_FIRST_IN_PAIR) == (find_it->second.first.flags & IS_FIRST_IN_PAIR))
    {
      BOOST_LOG_TRIVIAL(error) << "Reads with name '" << read_name << "' both have IS_FIRST_IN_PAIR="
                               << ((geno_paths.first.flags & IS_FIRST_IN_PAIR) != 0);
      std::exit(1);
    }

    further_update_paths_for_discovery(geno_paths, hts_rec.record);

    // Find the better pair of geno paths
    std::pair<GenotypePaths *, GenotypePaths *> better_paths =
      get_better_paths(find_it->second, geno_paths);

    if (better_paths.first)
    {
      // Discover new variants
      {
        std::vector<VariantCandidate> new_vars = better_paths.first->find_new_variants();

        if (new_vars.size() > 0)
          varmap.add_variants(std::move(new_vars), sample_i);

        new_vars = better_paths.second->find_new_variants();

        if (new_vars.size() > 0)
          varmap.add_variants(std::move(new_vars), sample_i);
      }

      // Add reference depth
      reference_depth.add_genotype_paths(*better_paths.first, sample_i);
      reference_depth.add_genotype_paths(*better_paths.second, sample_i);

      writer.update_haplotype_scores_geno(*better_paths.first, sample_i, primers);
      writer.update_haplotype_scores_geno(*better_paths.second, sample_i, primers);
    }

    map_gpaths.erase(find_it); // Remove the read name afterwards
  }
}


void
parallel_reader_genotype_only(std::string * out_path,
                              std::vector<std::string> const * hts_paths_ptr,
                              std::vector<double> const * avg_cov_ptr,
                              std::string const * output_dir_ptr,
                              std::string const * reference_fn_ptr,
                              std::string const * region_ptr,
                              PHIndex const * ph_index_ptr,
                              Primers const * primers,

                              std::map<std::pair<uint16_t, uint16_t>,
                                       std::map<std::pair<uint16_t, uint16_t>, int8_t> > *
#ifdef GT_DEV
                              ph_ptr,
#else
                             ,
#endif
                              bool const is_writing_calls_vcf,
                              bool const is_writing_hap)
{
  assert(hts_paths_ptr);
  assert(avg_cov_ptr);
  assert(ph_index_ptr);
  assert(output_dir_ptr);
  assert(reference_fn_ptr);
  assert(region_ptr);

  auto const & hts_paths = *hts_paths_ptr;
  auto const & avg_cov_by_readlen = *avg_cov_ptr;
  auto const & ph_index = *ph_index_ptr;
  auto const & output_dir = *output_dir_ptr;
  auto const & reference = *reference_fn_ptr;
  auto const & region = *region_ptr;

  // Inititalize the HTS parallel reader
  HtsParallelReader hts_preader;
  hts_preader.open(hts_paths, reference, region);

  // Set up VcfWriter
  VcfWriter writer(SPLIT_VAR_THRESHOLD - 1);
  writer.set_samples(hts_preader.get_samples());

  if (writer.pns.size() == 0)
  {
    BOOST_LOG_TRIVIAL(error) << __HERE__ << " No samples were extracted.";
    std::exit(1);
  }

  std::string const & first_sample = writer.pns[0];

  // Set up reference depth tracks and bin counts if we are SV calling
  ReferenceDepth reference_depth;

  if (graph.is_sv_graph)
    reference_depth.set_depth_sizes(writer.pns.size());

  std::vector<TMapGPaths> maps; // One map for each file. Each map relates read names to their graph alignments
  maps.resize(hts_preader.get_num_rg());

#ifndef NDEBUG
  if (Options::const_instance()->stats.size() > 0)
  {
    writer.print_statistics_headers();
    writer.print_variant_details();
    writer.print_variant_group_details();
  }
#endif // ifndef NDEBUG

  long num_records{0};
  long num_duplicated_records{0};
  std::pair<GenotypePaths, GenotypePaths> prev_paths;
  HtsRecord prev;

  // returns true if the read is good enough to attempt graph alignment and use in genotyping
  auto is_good_read =
    [](bam1_t * record) -> bool
    {
      auto const & core = record->core;

      if ((core.flag & IS_UNMAPPED) != 0u)
        return false;

      std::array<char, 16> static constexpr CIGAR_MAP = {{
        'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', 'B', '*', '*', '*', '*', '*', '*'
      }};

      auto cigar_it = record->data + core.l_qname;
      long const N_CIGAR = core.n_cigar;
      bool const is_mate_far_away = core.tid != core.mtid || std::abs(core.pos - core.mpos) > 200000;

      if (core.qual <= 15 && is_mate_far_away)
        return false;

      // Check if both ends are clipped
      if (N_CIGAR >= 2)
      {
        uint32_t opAndCnt;
        memcpy(&opAndCnt, cigar_it, sizeof(uint32_t));
        char cigar_operation_front = CIGAR_MAP[opAndCnt & 15];
        uint32_t count_front = opAndCnt >> 4;

        cigar_it += sizeof(uint32_t) * (N_CIGAR - 1);
        memcpy(&opAndCnt, cigar_it, sizeof(uint32_t));
        char cigar_operation_back = CIGAR_MAP[opAndCnt & 15];
        uint32_t count_back = opAndCnt >> 4;

        bool const is_one_clipped = (cigar_operation_front == 'S' && count_front >= 12) ||
                                    (cigar_operation_back == 'S' && count_back >= 12);

        bool const are_both_clipped = cigar_operation_front == 'S' && cigar_operation_back == 'S';

        if (are_both_clipped || (core.qual <= 15 && is_one_clipped))
          return false;
      }

      return true;
    };

  // Read the first record
  bool is_done = !hts_preader.read_record(prev);
  assert(is_done || prev.record);
  bool const IS_SV_CALLING{graph.is_sv_graph};

  while (!is_done && ((prev.record->core.flag & Options::const_instance()->sam_flag_filter) != 0u ||
                      (IS_SV_CALLING && !is_good_read(prev.record))))
  {
    BOOST_LOG_TRIVIAL(debug) << __HERE__ << " Reread first record";
    is_done = !hts_preader.read_record(prev);
  }

  if (is_done)
  {
    BOOST_LOG_TRIVIAL(debug) << __HERE__ << " No reads found in BAM.";
  }
  else
  {
    // bin counts for the (extreme) coverage filter
    long const first_pos{prev.record->core.pos};

    std::vector<std::vector<uint16_t> > bin_counts;
    bool const IS_COVERAGE_FILTER{IS_SV_CALLING && !Options::const_instance()->no_filter_on_coverage};

    if (IS_COVERAGE_FILTER)
      bin_counts.resize(hts_preader.get_samples().size());

    // Lambda function to update bin counts for the first record. Returns false when there are too many
    // reads in bin
    auto update_bin_count =
      [&](HtsRecord const & hts_rec) -> bool
      {
        if (!IS_COVERAGE_FILTER)
          return true;

        long sample_i{0};
        long rg_i{0};
        hts_preader.get_sample_and_rg_index(sample_i, rg_i, hts_rec);

        if (sample_i >= static_cast<long>(avg_cov_by_readlen.size()) || avg_cov_by_readlen[sample_i] <= 0.0)
          return true;

        uint16_t max_bin_count = std::min(static_cast<long>(std::numeric_limits<uint16_t>::max()),
                                          static_cast<long>(avg_cov_by_readlen[sample_i] * 50.0 * 3.0 + 0.5));

        assert(sample_i < static_cast<long>(bin_counts.size()));
        auto & sample_bin_counts = bin_counts[sample_i];
        long const bin = (hts_rec.record->core.pos - first_pos) / 50l;

        if (bin >= static_cast<long>(sample_bin_counts.size()))
        {
          sample_bin_counts.resize(bin + 1, 0u);
          ++sample_bin_counts[bin]; // increment read count in new bin
          return true;
        }
        else if (sample_bin_counts[bin] > max_bin_count)
        {
          return false;
        }
        else
        {
          ++sample_bin_counts[bin]; // increment read count in bin
          return true;
        }
      };

    update_bin_count(prev);
    ++num_records;
    seqan::IupacString seq;
    seqan::IupacString rseq;

    genotype_only(hts_preader, writer, reference_depth, maps, prev_paths, ph_index, primers,
                  prev, seq, rseq, true /*update prev_geno_paths*/, IS_SV_CALLING);

    HtsRecord curr;

    while (hts_preader.read_record(curr))
    {
      // Ignore reads with the following bits set
      if ((curr.record->core.flag & Options::const_instance()->sam_flag_filter) != 0u ||
          (IS_SV_CALLING && !is_good_read(curr.record)))
      {
        continue;
      }

      ++num_records;

      if (equal_pos_seq(prev.record, curr.record))
      {
        // The two records are equal
        update_bin_count(curr);
        ++num_duplicated_records;

        genotype_only(hts_preader, writer, reference_depth, maps, prev_paths, ph_index, primers,
                      curr, seq, rseq, false /*update prev_geno_paths*/, IS_SV_CALLING);
      }
      else
      {
        if (!update_bin_count(curr))
        {
          --num_records; // Skip record
          continue;
        }

        genotype_only(hts_preader, writer, reference_depth, maps, prev_paths, ph_index, primers,
                      curr, seq, rseq, true /*update prev_geno_paths*/, IS_SV_CALLING);

        hts_preader.move_record(prev, curr); // move curr to prev
      }
    }

    BOOST_LOG_TRIVIAL(debug) << __HERE__ << " Num of duplicated records: "
                             << num_duplicated_records << " / " << num_records;
  }

  // align leftover reads
  if (IS_SV_CALLING)
  {
    auto process_leftovers =
      [&](long sample_i, long rg_i)
      {
        auto & map_gpaths = maps[rg_i];

        for (auto read_it = map_gpaths.begin(); read_it != map_gpaths.end(); ++read_it)
        {
          //BOOST_LOG_TRIVIAL(info) << __HERE__ << " " << read_it->first;
          std::pair<GenotypePaths, GenotypePaths> copy(read_it->second);

          toggle_bit(copy.first.flags, IS_FIRST_IN_PAIR | IS_SEQ_REVERSED);
          toggle_bit(copy.second.flags, IS_FIRST_IN_PAIR | IS_SEQ_REVERSED);

          std::pair<GenotypePaths *, GenotypePaths *> better_paths =
            get_better_paths(read_it->second, copy);

          if (better_paths.first)
          {
            assert(better_paths.second);

            // Add reference depth
            reference_depth.add_genotype_paths(*better_paths.first, sample_i);
            //reference_depth.add_genotype_paths(*better_paths.second, sample_i);

            writer.update_haplotype_scores_geno(*better_paths.first, sample_i, primers);
            //writer.update_haplotype_scores_geno(*better_paths.second, sample_i, primers);
          }
        }
      };

    long const NUM_FILES = hts_preader.hts_files.size();

    for (long file_i{0}; file_i < NUM_FILES; ++file_i)
    {
      auto const & hts_f = hts_preader.hts_files[file_i];

      if (hts_f.rg2sample_i.size() <= 1)
      {
        long const sample_i = hts_f.sample_index_offset;
        long const rg_i = hts_f.rg_index_offset;
        process_leftovers(sample_i, rg_i);
      }
      else
      {
        for (long rg_i{0}; rg_i < static_cast<long>(hts_f.rg2sample_i.size()); ++rg_i)
        {
          long const sample_i = hts_f.sample_index_offset + hts_f.rg2sample_i[rg_i];
          process_leftovers(sample_i, hts_f.rg_index_offset + rg_i);
        }
      }
    }

    maps.clear();
  }

#ifndef NDEBUG
  //if (!gyper::graph.is_sv_graph)
  {
    check_if_maps_are_empty(maps);
  }
#endif // NDEBUG

  // Write haplotype calls
  if (is_writing_hap)
  {
#ifdef GT_DEV
    assert(ph_ptr);
    std::map<std::pair<uint16_t, uint16_t>, std::map<std::pair<uint16_t, uint16_t>, int8_t> > & ph = *ph_ptr;

    for (long ps1{0}; ps1 < (static_cast<long>(writer.haplotypes.size()) - 1l); ++ps1)
    {
      auto const & hap1 = writer.haplotypes[ps1];
      assert(hap1.gts.size() == 1);
      long const order1 = hap1.gts[0].id;

      for (long ps2{ps1 + 1l}; ps2 < static_cast<long>(writer.haplotypes.size()); ++ps2)
      {
        auto const & hap2 = writer.haplotypes[ps2];
        assert(hap2.gts.size() == 1);
        long const order2 = hap2.gts[0].id;

        if (order2 >= order1 + 100)
          break;

        for (long s{0}; s < static_cast<long>(hap1.hap_samples.size()); ++s)
        {
          auto const & hap_samples1 = hap1.hap_samples[s];
          auto const & hap_samples2 = hap2.hap_samples[s];
          auto const & conn_vec = hap_samples1.connections;
          assert(conn_vec.size() == hap1.gts[0].num);

          // skip reference
          for (long cov1{1}; cov1 < static_cast<long>(hap1.gts[0].num); ++cov1)
          {
            // only check for phasing if at least N% of the reads map to this allele
            assert(hap_samples1.gt_coverage.size() == 1);
            assert(cov1 < static_cast<long>(hap_samples1.gt_coverage[0].size()));

            auto const & conn = conn_vec[cov1];
            auto find_it = conn.find(ps2);

            if (find_it == conn.end()) // no info between these haplotypes
              continue;

            bool const is_clearly_seen1 = hap_samples1.gt_coverage[0][cov1] >= 4 ||
                                          (static_cast<double>(hap_samples1.gt_coverage[0][cov1]) /
                                           std::accumulate(hap_samples1.gt_coverage[0].begin(),
                                                           hap_samples1.gt_coverage[0].end(),
                                                           0.0)) >= 0.28;

            bool const is_not_seen1 = hap_samples1.gt_coverage[0][cov1] <= 2 ||
                                      (static_cast<double>(hap_samples1.gt_coverage[0][cov1]) /
                                       std::accumulate(hap_samples1.gt_coverage[0].begin(),
                                                       hap_samples1.gt_coverage[0].end(),
                                                       0.0)) < 0.22;

            auto find1_it = ph.find({ps1, cov1});

            if (find1_it == ph.end())
            {
              auto insert_it = ph.insert({{ps1, cov1}, std::map<std::pair<uint16_t, uint16_t>, int8_t>()});
              assert(insert_it.second);
              find1_it = insert_it.first;
            }

            std::vector<uint16_t> const & support_vec = find_it->second;
            assert(support_vec.size() == hap2.gts[0].num);

            // total support between any alleles in these two variant groups
            long total_support = std::accumulate(support_vec.begin(), support_vec.end(), 0l);

            // skip reference
            for (long cov2{1}; cov2 < static_cast<long>(support_vec.size()); ++cov2)
            {
              double const support = static_cast<double>(support_vec[cov2]);

              //if (total_support <= 2)
              //  continue;

              int8_t is_good{0};

              bool const is_clearly_seen2 = hap_samples2.gt_coverage[0][cov2] >= 4 ||
                                            (static_cast<double>(hap_samples2.gt_coverage[0][cov2]) /
                                             std::accumulate(hap_samples2.gt_coverage[0].begin(),
                                                             hap_samples2.gt_coverage[0].end(),
                                                             0.0)) >= 0.28;

              bool const is_not_seen2 = hap_samples2.gt_coverage[0][cov2] <= 2 ||
                                        (static_cast<double>(hap_samples2.gt_coverage[0][cov2]) /
                                         std::accumulate(hap_samples2.gt_coverage[0].begin(),
                                                         hap_samples2.gt_coverage[0].end(),
                                                         0.0)) < 0.22;


              //if (total_support <= 2)
              //{
              //  // is_good = 0;
              //}
              //else if (hap_samples2.gt_coverage[0][cov2] <= 2 ||
              //         (static_cast<double>(hap_samples2.gt_coverage[0][cov2]) /
              //          std::accumulate(hap_samples2.gt_coverage[0].begin(),
              //                          hap_samples2.gt_coverage[0].end(),
              //                          0.0)) < 0.22)
              //{
              //  is_good = IS_ANY_ANTI_HAP_SUPPORT;
              //}
              //if (!is_clearly_seen1 && !is_clearly_seen2)
              if (is_not_seen1 && is_not_seen2)
              {
                // is_good==0
                continue; // either needs to be clearly seen to determine anything
              }

              if ((is_not_seen1 && is_clearly_seen2) || (is_not_seen2 && is_clearly_seen1))
              {
                // if only one of the alleles is seen, add anti hap support to hinder them from false being in hap
                is_good = IS_ANY_ANTI_HAP_SUPPORT;
              }
              else
              {
                // both alleles are seen, but are they on the same haplotype?
                if (total_support <= 2)
                  continue; // cannot determine - too low support between any two

                if (is_clearly_seen1 && is_clearly_seen2 && support / static_cast<double>(total_support) > 0.78)
                {
                  is_good = IS_ANY_HAP_SUPPORT;
                }
                else if (support / static_cast<double>(total_support) < 0.22)
                {
                  is_good = IS_ANY_ANTI_HAP_SUPPORT;
                }
                else
                {
                  // is_good==0, ambigous
                  continue;
                }
              }

              //auto insert_it = find1_it->second.insert({{ps2, cov2}, 0});
              //insert_it.first->second |= is_good;
              find1_it->second[{ps2, cov2}] |= is_good;

              //if (/*is_good == IS_ANY_ANTI_HAP_SUPPORT &&*/ hap1.gts[0].id >= 86840 && hap1.gts[0].id <= 86850 &&
              //    hap2.gts[0].id <= 86869)
              //{
              //  BOOST_LOG_TRIVIAL(info) << __HERE__
              //                          << " s=" << s
              //                          << " order1=" << hap1.gts[0].id
              //                          << " var1=" << ps1 << "." << cov1
              //                          << " raw_supp1=" << hap_samples1.gt_coverage[0][cov1]
              //                          << " order2=" << hap2.gts[0].id
              //                          << " var2=" << ps2 << "." << cov2
              //                          << " raw_supp2=" << hap_samples2.gt_coverage[0][cov2]
              //                          << " ph_supp=" << support << "/" << total_support
              //                          << " is_good=" << static_cast<long>(is_good);
              //}
            }
          }
        }
      }
    }

    // check ph
    /*
    for (auto const & ph1 : ph)
    {
      std::pair<uint16_t, uint16_t> const & hap1 = ph1.first;

      for (auto const & ph2 : ph1.second)
      {
        std::pair<uint16_t, uint16_t> const & hap2 = ph2.first;

        if (hap1.first == 249 && hap2.first == 250)
        {
          if (ph2.second == IS_ANY_ANTI_HAP_SUPPORT || ph2.second == IS_ANY_HAP_SUPPORT)
          {
            BOOST_LOG_TRIVIAL(info) << __HERE__ << " " << hap1.first << "." << hap1.second
                                    << " o1=" << writer.haplotypes[hap1.first].gts[0].id << " "
                                    << hap2.first << "." << hap2.second
                                    << " o2=" << writer.haplotypes[hap2.first].gts[0].id
                                    << " is_good=" << static_cast<long>(ph2.second);
          }
          else if (ph2.second == (IS_ANY_ANTI_HAP_SUPPORT | IS_ANY_HAP_SUPPORT))
          {
            BOOST_LOG_TRIVIAL(info) << __HERE__ << " " << hap1.first << "." << hap1.second
                                    << " o1=" << writer.haplotypes[hap1.first].gts[0].id << " "
                                    << hap2.first << "." << hap2.second
                                    << " o2=" << writer.haplotypes[hap2.first].gts[0].id
                                    << " is_good=" << static_cast<long>(ph2.second);
          }
        }
      }
    }
    */

    if (!is_writing_calls_vcf)
    {
      assert(!graph.is_sv_graph);

      // No need to create a VCF here if we are also creating a calls VCF
      Vcf vcf;

      // Set sample names
      vcf.sample_names = writer.pns;

      for (long ps{0}; ps < static_cast<long>(writer.haplotypes.size()); ++ps)
      {
        Haplotype & hap = writer.haplotypes[ps];
        vcf.add_haplotype(writer.haplotypes[ps], static_cast<int32_t>(ps));
        hap.clear();
      }

      save_vcf(vcf, output_dir + "/" + first_sample);
    }

#else
    // Write genotype calls. No need to write haplotype calls in SV genotyping
    std::ostringstream hap_calls_path;
    hap_calls_path << output_dir << "/" << first_sample << ".hap";
    HaplotypeCalls hap_calls(writer.get_haplotype_calls());

    BOOST_LOG_TRIVIAL(debug) << __HERE__ << " Writing haplotype calls to '"
                             << hap_calls_path.str() << "'";
    save_calls(hap_calls, hap_calls_path.str());
#endif
  }

  // Optionally we can write _calls.vcf.gz files for each pool of samples
  if (is_writing_calls_vcf)
  {
    BOOST_LOG_TRIVIAL(debug) << __HERE__ << " Writing calls to '"
                             << output_dir << "/"
                             << first_sample << "_*'";

    Vcf vcf;

    // Set sample names
    vcf.sample_names = writer.pns;

    for (long ps{0}; ps < static_cast<long>(writer.haplotypes.size()); ++ps)
    {
      Haplotype & hap = writer.haplotypes[ps];
      vcf.add_haplotype(writer.haplotypes[ps], static_cast<int32_t>(ps));
      hap.clear(); // clear haplotype
    }

    //for (auto & var : vcf.variants)
    //  var.trim_sequences(false);   // Don't keep one match

    if (graph.is_sv_graph)
    {
      reformat_sv_vcf_records(vcf.variants, reference_depth);

      if (vcf.sample_names.size() > 0)
      {
        for (auto & var : vcf.variants)
          var.generate_infos();
      }

      // Non-SVs are at the end from the reformat_sv_vcf_records function, so this is needed
      std::sort(vcf.variants.begin(),
                vcf.variants.end(),
                [](Variant const & a, Variant const & b){
          return a.abs_pos < b.abs_pos || (a.abs_pos == b.abs_pos && a.seqs < b.seqs);
        });
    }

    save_vcf(vcf, output_dir + "/" + first_sample);
  }

  assert(out_path);
  *out_path = output_dir + "/" + first_sample;
}


void
parallel_reader_with_discovery(std::string * out_path,
                               std::vector<std::string> const * hts_paths_ptr,
                               std::string const * output_dir_ptr,
                               std::string const * reference_fn_ptr,
                               std::string const * region_ptr,
                               PHIndex const * ph_index_ptr,
                               Primers const * primers,
                               long const minimum_variant_support,
                               double const minimum_variant_support_ratio,
                               bool const is_writing_calls_vcf,
                               bool const is_writing_hap)
{
  assert(hts_paths_ptr);
  assert(ph_index_ptr);
  assert(output_dir_ptr);
  assert(reference_fn_ptr);
  assert(region_ptr);

  auto const & hts_paths = *hts_paths_ptr;
  auto const & ph_index = *ph_index_ptr;
  auto const & output_dir = *output_dir_ptr;
  auto const & reference = *reference_fn_ptr;
  auto const & region = *region_ptr;

  // Initialize the HTS parallel reader
  HtsParallelReader hts_preader;
  hts_preader.open(hts_paths, reference, region);

  // Set up VcfWriter
  VcfWriter writer(SPLIT_VAR_THRESHOLD - 1);
  writer.set_samples(hts_preader.get_samples());
  std::string const & first_sample = writer.pns[0];

  ReferenceDepth reference_depth;
  reference_depth.set_depth_sizes(writer.pns.size());

  VariantMap varmap;
  varmap.set_samples(hts_preader.get_samples());
  varmap.minimum_variant_support = minimum_variant_support;
  varmap.minimum_variant_support_ratio = minimum_variant_support_ratio;

  std::vector<TMapGPaths> maps; // One map for each file. Each map relates read names to their graph alignments
  maps.resize(hts_preader.get_num_rg());

#ifndef NDEBUG
  if (Options::const_instance()->stats.size() > 0)
  {
    writer.print_statistics_headers();
    writer.print_variant_details();
    writer.print_variant_group_details();
  }
#endif // ifndef NDEBUG

  long num_records{0};
  long num_duplicated_records{0};
  std::pair<GenotypePaths, GenotypePaths> prev_paths;
  HtsRecord prev;

  // Read the first record
  bool is_done = !hts_preader.read_record(prev);
  assert(is_done || prev.record);

  while (!is_done && (prev.record->core.flag & Options::const_instance()->sam_flag_filter) != 0u)
  {
    BOOST_LOG_TRIVIAL(debug) << __HERE__ << " Reread first record";
    is_done = !hts_preader.read_record(prev);
  }

  if (is_done)
  {
    BOOST_LOG_TRIVIAL(debug) << __HERE__ << " No reads found in BAM.";
  }
  else
  {
    ++num_records;
    seqan::IupacString seq;
    seqan::IupacString rseq;
    genotype_and_discover(hts_preader, writer, reference_depth, varmap, maps, prev_paths, ph_index, primers, prev,
                          seq, rseq, true /*update prev_geno_paths*/);
    HtsRecord curr;

    while (hts_preader.read_record(curr))
    {
      // Ignore reads with the following bits set
      if ((curr.record->core.flag & Options::const_instance()->sam_flag_filter) != 0u)
        continue;

      ++num_records;

      if (equal_pos_seq(prev.record, curr.record))
      {
        // The two records are equal
        ++num_duplicated_records;
        genotype_and_discover(hts_preader, writer, reference_depth, varmap, maps, prev_paths, ph_index, primers, curr,
                              seq, rseq, false /*update prev_geno_paths*/);
      }
      else
      {
        genotype_and_discover(hts_preader, writer, reference_depth, varmap, maps, prev_paths, ph_index, primers, curr,
                              seq, rseq, true /*update prev_geno_paths*/);
        hts_preader.move_record(prev, curr); // move curr to prev
      }
    }

    BOOST_LOG_TRIVIAL(debug) << "[graphtyper::hts_parallel_reader] Num of duplicated records: "
                             << num_duplicated_records << " / " << num_records;
  }

#ifndef NDEBUG
  check_if_maps_are_empty(maps);
#endif // NDEBUG

  // Output variants
  varmap.create_varmap_for_all(reference_depth);

#ifndef NDEBUG
  if (Options::instance()->stats.size() > 0)
    varmap.write_stats("2");
#endif // NDEBUG

  std::ostringstream variant_map_path;
  variant_map_path << output_dir << "/" << first_sample << "_variant_map";
  BOOST_LOG_TRIVIAL(debug) << __HERE__ << " Writing variant map to '" << variant_map_path.str();
  save_variant_map(variant_map_path.str(), varmap);

  // Write haplotype calls
  if (is_writing_hap)
  {
    // Write genotype calls. No need to write haplotype calls in SV genotyping
    std::ostringstream hap_calls_path;
    hap_calls_path << output_dir << "/" << writer.pns[0] << ".hap";
    HaplotypeCalls hap_calls(writer.get_haplotype_calls());

    BOOST_LOG_TRIVIAL(debug) << __HERE__ << " Writing haplotype calls to '"
                             << hap_calls_path.str() << "'";
    save_calls(hap_calls, hap_calls_path.str());
  }

  // Optionally we can write _calls.vcf.gz files for each pool of samples
  if (is_writing_calls_vcf)
  {
    BOOST_LOG_TRIVIAL(debug) << __HERE__ << " Writing calls to '"
                             << output_dir << "/"
                             << first_sample << "_*'";

    // Handle SNP/indel graphs
    Vcf vcf;

    // Set sample names
    vcf.sample_names = writer.pns;

    for (long ps = 0; ps < static_cast<long>(writer.haplotypes.size()); ++ps)
    {
      Haplotype & hap = writer.haplotypes[ps];
      vcf.add_haplotype(hap, static_cast<int32_t>(ps));
      hap.clear();
    }

    for (auto & var : vcf.variants)
      var.trim_sequences(false);   // Don't keep one match

    if (graph.is_sv_graph)
    {
      reformat_sv_vcf_records(vcf.variants, reference_depth);

      if (vcf.sample_names.size() > 0)
      {
        for (auto & var : vcf.variants)
          var.generate_infos();
      }

      // Non-SVs are at the end from the reformat_sv_vcf_records function, so this is needed
      std::sort(vcf.variants.begin(),
                vcf.variants.end(),
                [](Variant const & a, Variant const & b){
          return a.abs_pos < b.abs_pos || (a.abs_pos == b.abs_pos && a.seqs < b.seqs);
        });
    }

    save_vcf(vcf, output_dir + "/" + first_sample);
  }

  assert(out_path);
  *out_path = output_dir + "/" + first_sample;
}


void
sam_merge(std::string const & output_sam, std::vector<std::string> const & input_sams)
{
  {
    HtsParallelReader hts_preader;
    hts_preader.open(input_sams);

    HtsRecord hts_rec;
    HtsWriter hts_writer;
    hts_writer.open(output_sam);
    hts_writer.copy_header(hts_preader);

    while (hts_preader.read_record(hts_rec))
    {
      assert(hts_rec.record);
      hts_writer.write(hts_rec.record);
    }

    hts_writer.close();
    hts_preader.close();
  }

  // Remove old files
  for (std::string const & input_sam : input_sams)
  {
    int ret = unlink(input_sam.c_str());

    if (ret < 0)
      BOOST_LOG_TRIVIAL(warning) << __HERE__ << " Unable to remove " << input_sam;
  }
}


} // namespace gyper
