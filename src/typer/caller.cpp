#include <cassert> // assert
#include <memory> // std::unique_ptr
#include <sstream> // std::ostringstream
#include <string> // std::string
#include <unordered_map> // std::unordered_map
#include <vector> // std::vector

#include <paw/station.hpp>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/hts_io.h>

#include <boost/log/trivial.hpp> // BOOST_LOG_TRIVIAL

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/absolute_position.hpp> // gyper::absolute_pos
#include <graphtyper/graph/graph_serialization.hpp> // gyper::load_graph
#include <graphtyper/graph/haplotype_calls.hpp>
#include <graphtyper/graph/haplotype_extractor.hpp>
#include <graphtyper/graph/reference_depth.hpp>
#include <graphtyper/index/indexer.hpp> // gyper::index (global)
#include <graphtyper/index/ph_index.hpp>
#include <graphtyper/typer/alignment.hpp>
#include <graphtyper/typer/caller.hpp>
#include <graphtyper/typer/primers.hpp> // gyper::Primers
#include <graphtyper/typer/variant_candidate.hpp>
#include <graphtyper/typer/variant_map.hpp>
#include <graphtyper/utilities/hts_parallel_reader.hpp> // gyper::HtsParallelReader
#include <graphtyper/utilities/hts_reader.hpp> // gyper::HtsReader
#include <graphtyper/utilities/io.hpp>
#include <graphtyper/utilities/options.hpp> // gyper::Options


namespace
{

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

  jobs = Options::const_instance()->threads;
  num_parts = jobs;
  long const MAX_FILES_OPEN = Options::const_instance()->max_files_open;

  if (jobs >= NUM_SAMPLES || jobs >= MAX_FILES_OPEN)
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
     std::string const & graph_path,
     PHIndex const & ph_index,
     std::string const & output_dir,
     std::string const & reference_fn,
     std::string const & region,
     Primers const * primers,
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
  assert(Options::const_instance()->max_files_open > 0);

  long jobs = 1; // Running jobs
  long num_parts = 1; // Number of parts to split the work into
  long const NUM_SAMPLES = hts_paths.size();

  _determine_num_jobs_and_num_parts(jobs, num_parts, NUM_SAMPLES);

  {
    auto it = hts_paths.begin();

    for (long i = 0; i < num_parts; ++i)
    {
      auto end_it = it + NUM_SAMPLES / num_parts + (i < (NUM_SAMPLES % num_parts));
      assert(std::distance(hts_paths.begin(), end_it) <= NUM_SAMPLES);
      spl_hts_paths.emplace_back(new std::vector<std::string>(it, end_it));
      it = end_it;
    }
  }

  BOOST_LOG_TRIVIAL(debug) << "[" << __HERE__ << "] Number of pools = " << spl_hts_paths.size();
  long const NUM_POOLS = spl_hts_paths.size();
  paths.resize(NUM_POOLS);

  // Run in parallel
  {
    paw::Station call_station(jobs); // last parameter is queue_size

    // Push all but the last pool to the thread pool
    if (!is_discovery)
    {
      for (long i = 0; i < (NUM_POOLS - 1l); ++i)
      {
        call_station.add_work(parallel_reader_genotype_only,
                              &paths[i],
                              spl_hts_paths[i].get(),
                              &output_dir,
                              &reference_fn,
                              &region,
                              &ph_index,
                              primers,
                              is_writing_calls_vcf,
                              is_writing_hap);
      }

      // Do the last pool on the current thread
      call_station.add_to_thread(jobs - 1,
                                 parallel_reader_genotype_only,
                                 &paths[NUM_POOLS - 1],
                                 spl_hts_paths[NUM_POOLS - 1].get(),
                                 &output_dir,
                                 &reference_fn,
                                 &region,
                                 &ph_index,
                                 primers,
                                 is_writing_calls_vcf,
                                 is_writing_hap);
    }
    else
    {
      for (long i = 0; i < (NUM_POOLS - 1l); ++i)
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
                                graph.ref_nodes[0].get_label().order :
                                0;

  long read_pos = 0;
  std::string read(seqan::begin(record.seq), seqan::end(record.seq));
  long begin_clipping_size = 0;
  bool is_clipped = false;
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
          BOOST_LOG_TRIVIAL(error) << "[graphtyper::typer::caller] Unable to find RG tag in read.";
          std::exit(1);
        }

        std::string read_group(reinterpret_cast<char *>(rg_tag + 1)); // Skip 'Z'

        auto find_rg_it = rg2sample_i.find(read_group);

        if (find_rg_it == rg2sample_i.end())
        {
          BOOST_LOG_TRIVIAL(error) << "[graphtyper::typer::caller] Unable to find read group. " << read_group;
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

  BOOST_LOG_TRIVIAL(debug) << "[graphtyper::caller] Writing variant map to '"
                           << variant_map_path.str();

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
  BOOST_LOG_TRIVIAL(debug) << "[graphtyper::caller] Number of pools = " << NUM_POOLS;
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


} // namespace gyper
