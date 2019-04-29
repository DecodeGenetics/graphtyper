#include <cassert> // assert
#include <list> // std::list
#include <memory> // std::shared_ptr
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
#include <graphtyper/graph/absolute_position.hpp>
#include <graphtyper/graph/graph_serialization.hpp> // gyper::load_graph
#include <graphtyper/graph/haplotype_calls.hpp>
#include <graphtyper/graph/haplotype_extractor.hpp>
#include <graphtyper/graph/reference_depth.hpp>
#include <graphtyper/index/rocksdb.hpp> // gyper::index (global)
#include <graphtyper/index/indexer.hpp> // gyper::index (global)
#include <graphtyper/index/mem_index.hpp> // gyper::mem_index (global)
#include <graphtyper/typer/alignment.hpp>
#include <graphtyper/typer/caller.hpp>
#include <graphtyper/typer/discovery.hpp>
#include <graphtyper/typer/graph_swapper.hpp>
#include <graphtyper/typer/genotype_paths.hpp> // gyper::GenotypePaths
#include <graphtyper/typer/segment_calling.hpp>
#include <graphtyper/typer/variant_support.hpp>
#include <graphtyper/typer/vcf.hpp>
#include <graphtyper/typer/vcf_writer.hpp> // gyper::VcfWriter
#include <graphtyper/typer/variant_map.hpp>
#include <graphtyper/utilities/io.hpp>
#include <graphtyper/utilities/options.hpp> // gyper::Options


// gyper::TReads is defined in sam_reader.hpp

namespace gyper
{

void
caller(std::shared_ptr<TReads> reads,
       std::shared_ptr<VcfWriter> writer,
       std::shared_ptr<std::size_t const> pn_index
       )
{
  assert(reads->size() > 0);

  // Check if the reads are paired
  if (seqan::length((*reads)[0].second.seq) == 0)
  {
    // Do not process unpaired reads when calling segments
    if (!Options::instance()->is_segment_calling &&
        !Options::instance()->hq_reads
        )
    {
      // The reads are unpaired
      std::vector<GenotypePaths> genos;
      align_unpaired_read_pairs(*reads, genos);

      if (!Options::instance()->no_new_variants)
      {
        std::vector<VariantCandidate> variants = discover_variants(genos);
        global_varmap.add_variants(std::move(variants), *pn_index);
      }

      if (!Options::instance()->no_new_variants || graph.is_sv_graph)
      {
        // Add reference depth
        ReferenceDepth reference_depth;

        for (auto const & geno : genos)
          reference_depth.add_genotype_paths(geno);

        global_reference_depth.add_reference_depths_from(reference_depth, *pn_index);
      }

      writer->update_haplotype_scores_from_paths(genos, *pn_index);
    }
  }
  else
  {
    std::vector<std::pair<GenotypePaths, GenotypePaths> > geno_pairs =
      align_paired_reads(*reads);

    if (!Options::instance()->no_new_variants)
    {
      std::vector<VariantCandidate> variants = discover_variants(geno_pairs);
      global_varmap.add_variants(std::move(variants), *pn_index);
    }

    if (!Options::instance()->no_new_variants || graph.is_sv_graph)
    {
      // Add reference depth
      ReferenceDepth reference_depth;

      for (auto const & geno : geno_pairs)
      {
        reference_depth.add_genotype_paths(geno.first);
        reference_depth.add_genotype_paths(geno.second);
      }

      global_reference_depth.add_reference_depths_from(reference_depth, *pn_index);
    }

    writer->update_haplotype_scores_from_paths(geno_pairs, *pn_index);
  }

  // The reads are no longer needed, let the memory free
  reads->clear();
  reads->shrink_to_fit();
}


void
read_samples(std::unordered_map<std::string, std::string> & rg2sample,
             std::vector<std::string> & samples,
             std::vector<std::string> const & hts_paths
             )
{
  for (auto const & hts_path : hts_paths)
  {
    std::vector<std::string> new_samples;

    if (!Options::instance()->get_sample_names_from_filename)
      new_samples = gyper::get_sample_names_from_bam_header(hts_path, rg2sample);

    if (new_samples.size() == 0)
    {
      std::string sample_name = hts_path.substr(hts_path.rfind('/') + 1, hts_path.rfind('.'));

      if (std::count(sample_name.begin(), sample_name.end(), '.') > 0)
        sample_name = sample_name.substr(0, sample_name.find('.'));

      new_samples.push_back(sample_name);
    }
    else if (new_samples.size() > 1)
    {
      BOOST_LOG_TRIVIAL(error) << "[graphtyper::caller] Sorry, currently Graphtyper does not "
                               << "support merged SAM/BAM files. Only one sample per file, please.";
      std::exit(1);
    }

    std::move(new_samples.begin(), new_samples.end(), std::back_inserter(samples));
  }
}


void
call(std::vector<std::string> const & hts_paths,
     std::string const & graph_path,
     std::string const & index_path,
     std::vector<std::string> const & regions,
     std::vector<std::string> const & segment_fasta_files,
     std::string const & output_dir
     )
{
  if (hts_paths.size() == 0)
  {
    BOOST_LOG_TRIVIAL(error) << "[graphtyper::caller] No input BAM files.";
    std::exit(1);
  }

  assert(regions.size() > 0);

  // Extract sample names from SAM.
  std::unordered_map<std::string, std::string> rg2sample;
  std::vector<std::string> samples;

  // Gather all the sample names
  read_samples(rg2sample, samples, hts_paths);
  assert(samples.size() > 0);
  std::string const & pn = samples[0];
  std::shared_ptr<Vcf> vcf;

  {
    std::ostringstream ss;
    ss << output_dir << "/" << pn << "_calls.vcf.gz";
    vcf = std::make_shared<Vcf>(WRITE_BGZF_MODE, ss.str());

    // Set sample names
    vcf->sample_names = samples;
  }

  load_graph(graph_path); // Loads the graph into the global variable 'graph'
  load_index(index_path); // Loads the index into the global variable 'index'
  mem_index.load(); // Loads the in-memory index
  //mem_index.generate_hamming1_hash_map(); // Generate a hashmap with edit distance 1
  index.close(); // Close the RocksDB index, we will use the in-memory index for querying reads

  // Increasing variant distance can increase computational time and file sizes of *.hap files.
  std::shared_ptr<VcfWriter> writer;

  if (Options::instance()->phased_output)
  {
    writer = std::make_shared<VcfWriter>(samples, 60 /*variant distance*/);
  }
  else if (Options::instance()->is_segment_calling)
  {
    writer = std::make_shared<VcfWriter>(samples, 1 /*variant distance*/);
  }
  else
  {
    writer = std::make_shared<VcfWriter>(samples,
                                         Options::instance()->max_merge_variant_dist
                                         /*variant distance*/);
  }

  if (Options::instance()->stats.size() > 0)
  {
    writer->print_statistics_headers();
    writer->print_variant_details();
    writer->print_variant_group_details();
  }

  {
    std::size_t const READ_BATCH_SIZE = Options::instance()->read_chunk_size;
    std::list<std::shared_ptr<TReads> > reads;

    if (!Options::instance()->no_new_variants)
    {
      // Only use varmap when discovering new variants
      global_varmap.set_samples(samples);
    }

    if (!Options::instance()->no_new_variants || graph.is_sv_graph)
      global_reference_depth.set_pn_count(samples.size());

    std::vector<std::shared_ptr<std::size_t> > pn_indexes;

    {
      paw::Station caller_station(Options::instance()->threads, 3 /*queue size*/);
      caller_station.options.verbosity = 2; // Print messages

      for (long i = 0; i < static_cast<long>(hts_paths.size()); ++i)
      {
        SamReader sam_reader(hts_paths[i], regions);

        // Flush logs
        if (Options::instance()->sink)
          Options::instance()->sink->flush();

        while (true)
        {
          pn_indexes.push_back(std::make_shared<std::size_t>(i));
          reads.push_back(std::make_shared<TReads>(sam_reader.read_N_reads(READ_BATCH_SIZE)));

          if (reads.back()->size() > 0)
          {
            caller_station.add(caller, reads.back(), writer, pn_indexes.back());
          }
          else
          {
            // Remove the batch with no reads
            reads.pop_back();
            pn_indexes.pop_back();
            break;
          }
        }
      }

      caller_station.join();
    }
  }

  // Check segments
  if (segment_fasta_files.size() > 0)
  {
    BOOST_LOG_TRIVIAL(info) << "[graphtyper::caller]"
                            << "Estimating the likelihoods of few different segments.";

    std::ostringstream segment_calls_path;
    segment_calls_path << output_dir << "/" << pn << "_segments.vcf.gz";
    segment_calling(segment_fasta_files, *writer, segment_calls_path.str(), samples);
  }

  // Write genotype calls
  if (!graph.is_sv_graph)
  {
    // No need to write haplotype calls in SV genotyping
    std::ostringstream hap_calls_path;
    hap_calls_path << output_dir << "/" << pn << ".hap";
    HaplotypeCalls hap_calls(writer->get_haplotype_calls());
    BOOST_LOG_TRIVIAL(info) << "[graphtyper::caller] Writing haplotype calls to '"
                            << hap_calls_path.str() << "'";
    save_calls(hap_calls, hap_calls_path.str());
    BOOST_LOG_TRIVIAL(info) << "[graphtyper::caller] Writing calls to '"
                            << output_dir
                            << "/"
                            << pn
                            << "_calls.vcf.gz'";
  }


  for (long ps = 0; ps < static_cast<long>(writer->haplotypes.size()); ++ps)
  {
    vcf->add_haplotype(writer->haplotypes[ps],
                       true /*clear haplotypes*/,
                       static_cast<uint32_t>(ps)
                       );
  }

  vcf->post_process_variants(false /*normalize variants?*/, true /*trim variant sequences?*/);
  vcf->write();

  // Write a VCF with all the new variants
  if (!gyper::Options::instance()->no_new_variants)
  {
    std::ostringstream discovery_vcf_path;
    discovery_vcf_path << output_dir << "/" << pn << "_variants.vcf.gz";

    BOOST_LOG_TRIVIAL(info) << "[graphtyper::caller] Writing discovery results to '"
                            << discovery_vcf_path.str();

    global_varmap.create_varmap_for_all();

    std::ostringstream variant_map_path;
    variant_map_path << output_dir << "/" << pn << "_variant_map";

    BOOST_LOG_TRIVIAL(info) << "[graphtyper::caller] Writing variant map to '"
                            << variant_map_path.str();

    save_variant_map(variant_map_path.str(), global_varmap);
    global_varmap.filter_varmap_for_all();
    global_varmap.write_vcf(discovery_vcf_path.str());
  }

  BOOST_LOG_TRIVIAL(info) << "[graphtyper::caller] Finished.";
}


std::vector<VariantCandidate>
find_variants_in_cigar(seqan::BamAlignmentRecord const & record,
                       GenomicRegion const & region,
                       std::string const & ref
                       )
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
    return new_var_candidates;

  long ref_to_seq_offset = 0;

  // Reset to original ref_abs_pos
  ref_abs_pos = absolute_pos.get_absolute_position(region.chr, record.beginPos + 1);

  Variant new_var =
    make_variant_of_gapped_strings(reference_seq, read_seq, ref_abs_pos, ref_to_seq_offset);

  if (new_var.abs_pos == 0)
    return new_var_candidates;

  std::vector<Variant> new_vars =
    extract_sequences_from_aligned_variant(std::move(new_var), SPLIT_VAR_THRESHOLD);

  // Too many candidates is a smell of a problem
  if (new_vars.size() >= 5)
    return new_var_candidates;

  new_var_candidates.resize(new_vars.size());

  for (unsigned i = 0; i < new_vars.size(); ++i)
  {
    assert(i < new_var_candidates.size());
    auto & new_var = new_vars[i];
    auto & new_var_candidate = new_var_candidates[i];
    assert(new_var.seqs.size() >= 2);
    new_var.normalize();

    //// Check if high or low quality
    long const r = new_var.abs_pos - ref_to_seq_offset;
    long const r_end = r + new_var.seqs[1].size();
    ref_to_seq_offset += new_var.seqs[0].size() - new_var.seqs[1].size();

    //// In case no qual is available
    if (r_end < static_cast<long>(seqan::length(record.qual)))
    {
      int const MAX_QUAL = *std::max_element(seqan::begin(record.qual) + r,
                                             seqan::begin(record.qual) + r_end
                                             ) - 33;

      new_var_candidate.is_low_qual = MAX_QUAL < 25;
    }

    new_var_candidate.abs_pos = new_var.abs_pos;
    new_var_candidate.seqs = std::move(new_var.seqs);
    new_var_candidate.is_low_qual = false;
    new_var_candidate.is_in_proper_pair = seqan::hasFlagAllProper(record);
    new_var_candidate.is_mapq0 = record.mapQ == 0;
    new_var_candidate.is_unaligned = record.beginPos == -1;
    new_var_candidate.is_clipped = is_clipped;
    new_var_candidate.is_first_in_pair = seqan::hasFlagFirst(record);
    new_var_candidate.is_seq_reversed = seqan::hasFlagRC(record);
    new_var_candidate.original_pos =
      static_cast<uint32_t>(record.beginPos + 1 - begin_clipping_size);
  }

  return new_var_candidates;
}


void
discover_directly_from_bam(std::string const & graph_path,
                           std::vector<std::string> const & sams,
                           std::vector<std::string> const & regions,
                           std::string const & output_dir
                           )
{
  load_graph(graph_path);
  assert(regions.size() > 0);
  assert(regions.size() == 1); // Only support for one region right now

  GenomicRegion region(regions[0]);
  // auto const begin_pos = region.begin + 1; // Change to 1-based
  // auto const end_pos = region.end + 1;
  std::size_t const REGION_SIZE = region.end - region.begin;
  graph.generate_reference_genome();
  std::string ref_str(graph.reference.begin(), graph.reference.end());

  // Extract sample names from SAM
  std::unordered_map<std::string, std::string> rg2sample; // Note: Not used yet
  std::vector<std::string> samples;

  // Gather all the sample names
  read_samples(rg2sample, samples, sams);
  assert(samples.size() > 0);
  std::string const & pn = samples[0];

  global_varmap.set_samples(samples);
  global_reference_depth.set_pn_count(samples.size());

  for (long s = 0; s < static_cast<long>(sams.size()); ++s)
  {
    auto const & sam = sams[s];
    seqan::HtsFile hts(sam.c_str(), "r");
    seqan::BamAlignmentRecord record;
    gyper::ReferenceDepth ref_depth;
    ref_depth.depth.resize(REGION_SIZE);

    while (seqan::readRecord(record, hts))
    {
      // Skip supplementary, secondary and QC fail, duplicated and unmapped reads
      if (seqan::hasFlagSupplementary(record) ||
          seqan::hasFlagSecondary(record) ||
          seqan::hasFlagQCNoPass(record) ||
          seqan::hasFlagDuplicate(record) ||
          seqan::hasFlagUnmapped(record)
          )
      {
        continue;
      }

      // Skip reads that are clipped on both ends
      if (length(record.cigar) == 0 ||
          (record.cigar[0].operation == 'S' &&
           record.cigar[length(record.cigar) - 1].operation == 'S'))
      {
        continue;
      }

      assert(record.beginPos >= 0);

      // 0-based positions
      int64_t begin_pos = absolute_pos.get_absolute_position(region.chr, record.beginPos + 1);
      int64_t end_pos = begin_pos + seqan::getAlignmentLengthInRef(record);

      // Check if read is within region
      if (begin_pos < region.get_absolute_begin_position() ||
          end_pos > region.get_absolute_end_position()
          )
      {
        continue;
      }

      assert(begin_pos >= 0);
      assert(end_pos >= 0);

      // Update reference depth (very arbitrarily)
      if (end_pos - begin_pos < 50)
        ref_depth.add_depth(begin_pos, end_pos);
      else
        ref_depth.add_depth(begin_pos + 4, end_pos - 4);

      // Add variant candidates
      std::vector<VariantCandidate> var_candidates =
        find_variants_in_cigar(record, region, ref_str);

      global_varmap.add_variants(std::move(var_candidates), s);
    }

    global_reference_depth.add_reference_depths_from(ref_depth, s);
  }

  // Write variant map
  {
    std::ostringstream variant_map_path;
    variant_map_path << output_dir << "/" << pn << "_variant_map";

    BOOST_LOG_TRIVIAL(info) << "[graphtyper::caller] Writing variant map to '"
                            << variant_map_path.str();

    global_varmap.create_varmap_for_all();
    save_variant_map(variant_map_path.str(), global_varmap);

    // Uncomment for debbuging purposes
    //global_varmap.filter_varmap_for_all();
    //variant_map_path << ".vcf.gz";
    //global_varmap.write_vcf(variant_map_path.str());
  }
}


} // namespace gyper
