#include <string> // std::string
#include <memory> // std::shared_ptr

#include <stations/join.hpp>
#include <stations/split.hpp>
#include <stations/station.hpp>
#include <stations/worker_queue.hpp>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/hts_io.h>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/log/trivial.hpp>

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/graph_serialization.hpp> // gyper::load_graph
#include <graphtyper/graph/haplotype_calls.hpp>
#include <graphtyper/graph/reference_depth.hpp>
#include <graphtyper/index/rocksdb.hpp>
#include <graphtyper/index/indexer.hpp>
#include <graphtyper/index/mem_index.hpp>
#include <graphtyper/typer/alignment.hpp>
#include <graphtyper/typer/caller.hpp>
#include <graphtyper/typer/discovery.hpp>
#include <graphtyper/typer/graph_swapper.hpp>
#include <graphtyper/typer/genotype_paths.hpp>
#include <graphtyper/typer/vcf_writer.hpp>
#include <graphtyper/typer/segment_calling.hpp>
#include <graphtyper/typer/vcf.hpp>
#include <graphtyper/typer/variant_map.hpp>
#include <graphtyper/utilities/io.hpp>
#include <graphtyper/utilities/kmer_help_functions.hpp>
#include <graphtyper/utilities/type_conversions.hpp>
#include <graphtyper/utilities/hash_seqan.hpp>
#include <graphtyper/utilities/options.hpp>


namespace gyper
{

void
caller(std::shared_ptr<TReads> reads,
       std::shared_ptr<VcfWriter> writer,
       std::shared_ptr<std::size_t const> pn_index
       )
{
  // Check if the reads are paired
  if (reads->size() == 0)
  {
    BOOST_LOG_TRIVIAL(debug) << "[graphtyper::caller] This thread has unexpectedly no reads.";
    return;
  }

  if (seqan::length((*reads)[0].second.seq) == 0)
  {
    // Do not process unpaired reads when calling segments
    if (!Options::instance()->is_segment_calling && !Options::instance()->hq_reads)
    {
      // The reads are unpaired
      std::vector<GenotypePaths> genos;
      align_unpaired_read_pairs(*reads, genos);

      {
        ReferenceDepth reference_depth;

        std::vector<VariantCandidate> variants = discover_variants(genos);
        global_varmap.add_variants(std::move(variants), *pn_index);

        for (auto const & geno : genos)
          reference_depth.add_genotype_paths(geno);

        global_reference_depth.add_reference_depths_from(reference_depth, *pn_index);
      }

      writer->update_haplotype_scores_from_paths(genos, *pn_index);
    }
  }
  else
  {
    std::vector<std::pair<GenotypePaths, GenotypePaths> > geno_pairs = get_best_genotype_paths(*reads);

    if (!Options::instance()->no_new_variants)
    {
      // We only use reference depth for paired reads (since we do not discovery new variants using the unpaired reads)
      ReferenceDepth reference_depth;

      std::vector<VariantCandidate> variants = discover_variants(geno_pairs);
      global_varmap.add_variants(std::move(variants), *pn_index);

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
call(std::vector<std::string> hts_paths,
     std::string graph_path,
     std::string index_path,
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
  std::shared_ptr<Vcf> vcf;
  //Options::instance()->max_index_labels = 2048;

  // Extract sample names from SAM. If that is not possible, try to get a sample name from the filename.
  std::unordered_map<std::string, std::string> rg2sample;
  std::vector<std::string> samples;

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
      BOOST_LOG_TRIVIAL(error) << "[graphtyper::caller] Sorry, currently Graphtyper does not support merged SAM/BAM files. Only one sample per file, please.";
      std::exit(1);
    }

    samples.insert(samples.end(), new_samples.begin(), new_samples.end());

    if (!vcf)
      vcf = std::make_shared<Vcf>(WRITE_BGZF_MODE, output_dir + "/" + samples[0] + "_calls.vcf.gz");

    std::move(new_samples.begin(), new_samples.end(), std::back_inserter(vcf->sample_names));
  }

  assert(samples.size() > 0);
  std::string const & pn = samples[0];
  load_graph(graph_path); // Loads the graph into the global variable 'graph'
  load_index(index_path); // Loads the index into the global variable 'index'
  mem_index.load(); // Loads the in-memory index
  //mem_index.generate_hamming1_hash_map(); // Generate a hashmap with edit distance 1
  index.close(); // Close the RocksDB index, we will use the in-memory index for querying reads

  // Increasing variant distance can increase computational time and file sizes of *.hap files.
  std::shared_ptr<VcfWriter> writer;

  if (Options::instance()->phased_output)
    writer = std::make_shared<VcfWriter>(samples, 60 /*variant distance*/);
  else if (Options::instance()->is_segment_calling)
    writer = std::make_shared<VcfWriter>(samples, 1 /*variant distance*/);
  else
    writer = std::make_shared<VcfWriter>(samples, Options::instance()->max_merge_variant_dist /*variant distance*/);

  std::vector<ReferenceDepth> reference_depth_per_sample;

  {
    std::size_t const READ_BATCH_SIZE = Options::instance()->read_chunk_size;
    std::list<std::shared_ptr<TReads> > reads;
    global_varmap.set_pn_count(samples.size());
    global_reference_depth.set_pn_count(samples.size());
    std::vector<std::shared_ptr<std::size_t> > pn_indexes;

    {
      stations::Station caller_station(Options::instance()->threads, 3 /*queue size*/);

      for (std::size_t i = 0; i < hts_paths.size(); ++i)
      {
        //std::cerr << "PN = " << i << "\n";
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
    BOOST_LOG_TRIVIAL(info) << "[graphtyper::caller] Estimating the likelihoods of few different segments.";
    std::ostringstream segment_calls_path;
    segment_calls_path << output_dir << "/" << pn << "_segments.vcf.gz";
    segment_calling(segment_fasta_files, *writer, segment_calls_path.str(), samples);
  }

  // Generate statistics
  writer->generate_statistics(samples);

  // Write genotype calls
  {
    std::ostringstream hap_calls_path;
    hap_calls_path << output_dir << "/" << pn << ".hap";
    HaplotypeCalls hap_calls(writer->get_haplotype_calls());
    BOOST_LOG_TRIVIAL(info) << "[graphtyper::caller] Writing haplotype calls to '" << hap_calls_path.str() << "'";
    save_calls(hap_calls, hap_calls_path.str());
    BOOST_LOG_TRIVIAL(info) << "[graphtyper::caller] Writing calls to '" << output_dir << "/" << pn << "_calls.vcf.gz'";

    for (unsigned ps = 0; ps < writer->haplotypes.size(); ++ps)
      vcf->add_haplotype(writer->haplotypes[ps], true /*clear haplotypes*/, ps);

    vcf->post_process_variants();
    vcf->write();
  }

  // Write a VCF with all the new variants
  if (!gyper::Options::instance()->no_new_variants)
  {
    std::ostringstream discovery_vcf_path;
    BOOST_LOG_TRIVIAL(info) << "[graphtyper::caller] Writing discovery results to '" << output_dir << "/" << pn << "_variants.vcf.gz'";
    discovery_vcf_path << output_dir << "/" << pn << "_variants.vcf.gz";
    global_varmap.create_varmap_for_all();
    global_varmap.filter_varmap_for_all();
    global_varmap.write_vcf(discovery_vcf_path.str());
  }

  BOOST_LOG_TRIVIAL(info) << "[graphtyper::caller] Finished.";
}


} // namespace gyper
