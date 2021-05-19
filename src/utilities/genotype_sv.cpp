#include <fstream>
#include <string>
#include <sstream>
#include <vector>

#include <graphtyper/graph/absolute_position.hpp>
#include <graphtyper/graph/constructor.hpp>
#include <graphtyper/graph/genomic_region.hpp>
#include <graphtyper/graph/graph_serialization.hpp>
#include <graphtyper/graph/haplotype_extractor.hpp>
#include <graphtyper/index/indexer.hpp>
#include <graphtyper/index/ph_index.hpp>
#include <graphtyper/typer/caller.hpp>
#include <graphtyper/typer/variant_map.hpp>
#include <graphtyper/typer/vcf.hpp>
#include <graphtyper/typer/vcf_operations.hpp>
#include <graphtyper/utilities/filesystem.hpp>
#include <graphtyper/utilities/options.hpp>
#include <graphtyper/utilities/genotype.hpp>
#include <graphtyper/utilities/hts_parallel_reader.hpp>
#include <graphtyper/utilities/system.hpp>

#include <boost/log/trivial.hpp>

namespace gyper
{

void
genotype_sv(std::string ref_path,
            std::string const & sv_vcf,
            std::vector<std::string> const & sams,
            std::vector<double> const & avg_cov_by_readlen,
            GenomicRegion const & genomic_region,
            std::string const & output_path,
            bool const is_copy_reference)
{
  // TODO: If the reference is only Ns then output an empty vcf with the sample names
  // TODO: Extract the reference sequence and use that to discover directly from BAM
  gyper::Options const & copts = *(Options::const_instance());
  bool constexpr is_writing_calls_vcf{true};
  bool constexpr is_writing_hap{false};

  long const NUM_SAMPLES = sams.size();

  BOOST_LOG_TRIVIAL(info) << "SV genotyping region " << genomic_region.to_string();
  BOOST_LOG_TRIVIAL(info) << "Path to genome is '" << ref_path << "'";
  BOOST_LOG_TRIVIAL(info) << "Running with up to " << copts.threads << " threads.";
  BOOST_LOG_TRIVIAL(info) << "Copying data from " << NUM_SAMPLES << " input SAM/BAM/CRAMs to local disk.";
  std::string tmp = create_temp_dir(genomic_region);

  BOOST_LOG_TRIVIAL(info) << "Temporary folder is " << tmp;

  // Create directories
  mkdir(output_path.c_str(), 0755);
  mkdir((output_path + "/" + genomic_region.chr).c_str(), 0755);

  // Copy reference genome to temporary directory
  if (is_copy_reference)
  {
    BOOST_LOG_TRIVIAL(info) << "Copying reference genome FASTA and its index to temporary folder.";

    filesystem::copy_file(ref_path, tmp + "/genome.fa");
    filesystem::copy_file(ref_path + ".fai", tmp + "/genome.fa.fai");

    ref_path = tmp + "/genome.fa";
  }

  //std::vector<std::string> shrinked_sams = std::move(sams);
  GenomicRegion padded_region(genomic_region);
  padded_region.pad_end(200000l);
  GenomicRegion unpadded_region(padded_region);
  padded_region.pad(1000l);

  // Iteration 1 out of 1
  {
    BOOST_LOG_TRIVIAL(info) << "Genotype calling step starting.";
    std::string const output_vcf = tmp + "/it1/final.vcf.gz";
    std::string const out_dir = tmp + "/it1";
    mkdir(out_dir.c_str(), 0755);
    BOOST_LOG_TRIVIAL(info) << "Padded region is: " << padded_region.to_string();

    {
      bool constexpr is_sv_graph{true};
      bool constexpr use_index{true};

      BOOST_LOG_TRIVIAL(info) << "Constructing graph.";

      gyper::construct_graph(ref_path,
                             sv_vcf,
                             padded_region.to_string(),
                             is_sv_graph,
                             use_index);

      BOOST_LOG_TRIVIAL(info) << "Calculating contig offsets.";

      absolute_pos.calculate_offsets(gyper::graph.contigs);
    }

#ifndef NDEBUG
    // Save graph in debug mode
    save_graph(out_dir + "/graph");
#endif // NDEBUG

    PHIndex ph_index = index_graph(gyper::graph);
    std::string reference_fn{}; // empty by default
    std::map<std::pair<uint16_t, uint16_t>, std::map<std::pair<uint16_t, uint16_t>, int8_t> > ph;

    if (Options::const_instance()->force_use_input_ref_for_cram_reading)
      reference_fn = ref_path;

    std::vector<std::string> paths = gyper::call(sams,
                                                 avg_cov_by_readlen,
                                                 "", // graph_path
                                                 ph_index,
                                                 out_dir,
                                                 reference_fn, // reference
                                                 unpadded_region.to_string(), // region
                                                 nullptr,
                                                 ph,
                                                 is_writing_calls_vcf,
                                                 is_writing_hap);

    BOOST_LOG_TRIVIAL(info) << "Merging output VCFs.";

    // VCF merge
    {
      // Append _calls.vcf.gz
      //for (auto & path : paths)
      //  path += "_calls.vcf.gz";

      //> FILTER_ZERO_QUAL, force_no_variant_overlapping
      vcf_merge_and_break(paths, tmp + "/graphtyper.vcf.gz", genomic_region.to_string(), false, false, true);
    }
  }

  // Copy final VCFs
  auto copy_vcf_to_system =
    [&](std::string const & extension) -> void
    {
      filesystem::path src = tmp + "/graphtyper.vcf.gz" + extension;
      filesystem::path dest = output_path + "/" + genomic_region.to_file_string() + ".vcf.gz" + extension;

      filesystem::copy_file(src, dest, filesystem::copy_options::overwrite_existing);
    };

  std::string const index_ext = copts.is_csi ? ".csi" : ".tbi";

  copy_vcf_to_system(""); // Copy final VCF
  copy_vcf_to_system(index_ext); // Copy tabix index for final VCF

  if (!copts.no_cleanup)
  {
    BOOST_LOG_TRIVIAL(info) << "Cleaning up temporary files.";
    remove_file_tree(tmp.c_str());
  }
  else
  {
    BOOST_LOG_TRIVIAL(info) << "Temporary files left: " << tmp;
  }

  {
    std::ostringstream ss;

    ss << output_path << "/" << genomic_region.chr << "/"
       << std::setw(9) << std::setfill('0') << (genomic_region.begin + 1)
       << '-'
       << std::setw(9) << std::setfill('0') << genomic_region.end
       << ".vcf.gz";

    BOOST_LOG_TRIVIAL(info) << "Finished! Output written at: " << ss.str();
  }

  // free memory
  graph = Graph();
}


void
genotype_sv_regions(std::string ref_path,
                    std::string const & sv_vcf,
                    std::vector<std::string> const & sams,
                    std::vector<double> const & avg_cov_by_readlen,
                    std::vector<gyper::GenomicRegion> const & regions,
                    std::string const & output_path,
                    bool const is_copy_reference)
{
  for (auto const & region : regions)
  {
    genotype_sv(ref_path, sv_vcf, sams, avg_cov_by_readlen, region, output_path, is_copy_reference);
  }
}


} // namespace gyper
