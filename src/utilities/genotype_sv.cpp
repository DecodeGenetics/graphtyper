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
            GenomicRegion const & genomic_region,
            std::string const & output_path,
            bool const is_copy_reference)
{
  // TODO: If the reference is only Ns then output an empty vcf with the sample names
  // TODO: Extract the reference sequence and use that to discover directly from BAM
  bool constexpr is_writing_calls_vcf {
    true
  };
  bool constexpr is_writing_hap {
    false
  };
  bool constexpr is_discovery {
    false
  };

  long const NUM_SAMPLES = sams.size();

  BOOST_LOG_TRIVIAL(info) << "SV genotyping region " << genomic_region.to_string();
  BOOST_LOG_TRIVIAL(info) << "Path to genome is '" << ref_path << "'";
  BOOST_LOG_TRIVIAL(info) << "Running with up to " << Options::const_instance()->threads << " threads.";
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

    {
      std::ostringstream ss_cmd;
      ss_cmd << "cp " << ref_path << " " << tmp << "/genome.fa";

      int ret = system(ss_cmd.str().c_str());

      if (ret != 0)
      {
        BOOST_LOG_TRIVIAL(error) << "This command failed '" << ss_cmd.str() << "'";
        std::exit(ret);
      }
    }

    // Copy reference genome index
    {
      std::ostringstream ss_cmd;
      ss_cmd << "cp " << ref_path << ".fai " << tmp << "/genome.fa.fai";

      int ret = system(ss_cmd.str().c_str());

      if (ret != 0)
      {
        BOOST_LOG_TRIVIAL(error) << "This command failed '" << ss_cmd.str() << "'";
        std::exit(ret);
      }
    }

    ref_path = tmp + "/genome.fa";
  }

  //std::vector<std::string> shrinked_sams = std::move(sams);
  GenomicRegion padded_region(genomic_region);
  padded_region.pad_end(200000l);
  GenomicRegion unpadded_region(padded_region);
  padded_region.pad(1000l);

  // Iteration 1 out of 1
  {
    BOOST_LOG_TRIVIAL(info) << "Initial variant discovery step starting.";
    std::string const output_vcf = tmp + "/it1/final.vcf.gz";
    std::string const out_dir = tmp + "/it1";
    mkdir(out_dir.c_str(), 0755);
    BOOST_LOG_TRIVIAL(info) << "Padded region is: " << padded_region.to_string();

    {
      bool constexpr is_sv_graph{true};
      bool constexpr use_absolute_positions{true};
      bool constexpr check_index{true};

      BOOST_LOG_TRIVIAL(info) << "Constructing graph.";

      gyper::construct_graph(ref_path,
                             sv_vcf,
                             padded_region.to_string(),
                             is_sv_graph,
                             use_absolute_positions,
                             check_index);

      BOOST_LOG_TRIVIAL(info) << "Calculating contig offsets.";

      absolute_pos.calculate_offsets(gyper::graph.contigs);
    }

#ifndef NDEBUG
    // Save graph in debug mode
    save_graph(out_dir + "/graph");
#endif // NDEBUG

    PHIndex ph_index = index_graph(gyper::graph);

    std::vector<std::string> paths =
      gyper::call(sams,
                  "",   // graph_path
                  ph_index,
                  out_dir,
                  "", // reference
                  unpadded_region.to_string(), // region
                  nullptr,
                  5,//minimum_variant_support,
                  0.25,//minimum_variant_support_ratio,
                  is_writing_calls_vcf,
                  is_discovery,
                  is_writing_hap);

    BOOST_LOG_TRIVIAL(info) << "Merging output VCFs.";

    // VCF merge
    {
      // Append _calls.vcf.gz
      for (auto & path : paths)
        path += "_calls.vcf.gz";

      //> FILTER_ZERO_QUAL, force_no_variant_overlapping
      vcf_merge_and_break(paths, tmp + "/graphtyper.vcf.gz", genomic_region.to_string(), false, false);
    }
  }

  // Copy final VCFs
  auto copy_vcf_to_system =
    [&](std::string const & extension) -> void
    {
      std::ostringstream ss_cmd;
      ss_cmd << "cp -p " << tmp << "/graphtyper.vcf.gz" << extension << " "
             << output_path << "/" << genomic_region.chr << "/"
             << std::setw(9) << std::setfill('0') << (genomic_region.begin + 1)
             << '-'
             << std::setw(9) << std::setfill('0') << genomic_region.end
             << ".vcf.gz" << extension;

      int ret = system(ss_cmd.str().c_str());

      if (ret != 0)
      {
        BOOST_LOG_TRIVIAL(error) << "This command failed '" << ss_cmd.str() << "'";
        std::exit(ret);
      }
    };

  copy_vcf_to_system(""); // Copy final VCF
  copy_vcf_to_system(".tbi"); // Copy tabix index for final VCF

  if (!Options::instance()->no_cleanup)
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
                    std::vector<gyper::GenomicRegion> const & regions,
                    std::string const & output_path,
                    bool const is_copy_reference)
{
  for (auto const & region : regions)
  {
    genotype_sv(ref_path, sv_vcf, sams, region, output_path, is_copy_reference);
  }
}


} // namespace gyper
