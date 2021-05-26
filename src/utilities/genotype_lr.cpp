#include <graphtyper/graph/absolute_position.hpp>
#include <graphtyper/graph/constructor.hpp>
#include <graphtyper/graph/genomic_region.hpp>
#include <graphtyper/graph/graph_serialization.hpp>
#include <graphtyper/graph/haplotype_extractor.hpp>
#include <graphtyper/index/indexer.hpp>
#include <graphtyper/index/ph_index.hpp>
#include <graphtyper/typer/caller.hpp> // gyper::discover_directly_from_bam
#include <graphtyper/typer/primers.hpp>
#include <graphtyper/typer/variant_map.hpp>
#include <graphtyper/typer/vcf.hpp>
#include <graphtyper/typer/vcf_operations.hpp>
#include <graphtyper/utilities/bamshrink.hpp>
#include <graphtyper/utilities/filesystem.hpp>
#include <graphtyper/utilities/genotype.hpp>
#include <graphtyper/utilities/hts_parallel_reader.hpp>
#include <graphtyper/utilities/options.hpp>
#include <graphtyper/utilities/system.hpp>

#include <paw/station.hpp>

#include <graphtyper/utilities/logging.hpp>

#include <algorithm>
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <utility>
#include <vector>


namespace gyper
{

void
genotype_lr(std::string ref_path,
            std::vector<std::string> const & sams,
            GenomicRegion const & region,
            std::string const & output_path,
            bool const is_copy_reference)
{
  //bool is_writing_calls_vcf{true};
  //bool is_writing_hap{true};
  //bool is_discovery{true};

  //long minimum_variant_support = 5;
  //double minimum_variant_support_ratio = 0.25;
  gyper::Options const & copts = *(Options::const_instance());

  long const NUM_SAMPLES = sams.size();
  print_log(log_severity::info, "Genotyping region ", region.to_string());
  print_log(log_severity::info, "Path to genome is '", ref_path, "'");
  print_log(log_severity::info, "Running with up to ", copts.threads, " threads.");
  print_log(log_severity::info, "Copying data from ", NUM_SAMPLES, " input SAM/BAM/CRAMs to local disk.");

  std::string tmp = create_temp_dir(region);

  print_log(log_severity::info, "Temporary folder is ", tmp);

  // Create directories
  mkdir(output_path.c_str(), 0755);
  mkdir((output_path + "/" + region.chr).c_str(), 0755);
  mkdir((output_path + "/input_sites").c_str(), 0755);
  mkdir((output_path + "/input_sites/" + region.chr).c_str(), 0755);

  // Copy reference genome to temporary directory
  if (is_copy_reference)
  {
    print_log(log_severity::info, "Copying reference genome FASTA and its index to temporary folder.");

    filesystem::copy_file(ref_path, tmp + "/genome.fa");
    filesystem::copy_file(ref_path + ".fai", tmp + "/genome.fa.fai");

    ref_path = tmp + "/genome.fa";
  }

  GenomicRegion padded_region(region);
  //padded_region.pad(1000l);

  //if (copts.vcf.size() > 0)
  //{
  //  BOOST_LOG_TRIVIAL(info) << "Genotyping a input VCF";
  //  genotype_only_with_a_vcf(ref_path, shrinked_sams, region, padded_region, primers.get(), tmp);
  //}
  //else
  //std::vector<double> avg_cov_by_readlen(sams.size(), -1.0);

  // Iteration 1
  {
    print_log(log_severity::info, "Initial variant discovery step starting.");
    std::string const output_vcf = tmp + "/graphtyper.vcf.gz";
    std::string const out_dir = tmp + "/it1";
    mkdir(out_dir.c_str(), 0755);

    gyper::construct_graph(ref_path, "", padded_region.to_string(), false, false);

    Vcf variant_sites(WRITE_BGZF_MODE, output_vcf);
    gyper::streamlined_lr_genotyping(sams,
                                     ref_path,
                                     padded_region.to_string(),
                                     variant_sites);
    //variant_sites.gene
    variant_sites.write(".", copts.threads);
    variant_sites.write_tbi_index();
  }

  // Output the results
  std::string const index_ext = copts.is_csi ? ".vcf.gz.csi" : ".vcf.gz.tbi";

  // Copy final VCFs
  auto copy_to_results =
    [&](std::string const & basename_no_ext, std::string const & extension, std::string const & id)
    {
      filesystem::path src = tmp + "/" + basename_no_ext + extension;
      filesystem::path dest = output_path + "/" + region.to_file_string() + id + extension;

      filesystem::copy_file(src, dest, filesystem::copy_options::overwrite_existing);
    };

  std::string basename_no_ext{"graphtyper"};

  // Check if tabix file exists
  if (!is_file(tmp + "/graphtyper.vcf.gz.tbi") && !is_file(tmp + "/graphtyper.vcf.gz.csi"))
  {
    print_log(log_severity::warning, "Tabix creation appears to have failed, "
                              , "I will retry sorting the VCF by reading it in whole.");

    bool const no_sort{false};
    bool const sites_only{false};
    bool const write_tbi{true};

    std::string region{"."}; // already extracted region

    vcf_concatenate({tmp + "/graphtyper.vcf.gz"},
                    tmp + "/graphtyper.sorted.vcf.gz",
                    no_sort,
                    sites_only,
                    write_tbi,
                    region);

    basename_no_ext = "graphtyper.sorted";
  }

  copy_to_results(basename_no_ext, ".vcf.gz", ""); // Copy final VCF
  copy_to_results(basename_no_ext, index_ext, ""); // Copy tabix index for final VCF

  if (copts.normal_and_no_variant_overlapping)
  {
    copy_to_results("graphtyper.no_variant_overlapping", ".vcf.gz", ".no_variant_overlapping");
    copy_to_results("graphtyper.no_variant_overlapping", index_ext, ".no_variant_overlapping");
  }

  if (!copts.no_cleanup)
  {
    print_log(log_severity::info, "Cleaning up temporary files.");
    remove_file_tree(tmp.c_str());
  }
  else
  {
    print_log(log_severity::info, "Temporary files left: ", tmp);
  }

  {
    std::ostringstream ss;
    ss << output_path << "/" << region.chr << "/"
       << std::setw(9) << std::setfill('0') << (region.begin + 1)
       << '-'
       << std::setw(9) << std::setfill('0') << region.end
       << ".vcf.gz";

    print_log(log_severity::info, "Finished! Output written at: ", ss.str());
  }
}


void
genotype_lr_regions(std::string ref_path,
                    std::vector<std::string> const & sams,
                    std::vector<gyper::GenomicRegion> const & regions,
                    std::string const & output_path,
                    bool const is_copy_reference)
{
  for (auto const & region : regions)
  {
    genotype_lr(ref_path, sams, region, output_path, is_copy_reference);
  }
}


} // namespace gyper
