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
#include <graphtyper/index/mem_index.hpp>
#include <graphtyper/typer/caller.hpp>
#include <graphtyper/typer/variant_map.hpp>
#include <graphtyper/typer/vcf.hpp>
#include <graphtyper/typer/vcf_operations.hpp>
#include <graphtyper/utilities/options.hpp>
#include <graphtyper/utilities/genotype.hpp>
#include <graphtyper/utilities/hts_parallel_reader.hpp>
#include <graphtyper/utilities/system.hpp>


namespace gyper
{

void
genotype_camou(std::string const & interval_fn,
               std::string const & ref_fn,
               std::vector<std::string> const & sams,
               std::string & region,
               std::string const & output_path,
               std::vector<double> const & avg_cov_by_readlen)
{
  long const NUM_SAMPLES = sams.size();
  bool is_writing_calls_vcf{false};
  bool is_writing_hap{false};
  bool is_discovery{true};

  BOOST_LOG_TRIVIAL(info) << "Path to FASTA reference genome is '" << ref_fn << "'";
  BOOST_LOG_TRIVIAL(info) << "Path to interval BED file is '" << interval_fn << "'";
  BOOST_LOG_TRIVIAL(info) << "Running with up to " << Options::const_instance()->threads << " threads.";
  BOOST_LOG_TRIVIAL(info) << "Copying data from " << NUM_SAMPLES << " input SAM/BAM/CRAMs to local disk.";

  // Count the number of intervals
  long num_intervals = 0;

  {
    std::ifstream interval_f(interval_fn);

    if (region.size() == 0)
    {
      // ..then use the first interval from the interval-file
      std::string chrom;
      long begin{0};
      long end{0};

      interval_f >> chrom;
      interval_f >> begin;
      interval_f >> end;

      std::ostringstream ss;
      ss << chrom << ':' << begin << '-' << end;
      region = ss.str();
    }

    for (std::string line; std::getline(interval_f, line);)
      ++num_intervals;

    interval_f.close();
  }

  GenomicRegion genomic_region(region);
  BOOST_LOG_TRIVIAL(info) << "Camou genotyping region " << genomic_region.to_string()
                          << " from " << num_intervals << " intervals.";

  Options::instance()->ploidy = 2 * num_intervals;
  std::string tmp = create_temp_dir(genomic_region);

  BOOST_LOG_TRIVIAL(info) << "Temporary folder is " << tmp;

  // Create directories
  mkdir(output_path.c_str(), 0755);
  mkdir((output_path + "/" + genomic_region.chr).c_str(), 0755);
  std::vector<std::string> shrinked_sams;

  if (Options::const_instance()->no_bamshrink)
  {
    shrinked_sams = std::move(sams);
  }
  else
  {
    shrinked_sams = run_bamshrink(sams, ref_fn, interval_fn, avg_cov_by_readlen, tmp);
    std::sort(shrinked_sams.begin(), shrinked_sams.end()); // Sort by input filename
    run_samtools_merge(shrinked_sams, tmp);
  }

  GenomicRegion padded_genomic_region(genomic_region);
  padded_genomic_region.pad(1000l);

  // Iteration 1
  {
    BOOST_LOG_TRIVIAL(info) << "Camou variant discovery step starting.";
    std::string const out_dir = tmp + "/it1";
    std::string const index_path = out_dir + "/graph_gti";
    std::string const haps_output_vcf = out_dir + "/haps.vcf.gz";
    std::string const discovery_output_vcf = out_dir + "/discovery.vcf.gz";

    mkdir(out_dir.c_str(), 0755);
    construct_graph(ref_fn, "", padded_genomic_region.to_string(), false, true, false);
    absolute_pos.calculate_offsets(gyper::graph);
    BOOST_LOG_TRIVIAL(info) << "Graph construction complete.";

#ifndef NDEBUG
    // Save graph in debug mode
    save_graph(out_dir + "/graph");
#endif // NDEBUG

    index_graph(index_path);
    BOOST_LOG_TRIVIAL(info) << "Index construction complete.";

    long minimum_variant_support = 9;
    double minimum_variant_support_ratio = 0.32 / static_cast<double>(num_intervals);

    if (NUM_SAMPLES >= 500)
      ++minimum_variant_support;

    if (NUM_SAMPLES >= 10000)
      minimum_variant_support += 2;

    std::vector<std::string> paths =
      gyper::call(shrinked_sams,
                  "",   // graph_path
                  index_path,
                  out_dir,
                  minimum_variant_support,
                  minimum_variant_support_ratio,
                  is_writing_calls_vcf,
                  is_discovery,
                  is_writing_hap);

    BOOST_LOG_TRIVIAL(info) << "Variant calling complete.";

    Vcf haps_vcf;
    extract_to_vcf(haps_vcf,
                   paths,
                   haps_output_vcf,
                   false); // is_splitting_vars

    // Append _variant_map
    for (auto & path : paths)
      path += "_variant_map";

    VariantMap varmap;
    varmap.load_many_variant_maps(paths);
    varmap.filter_varmap_for_all();

    Vcf discovery_vcf;
    varmap.get_vcf(discovery_vcf, out_dir + "/final.vcf.gz");
    std::move(haps_vcf.variants.begin(), haps_vcf.variants.end(), std::back_inserter(discovery_vcf.variants));
    discovery_vcf.write();
#ifndef NDEBUG
    discovery_vcf.write_tbi_index(); // Write index in debug mode
#endif // NDEBUG

    // free memory
    graph = Graph();
    mem_index = MemIndex();
  }

  // Iteration 2
  {
    BOOST_LOG_TRIVIAL(info) << "Starting genotyping step.";

    is_writing_calls_vcf = true;
    is_discovery = false;

    std::string const out_dir = tmp + "/it2";
    std::string const index_path = out_dir + "/graph_gti";
    std::string const haps_output_vcf = out_dir + "/haps.vcf.gz";
    std::string const discovery_output_vcf = out_dir + "/discovery.vcf.gz";
    mkdir(out_dir.c_str(), 0755);

    construct_graph(ref_fn, tmp + "/it1/final.vcf.gz", padded_genomic_region.to_string(), false, true, false);

    #ifndef NDEBUG
    // Save graph in debug mode
    save_graph(out_dir + "/graph");
    #endif // NDEBUG

    index_graph(index_path);

    std::vector<std::string> paths =
      gyper::call(shrinked_sams,
                  "",   // graph_path
                  index_path,
                  out_dir,
                  5, // minimum_variant_support - will be ignored
                  0.25, // minimum_variant_support_ratio - will be ignored
                  is_writing_calls_vcf,
                  is_discovery,
                  is_writing_hap);

    for (auto & path : paths)
      path += "_calls.vcf.gz";

    vcf_merge_and_break(paths, tmp + "/graphtyper.vcf.gz", genomic_region.to_string(), false); //> FILTER_ZERO_QUAL
  }

  auto copy_camou_vcf_to_system =
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

  copy_camou_vcf_to_system(""); // Copy final camou VCF
  copy_camou_vcf_to_system(".tbi");

  if (!Options::instance()->no_cleanup)
  {
    BOOST_LOG_TRIVIAL(info) << "Cleaning up temporary files.";
    remove_file_tree(tmp.c_str());
  }
  else
  {
    BOOST_LOG_TRIVIAL(info) << "Temporary files left: " << tmp;
  }

  std::ostringstream ss;

  ss << output_path << "/" << genomic_region.chr << "/"
     << std::setw(9) << std::setfill('0') << (genomic_region.begin + 1)
     << '-'
     << std::setw(9) << std::setfill('0') << genomic_region.end
     << ".vcf.gz";

  BOOST_LOG_TRIVIAL(info) << "Finished! Output written at: " << ss.str();
}


} // namespace gyper
