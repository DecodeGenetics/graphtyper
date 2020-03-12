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


namespace
{

std::string
parse_interval(std::string const & line)
{
  std::string interval;
  int field = 0;

  for (int i{0}; i < static_cast<int>(line.size()); ++i)
  {
    char c = line[i];

    if (c == '\t')
    {
      if (field == 0)
      {
        c = ':';
        ++field;
      }
      else if (field == 1)
      {
        c = '-';
      }
      else
      {
        return interval;
      }
    }

    interval.append(1, c);
  }

  BOOST_LOG_TRIVIAL(info) << "Parsed interval: " << interval;
  return interval;
}


} // anon namespace


namespace gyper
{

void
genotype_camou(std::string const & interval_fn,
               std::string const & ref_fn,
               std::vector<std::string> const & sams,
               std::string const & output_path,
               std::vector<double> const & avg_cov_by_readlen)
{
  long const NUM_SAMPLES = sams.size();

  BOOST_LOG_TRIVIAL(info) << "Path to FASTA reference genome is '" << ref_fn << "'";
  BOOST_LOG_TRIVIAL(info) << "Path to interval BED file is '" << interval_fn << "'";
  BOOST_LOG_TRIVIAL(info) << "Running with up to " << Options::const_instance()->threads << " threads.";
  BOOST_LOG_TRIVIAL(info) << "Copying data from " << NUM_SAMPLES << " input SAM/BAM/CRAMs to local disk.";

  // Get intervals
  std::vector<std::string> intervals;

  {
    std::ifstream interval_f(interval_fn);

    if (!interval_f.is_open())
    {
      BOOST_LOG_TRIVIAL(error) << "Could not open BED file " << interval_fn;
      std::exit(1);
    }

    for (std::string line; std::getline(interval_f, line);)
    {
      intervals.push_back(parse_interval(line));
    }

    interval_f.close();
  }

  long num_intervals = intervals.size();

  if (num_intervals == 0)
  {
    BOOST_LOG_TRIVIAL(error) << "Found no intervals in " << interval_fn;
    std::exit(1);
  }

  GenomicRegion genomic_region_bams(intervals[0]);
  --genomic_region_bams.begin; // To make it unique
  BOOST_LOG_TRIVIAL(info) << "Camou genotyping from " << num_intervals << " intervals.";

  Options::instance()->ploidy = 2 * num_intervals;
  std::string tmp_bams = create_temp_dir(genomic_region_bams);

  BOOST_LOG_TRIVIAL(info) << "Temporary folder is " << tmp_bams;


  std::vector<std::string> shrinked_sams;

  if (Options::const_instance()->no_bamshrink)
  {
    shrinked_sams = std::move(sams);
  }
  else
  {
    shrinked_sams = run_bamshrink(sams, ref_fn, interval_fn, avg_cov_by_readlen, tmp_bams);
    std::sort(shrinked_sams.begin(), shrinked_sams.end()); // Sort by input filename
    run_samtools_merge(shrinked_sams, tmp_bams);
  }

  bool constexpr is_writing_hap{false};

  for (auto const & interval : intervals)
  {
    bool is_writing_calls_vcf{false};
    bool is_discovery{true};

    BOOST_LOG_TRIVIAL(info) << "Genotyping interval " << interval;
    GenomicRegion genomic_region(interval);
    std::string tmp = create_temp_dir(genomic_region);

    // Create directories
    mkdir(output_path.c_str(), 0755);
    mkdir((output_path + "/" + genomic_region.chr).c_str(), 0755);

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

      long minimum_variant_support = 14;
      double minimum_variant_support_ratio = 0.35 / static_cast<double>(num_intervals);

      std::vector<std::string> paths =
        gyper::call(shrinked_sams,
                    "", // graph_path
                    index_path,
                    out_dir,
                    "", // reference
                    ".",       // region
                    minimum_variant_support,
                    minimum_variant_support_ratio,
                    is_writing_calls_vcf,
                    is_discovery,
                    is_writing_hap);

      BOOST_LOG_TRIVIAL(info) << "Variant calling complete.";

      //Vcf haps_vcf;
      //extract_to_vcf(haps_vcf,
      //               paths,
      //               haps_output_vcf,
      //               false); // is_splitting_vars

      // Append _variant_map
      for (auto & path : paths)
        path += "_variant_map";

      VariantMap varmap;
      varmap.load_many_variant_maps(paths);
      varmap.filter_varmap_for_all();

      Vcf discovery_vcf;
      varmap.get_vcf(discovery_vcf, out_dir + "/final.vcf.gz");
      //std::move(haps_vcf.variants.begin(), haps_vcf.variants.end(), std::back_inserter(discovery_vcf.variants));
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
                    "", // graph_path
                    index_path,
                    out_dir,
                    "", // reference
                    ".",       // region
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
      BOOST_LOG_TRIVIAL(info) << "Cleaning up temporary files for this interval.";
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

    graph = Graph();
    mem_index = MemIndex();
    BOOST_LOG_TRIVIAL(info) << "Finished " << genomic_region.to_string() << "! Output written at: " << ss.str();
  }

  if (!Options::instance()->no_cleanup)
  {
    BOOST_LOG_TRIVIAL(info) << "Cleaning up temporary files for reads.";
    remove_file_tree(tmp_bams.c_str());
  }
  else
  {
    BOOST_LOG_TRIVIAL(info) << "Temporary files left: " << tmp_bams;
  }

  BOOST_LOG_TRIVIAL(info) << "Finished all " << num_intervals << " intervals.";
}


} // namespace gyper
