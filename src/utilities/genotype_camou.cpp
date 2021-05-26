#include <fstream>
#include <string>
#include <sstream>
#include <vector>

#include <graphtyper/utilities/logging.hpp>

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


namespace
{

std::string
parse_interval(std::string const & line)
{
  std::string interval;
  int field = 0;

  for (int i {0}; i < static_cast<int>(line.size()); ++i)
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

  print_log(gyper::log_severity::info, "Parsed interval: ", interval);
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

  print_log(log_severity::info, "Path to FASTA reference genome is '", ref_fn, "'");
  print_log(log_severity::info, "Path to interval BED file is '", interval_fn, "'");
  print_log(log_severity::info, "Running with up to ", Options::const_instance()->threads, " threads.");
  print_log(log_severity::info, "Copying data from ", NUM_SAMPLES, " input SAM/BAM/CRAMs to local disk.");

#ifdef GT_DEV
  Options::instance()->is_one_genotype_per_haplotype = true;
#endif // GT_DEV

  // Get intervals
  std::vector<std::string> intervals;

  {
    std::ifstream interval_f(interval_fn);

    if (!interval_f.is_open())
    {
      print_log(log_severity::error, "Could not open BED file ", interval_fn);
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
    print_log(log_severity::error, "Found no intervals in ", interval_fn);
    std::exit(1);
  }

  GenomicRegion genomic_region_bams(intervals[0]);
  --genomic_region_bams.begin; // To make it unique
  print_log(log_severity::info, "Camou genotyping from ", num_intervals, " intervals.");

  Options::instance()->ploidy = 2 * num_intervals;
  std::string tmp_bams = create_temp_dir(genomic_region_bams);

  print_log(log_severity::info, "Temporary folder is ", tmp_bams);


  std::vector<std::string> shrinked_sams;

  if (Options::const_instance()->no_bamshrink)
  {
    shrinked_sams = std::move(sams);
  }
  else
  {
    shrinked_sams = run_bamshrink_multi(sams, ref_fn, interval_fn, avg_cov_by_readlen, tmp_bams);
    std::sort(shrinked_sams.begin(), shrinked_sams.end()); // Sort by input filename
    run_samtools_merge(shrinked_sams, tmp_bams);
  }

  bool constexpr is_writing_hap{false};

  for (auto const & interval : intervals)
  {
    bool is_writing_calls_vcf{false};

    print_log(log_severity::info, "Genotyping interval ", interval);
    GenomicRegion genomic_region(interval);
    std::string tmp = create_temp_dir(genomic_region);

    // Create directories
    mkdir(output_path.c_str(), 0755);
    mkdir((output_path + "/" + genomic_region.chr).c_str(), 0755);

    GenomicRegion padded_genomic_region(genomic_region);
    padded_genomic_region.pad(1000l);

    if (Options::const_instance()->vcf.size() == 0)
    {
      {
        // Iteration 1
        print_log(log_severity::info, "Camou variant discovery step starting.");
        std::string const out_dir = tmp + "/it1";
        std::string const haps_output_vcf = out_dir + "/haps.vcf.gz";
        std::string const discovery_output_vcf = out_dir + "/discovery.vcf.gz";

        mkdir(out_dir.c_str(), 0755);
        construct_graph(ref_fn, "", padded_genomic_region.to_string(), false, false);
        absolute_pos.calculate_offsets(gyper::graph.contigs);
        print_log(log_severity::info, "Graph construction complete.");

#ifndef NDEBUG
        // Save graph in debug mode
        save_graph(out_dir + "/graph");
#endif // NDEBUG

        std::vector<std::string> paths;

        {
          PHIndex ph_index = index_graph(gyper::graph);
          print_log(log_severity::info, "Index construction complete.");
          std::map<std::pair<uint16_t, uint16_t>, std::map<std::pair<uint16_t, uint16_t>, int8_t> > ph;

          paths = gyper::call(shrinked_sams,
                              avg_cov_by_readlen,
                              "", // graph_path
                              ph_index,
                              out_dir,
                              "", // reference
                              ".", // region
                              nullptr,
                              ph,
                              is_writing_calls_vcf,
                              is_writing_hap);
        }

        print_log(log_severity::info, "Variant calling complete.");

        // Append _variant_map
        for (auto & path : paths)
          path += "_variant_map";

        VariantMap varmap;
        varmap.load_many_variant_maps(paths);
        varmap.filter_varmap_for_all();

        Vcf discovery_vcf;
        varmap.get_vcf(discovery_vcf, out_dir + "/final.vcf.gz");
        discovery_vcf.write();
#ifndef NDEBUG
        discovery_vcf.write_tbi_index(); // Write index in debug mode
#endif // NDEBUG

        // free memory
        graph = Graph();
      }

      // Iteration 2
      {
        print_log(log_severity::info, "Starting genotyping step.");

        is_writing_calls_vcf = true;

        std::string const out_dir = tmp + "/it2";
        std::string const haps_output_vcf = out_dir + "/haps.vcf.gz";
        std::string const discovery_output_vcf = out_dir + "/discovery.vcf.gz";
        mkdir(out_dir.c_str(), 0755);

        construct_graph(ref_fn, tmp + "/it1/final.vcf.gz", padded_genomic_region.to_string(), false, false);

    #ifndef NDEBUG
        // Save graph in debug mode
        save_graph(out_dir + "/graph");
    #endif // NDEBUG

        std::vector<std::string> paths;

        {
          PHIndex ph_index = index_graph(gyper::graph);
          std::map<std::pair<uint16_t, uint16_t>, std::map<std::pair<uint16_t, uint16_t>, int8_t> > ph;

          paths = gyper::call(shrinked_sams,
                              avg_cov_by_readlen,
                              "", // graph_path
                              ph_index,
                              out_dir,
                              "", // reference
                              ".", // region
                              nullptr,
                              ph,
                              is_writing_calls_vcf,
                              is_writing_hap);
        }

        bool constexpr filter_zero_qual{true};
        vcf_merge_and_break(paths,
                            tmp + "/graphtyper.vcf.gz",
                            genomic_region.to_string(),
                            filter_zero_qual,
                            false,
                            false);
      }
    }
    else
    {
      print_log(log_severity::info, "Starting genotyping step with input VCF.");

      is_writing_calls_vcf = true;

      std::string const out_dir = tmp + "/it2";
      std::string const haps_output_vcf = out_dir + "/haps.vcf.gz";
      std::string const discovery_output_vcf = out_dir + "/discovery.vcf.gz";
      mkdir(out_dir.c_str(), 0755);

      bool constexpr is_sv_graph{false};
      bool constexpr use_index{true};

      construct_graph(ref_fn,
                      Options::const_instance()->vcf,
                      padded_genomic_region.to_string(),
                      is_sv_graph,
                      use_index);

#ifndef NDEBUG
      // Save graph in debug mode
      save_graph(out_dir + "/graph");
#endif // NDEBUG

      std::vector<std::string> paths;

      {
        PHIndex ph_index = index_graph(gyper::graph);
        std::map<std::pair<uint16_t, uint16_t>, std::map<std::pair<uint16_t, uint16_t>, int8_t> > ph;

        paths = gyper::call(shrinked_sams,
                            avg_cov_by_readlen,
                            "",   // graph_path
                            ph_index,
                            out_dir,
                            "",   // reference
                            ".",   // region
                            nullptr,
                            ph,
                            is_writing_calls_vcf,
                            is_writing_hap);
      }

      bool constexpr filter_zero_qual{true};

      vcf_merge_and_break(paths,
                          tmp + "/graphtyper.vcf.gz",
                          genomic_region.to_string(),
                          filter_zero_qual,
                          false,
                          false);
    }

    auto copy_camou_vcf_to_system =
      [&](std::string const & extension) -> void
      {
        filesystem::path src = tmp + "/graphtyper.vcf.gz" + extension;
        filesystem::path dest = output_path + "/" + genomic_region.to_file_string() + ".vcf.gz" + extension;

        filesystem::copy_file(src, dest, filesystem::copy_options::overwrite_existing);

      };

    copy_camou_vcf_to_system(""); // Copy final camou VCF
    copy_camou_vcf_to_system(".tbi");

    if (!Options::instance()->no_cleanup)
    {
      print_log(log_severity::info, "Cleaning up temporary files for this interval.");
      remove_file_tree(tmp.c_str());
    }
    else
    {
      print_log(log_severity::info, "Temporary files left: ", tmp);
    }

    std::ostringstream ss;

    ss << output_path << "/" << genomic_region.chr << "/"
       << std::setw(9) << std::setfill('0') << (genomic_region.begin + 1)
       << '-'
       << std::setw(9) << std::setfill('0') << genomic_region.end
       << ".vcf.gz";

    graph = Graph();
    print_log(log_severity::info, "Finished ", genomic_region.to_string(), "! Output written at: ", ss.str());
  }

  if (!Options::instance()->no_cleanup)
  {
    print_log(log_severity::info, "Cleaning up temporary files for reads.");
    remove_file_tree(tmp_bams.c_str());
  }
  else
  {
    print_log(log_severity::info, "Temporary files left: ", tmp_bams);
  }

  print_log(log_severity::info, "Finished all ", num_intervals, " intervals.");
}


} // namespace gyper
