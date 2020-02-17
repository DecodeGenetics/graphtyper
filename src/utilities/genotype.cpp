#include <graphtyper/graph/absolute_position.hpp>
#include <graphtyper/graph/constructor.hpp>
#include <graphtyper/graph/genomic_region.hpp>
#include <graphtyper/graph/graph_serialization.hpp>
#include <graphtyper/graph/haplotype_extractor.hpp>
#include <graphtyper/index/indexer.hpp>
#include <graphtyper/index/mem_index.hpp>
#include <graphtyper/typer/caller.hpp> // gyper::discover_directly_from_bam
#include <graphtyper/typer/variant_map.hpp>
#include <graphtyper/typer/vcf.hpp>
#include <graphtyper/typer/vcf_operations.hpp>
#include <graphtyper/utilities/bamshrink.hpp>
#include <graphtyper/utilities/genotype.hpp>
#include <graphtyper/utilities/hts_parallel_reader.hpp>
#include <graphtyper/utilities/options.hpp>
#include <graphtyper/utilities/system.hpp>

#include <paw/station.hpp>

#include <boost/log/trivial.hpp>

#include <algorithm>
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <utility>
#include <vector>


namespace gyper
{

std::vector<std::string>
run_bamshrink(std::vector<std::string> const & sams,
              std::string const & ref_fn,
              GenomicRegion const & region,
              std::vector<double> const & avg_cov_by_readlen,
              std::string const & tmp)
{
  // Get SAM/BAM/CRAM files
  create_dir(tmp + "/bams");
  assert(sams.size() == avg_cov_by_readlen.size());
  GenomicRegion bs_region(region); // bs = bamshrink
  bs_region.pad(50);

  paw::Station bamshrink_station(Options::const_instance()->threads);
  std::vector<std::string> output_paths;
  output_paths.reserve(sams.size());

  auto get_basename_wo_ext =
    [](std::string const & sam) -> std::string
    {
      // Get basename without extension
      auto slash_it = std::find(sam.rbegin(), sam.rend(), '/');
      auto dot_it = std::find(sam.rbegin(), sam.rend(), '.');
      assert(slash_it != sam.rend());
      assert(dot_it != sam.rend());
      assert(std::distance(dot_it, slash_it) > 1);
      std::string reversed(dot_it + 1, slash_it);
      return std::string(reversed.rbegin(), reversed.rend());
    };

  for (long s = 0; s < static_cast<long>(sams.size()) - 1l; ++s)
  {
    auto const & sam = sams[s];
    auto const & avg_cov = avg_cov_by_readlen[s];
    std::string basename = get_basename_wo_ext(sam);
    std::ostringstream ss;
    ss << tmp << "/bams/" << basename << ".bam";
    std::string path_out = ss.str();
    output_paths.push_back(path_out);
    bamshrink_station.add_work(bamshrink,
                               bs_region.chr,
                               bs_region.begin,
                               bs_region.end,
                               sam,
                               path_out,
                               avg_cov,
                               ref_fn);
  }

  // Process the last sam on the main thread
  {
    long s = static_cast<long>(sams.size()) - 1l;
    auto const & sam = sams[s];
    auto const & avg_cov = avg_cov_by_readlen[s];
    std::string basename = get_basename_wo_ext(sam);
    std::ostringstream ss;
    ss << tmp << "/bams/" << basename << ".bam";
    std::string path_out = ss.str();
    output_paths.push_back(path_out);

    bamshrink_station.add_to_thread(Options::const_instance()->threads - 1,
                                    bamshrink,
                                    bs_region.chr,
                                    bs_region.begin,
                                    bs_region.end,
                                    sam,
                                    path_out,
                                    avg_cov,
                                    ref_fn);
  }

  std::string thread_info = bamshrink_station.join();
  BOOST_LOG_TRIVIAL(info) << "Finished copying data. Thread work: " << thread_info;
  return output_paths;
}


std::vector<std::string>
run_bamshrink(std::vector<std::string> const & sams,
              std::string const & ref_fn,
              std::string const & interval_fn,
              std::vector<double> const & avg_cov_by_readlen,
              std::string const & tmp)
{
  // Get SAM/BAM/CRAM files
  create_dir(tmp + "/bams");
  assert(sams.size() == avg_cov_by_readlen.size());

  paw::Station bamshrink_station(Options::const_instance()->threads);
  std::vector<std::string> output_paths;
  output_paths.reserve(sams.size());

  auto get_basename_wo_ext =
    [](std::string const & sam) -> std::string
    {
      // Get basename without extension
      auto slash_it = std::find(sam.rbegin(), sam.rend(), '/');
      auto dot_it = std::find(sam.rbegin(), sam.rend(), '.');
      assert(slash_it != sam.rend());
      assert(dot_it != sam.rend());
      assert(std::distance(dot_it, slash_it) > 1);
      std::string reversed(dot_it + 1, slash_it);
      return std::string(reversed.rbegin(), reversed.rend());
    };

  for (long s = 0; s < static_cast<long>(sams.size()) - 1l; ++s)
  {
    auto const & sam = sams[s];
    auto const & avg_cov = avg_cov_by_readlen[s];
    std::string basename = get_basename_wo_ext(sam);
    std::ostringstream ss;
    ss << tmp << "/bams/" << basename << ".bam";
    std::string path_out = ss.str();
    output_paths.push_back(path_out);
    bamshrink_station.add_work(bamshrink_multi, interval_fn, sam, path_out, avg_cov, ref_fn);
  }

  // Process the last sam on the main thread
  {
    long s = static_cast<long>(sams.size()) - 1l;
    auto const & sam = sams[s];
    auto const & avg_cov = avg_cov_by_readlen[s];
    std::string basename = get_basename_wo_ext(sam);
    std::ostringstream ss;
    ss << tmp << "/bams/" << basename << ".bam";
    std::string path_out = ss.str();
    output_paths.push_back(path_out);

    bamshrink_station.add_to_thread(Options::const_instance()->threads - 1,
                                    bamshrink_multi,
                                    interval_fn,
                                    sam,
                                    path_out,
                                    avg_cov,
                                    ref_fn);
  }

  std::string thread_info = bamshrink_station.join();
  BOOST_LOG_TRIVIAL(info) << "Finished copying data. Thread work: " << thread_info;
  return output_paths;
}


void
run_samtools_merge(std::vector<std::string> & shrinked_sams, std::string const & tmp)
{
  if (Options::const_instance()->max_files_open > static_cast<long>(shrinked_sams.size()) &&
      (static_cast<long>(shrinked_sams.size()) / static_cast<long>(Options::const_instance()->threads)) >= 200l)
  {
    BOOST_LOG_TRIVIAL(info) << "Merging input files.";

    long const CHUNK_SIZE =
      std::min(10l, static_cast<long>(shrinked_sams.size() / Options::const_instance()->threads / 100l));

    long const NUM_FILES = static_cast<long>(shrinked_sams.size());
    assert(CHUNK_SIZE > 1);
    std::vector<std::string> new_shrinked_sams;
    std::vector<std::vector<std::string> > all_input_sams;
    all_input_sams.resize(NUM_FILES / CHUNK_SIZE + 1);

    {
      paw::Station merge_station(Options::const_instance()->threads);

      for (long i = 0; (i * CHUNK_SIZE) < NUM_FILES; ++i)
      {
        assert(i < static_cast<long>(all_input_sams.size()));
        std::vector<std::string> & input_sams = all_input_sams[i];
        long const file_i = i * CHUNK_SIZE;
        long const next_file_i = file_i + CHUNK_SIZE;

        if (next_file_i >= NUM_FILES)
        {
          std::copy(shrinked_sams.begin() + file_i, shrinked_sams.end(), std::back_inserter(input_sams));
        }
        else
        {
          std::copy(shrinked_sams.begin() + file_i, shrinked_sams.begin() + next_file_i,
                    std::back_inserter(input_sams));
        }

        if (input_sams.size() == 1)
        {
          // No merging needed
          new_shrinked_sams.push_back(input_sams[0]);
        }
        else if (input_sams.size() > 1)
        {
          std::ostringstream ss;
          ss << tmp << "/bams/merged" << std::setw(5) << std::setfill('0') << i << ".bam";
          new_shrinked_sams.push_back(ss.str());

          if (next_file_i < NUM_FILES)
          {
            merge_station.add_work(sam_merge, new_shrinked_sams[new_shrinked_sams.size() - 1], all_input_sams[i]);
          }
          else
          {
            // Put the very last job to the main thread
            merge_station.add_to_thread(Options::const_instance()->threads - 1,
                                        sam_merge,
                                        new_shrinked_sams[new_shrinked_sams.size() - 1],
                                        all_input_sams[i]);
          }
        }
      }

      std::string thread_info = merge_station.join();
      BOOST_LOG_TRIVIAL(info) << "Finished merging. Thread work: " << thread_info;
    }

#ifndef NDEBUG
    BOOST_LOG_TRIVIAL(debug) << "Number of merged files are " << new_shrinked_sams.size() << "\n";
#endif // NDEBUG
    shrinked_sams = std::move(new_shrinked_sams);
  }
  else
  {
    BOOST_LOG_TRIVIAL(info) << "Skipping merging step. Max files open are " <<
      Options::const_instance()->max_files_open;

    BOOST_LOG_TRIVIAL(info) << "Number of bamShrinked files are " << shrinked_sams.size()
                            << " running accross " << Options::const_instance()->threads << " threads.";
  }
}


void
genotype_with_a_vcf(std::string const & ref_path,
                    std::vector<std::string> const & shrinked_sams,
                    GenomicRegion const & region,
                    GenomicRegion const & padded_region,
                    std::string const & tmp
                    )
{
  // Iteration 1
  BOOST_LOG_TRIVIAL(info) << "Genotyping using an input VCF.";
  std::string const output_vcf = tmp + "/it1/final.vcf.gz";
  std::string const out_dir = tmp + "/it1";
  mkdir(out_dir.c_str(), 0755);
  bool const is_sv_graph{false};
  bool const use_absolute_positions{true};
  bool const check_index{true};
  bool const is_writing_calls_vcf{true};
  bool const is_discovery{false};
  bool const is_writing_hap{false};

  gyper::construct_graph(ref_path,
                         Options::const_instance()->vcf,
                         padded_region.to_string(),
                         is_sv_graph,
                         use_absolute_positions,
                         check_index);

  absolute_pos.calculate_offsets(gyper::graph);

#ifndef NDEBUG
  // Save graph in debug mode
  save_graph(out_dir + "/graph");
#endif // NDEBUG

  std::string const index_path = out_dir + "/graph_gti";
  index_graph(index_path);
  std::vector<std::string> paths = gyper::call(shrinked_sams,
                                               "", // graph_path
                                               index_path,
                                               out_dir,
                                               5, //minimum_variant_support,
                                               0.25, //minimum_variant_support_ratio,
                                               is_writing_calls_vcf,
                                               is_discovery,
                                               is_writing_hap);

  BOOST_LOG_TRIVIAL(info) << "Merging output VCFs.";

  // VCF merge and break_down
  // Append _calls.vcf.gz
  for (auto & path : paths)
    path += "_calls.vcf.gz";

  vcf_merge_and_break(paths, tmp + "/graphtyper.vcf.gz", region.to_string(), true); //> FILTER_ZERO_QUAL

  // free memory
  graph = Graph();
}


void
genotype(std::string ref_path,
         std::vector<std::string> const & sams,
         GenomicRegion const & region,
         std::string const & output_path,
         std::vector<double> const & avg_cov_by_readlen,
         bool const is_copy_reference)
{
  // TODO: If the reference is only Ns then output an empty vcf with the sample names
  // TODO: Extract the reference sequence and use that to discover directly from BAM
  bool is_writing_calls_vcf{true};
  bool is_writing_hap{true};
  bool is_discovery{true};
  long minimum_variant_support = 5;
  double minimum_variant_support_ratio = 0.25;

  long const NUM_SAMPLES = sams.size();
  BOOST_LOG_TRIVIAL(info) << "Genotyping region " << region.to_string();
  BOOST_LOG_TRIVIAL(info) << "Path to genome is '" << ref_path << "'";
  BOOST_LOG_TRIVIAL(info) << "Running with up to " << Options::const_instance()->threads << " threads.";
  BOOST_LOG_TRIVIAL(info) << "Copying data from " << NUM_SAMPLES << " input SAM/BAM/CRAMs to local disk.";

  std::string tmp = create_temp_dir(region);

  BOOST_LOG_TRIVIAL(info) << "Temporary folder is " << tmp;

  // Create directories
  mkdir(output_path.c_str(), 0755);
  mkdir((output_path + "/" + region.chr).c_str(), 0755);
  mkdir((output_path + "/input_sites").c_str(), 0755);
  mkdir((output_path + "/input_sites/" + region.chr).c_str(), 0755);

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

  std::vector<std::string> shrinked_sams;

  if (Options::const_instance()->no_bamshrink)
  {
    shrinked_sams = std::move(sams);
  }
  else
  {
    shrinked_sams = run_bamshrink(sams, ref_path, region, avg_cov_by_readlen, tmp);
    std::sort(shrinked_sams.begin(), shrinked_sams.end()); // Sort by input filename
    run_samtools_merge(shrinked_sams, tmp);
  }

  GenomicRegion padded_region(region);
  padded_region.pad(1000l);

  if (Options::const_instance()->vcf.size() > 0)
  {
    BOOST_LOG_TRIVIAL(info) << "Genotyping a input VCF";
    genotype_with_a_vcf(ref_path, shrinked_sams, region, padded_region, tmp);
  }
  else
  {
    minimum_variant_support = 4;
    minimum_variant_support_ratio = 0.21;

    if (NUM_SAMPLES >= 500)
      minimum_variant_support_ratio += 0.03;

    if (NUM_SAMPLES >= 10000)
    {
      ++minimum_variant_support;
      minimum_variant_support_ratio += 0.03;
    }

#ifdef NDEBUG
    is_writing_calls_vcf = false; // Skip writing calls vcf in release mode in all iterations except the last one
#endif // NDEBUG

    // Iteration 1
    {
      BOOST_LOG_TRIVIAL(info) << "Initial variant discovery step starting.";
      std::string const output_vcf = tmp + "/it1/final.vcf.gz";
      std::string const out_dir = tmp + "/it1";
      mkdir(out_dir.c_str(), 0755);
      gyper::construct_graph(ref_path, "", padded_region.to_string(), false, true, false);
#ifndef NDEBUG
      // Save graph in debug mode
      save_graph(out_dir + "/graph");
#endif // NDEBUG
      absolute_pos.calculate_offsets(gyper::graph);
      auto output_paths = gyper::discover_directly_from_bam("",
                                                            shrinked_sams,
                                                            padded_region.to_string(),
                                                            out_dir,
                                                            minimum_variant_support,
                                                            minimum_variant_support_ratio);
      gyper::VariantMap varmap;
      varmap.load_many_variant_maps(output_paths);
      varmap.filter_varmap_for_all();
      Vcf final_vcf;
      varmap.get_vcf(final_vcf, output_vcf);
      final_vcf.write();
#ifndef NDEBUG
      final_vcf.write_tbi_index(); // Write index in debug mode
#endif // NDEBUG

      // free memory
      graph = Graph();
    }

#ifndef NDEBUG
    std::string stats;

    if (Options::const_instance()->stats.size() > 0)
    {
      stats = std::move(Options::instance()->stats);
      Options::instance()->stats.clear();
    }
#endif // NDEBUG

    // Iteration 2
    //if (false)
    {
      BOOST_LOG_TRIVIAL(info) << "Further variant discovery step starting.";
      std::string const out_dir = tmp + "/it2";
      std::string const index_path = out_dir + "/graph_gti";
      std::string const haps_output_vcf = out_dir + "/haps.vcf.gz";
      std::string const discovery_output_vcf = out_dir + "/discovery.vcf.gz";
      mkdir(out_dir.c_str(), 0755);
      construct_graph(ref_path, tmp + "/it1/final.vcf.gz", padded_region.to_string(), false, true, false);
#ifndef NDEBUG
      // Save graph in debug mode
      save_graph(out_dir + "/graph");
#endif // NDEBUG
      index_graph(index_path);

      minimum_variant_support = 9;
      minimum_variant_support_ratio = 0.32;

      if (NUM_SAMPLES >= 500)
        ++minimum_variant_support;

      if (NUM_SAMPLES >= 10000)
        minimum_variant_support += 2;

      std::vector<std::string> paths =
        gyper::call(shrinked_sams,
                    "", // graph_path
                    index_path,
                    out_dir,
                    minimum_variant_support,
                    minimum_variant_support_ratio,
                    is_writing_calls_vcf,
                    is_discovery,
                    is_writing_hap);

      Vcf haps_vcf;
      extract_to_vcf(haps_vcf,
                     paths,
                     haps_output_vcf,
                     true); // is_splitting_vars

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

    is_discovery = false; // No more discovery
    std::vector<std::string> paths;

    long constexpr FIRST_CALLONLY_ITERATION = 3;
    long constexpr LAST_ITERATION = 4;

    // Iteration FIRST_CALLONLY_ITERATION-LAST_ITERATION
    for (long i = FIRST_CALLONLY_ITERATION; i <= LAST_ITERATION; ++i)
    {
      BOOST_LOG_TRIVIAL(info) << "Call step " << (i - FIRST_CALLONLY_ITERATION + 1) << " starting.";

      if (i == LAST_ITERATION)
      {
        is_writing_calls_vcf = true; // Always write calls vcf in the last iteration
        is_writing_hap = false; // No need for writing .hap
      }

      std::string prev_out_vcf;
      std::string out_dir;

      {
        std::ostringstream ss_prev;
        ss_prev << tmp << "/it" << (i - 1) << "/final.vcf.gz";
        prev_out_vcf = ss_prev.str();

        std::ostringstream ss;
        ss << tmp << "/it" << i;
        out_dir = ss.str();
      }

      mkdir(out_dir.c_str(), 0755);
      std::string const index_path = out_dir + "/graph_gti";
      std::string const haps_output_vcf = out_dir + "/final.vcf.gz";
      construct_graph(ref_path, prev_out_vcf, padded_region.to_string(), false, true, false);

#ifndef NDEBUG
      // Save graph in debug mode
      save_graph(out_dir + "/graph");
#endif // NDEBUG

      index_graph(index_path);
      paths = gyper::call(shrinked_sams,
                          "", // graph_path
                          index_path,
                          out_dir,
                          minimum_variant_support,
                          minimum_variant_support_ratio,
                          is_writing_calls_vcf,
                          is_discovery,
                          is_writing_hap);

      if (i < LAST_ITERATION)
      {
        // Split variants unless its the next-to-last iteration
        bool const is_splitting_vars = (i + 1) < LAST_ITERATION;
        Vcf haps_vcf;
        extract_to_vcf(haps_vcf,
                       paths,
                       haps_output_vcf,
                       is_splitting_vars);

        haps_vcf.write();
#ifndef NDEBUG
        // Write index in debug mode
        haps_vcf.write_tbi_index();
#endif // NDEBUG

        // free memory
        graph = Graph();
        mem_index = MemIndex();
      }
    }

    BOOST_LOG_TRIVIAL(info) << "Merging output VCFs.";

    // VCF merge and break_down
    {
      // Append _calls.vcf.gz
      for (auto & path : paths)
        path += "_calls.vcf.gz";

      vcf_merge_and_break(paths, tmp + "/graphtyper.vcf.gz", region.to_string(), true); //> FILTER_ZERO_QUAL
    }


    BOOST_LOG_TRIVIAL(info) << "Copying results to output directory.";

    // Copy sites to system
    {
      std::ostringstream ss_cmd;
      ss_cmd << "cp -p " << tmp << "/it" << (LAST_ITERATION - 1) << "/final.vcf.gz" << " " // Change to (LAST_ITERATION - 1)
             << output_path << "/input_sites/" << region.chr << "/"
             << std::setw(9) << std::setfill('0') << (region.begin + 1)
             << '-'
             << std::setw(9) << std::setfill('0') << region.end
             << ".vcf.gz";

      int ret = system(ss_cmd.str().c_str());

      if (ret != 0)
      {
        BOOST_LOG_TRIVIAL(error) << "This command failed '" << ss_cmd.str() << "'";
        std::exit(ret);
      }
    }
  }

  // Copy final VCFs
  auto copy_vcf_to_system =
    [&](std::string const & extension) -> void
    {
      std::ostringstream ss_cmd;
      ss_cmd << "cp -p " << tmp << "/graphtyper.vcf.gz" << extension << " "
             << output_path << "/" << region.chr << "/"
             << std::setw(9) << std::setfill('0') << (region.begin + 1)
             << '-'
             << std::setw(9) << std::setfill('0') << region.end
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

  if (!Options::const_instance()->no_cleanup)
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
    ss << output_path << "/" << region.chr << "/"
       << std::setw(9) << std::setfill('0') << (region.begin + 1)
       << '-'
       << std::setw(9) << std::setfill('0') << region.end
       << ".vcf.gz";

    BOOST_LOG_TRIVIAL(info) << "Finished! Output written at: " << ss.str();
  }
}


void
genotype_regions(std::string const & ref_path,
                 std::vector<std::string> const & sams,
                 std::vector<GenomicRegion> const & regions,
                 std::string const & output_path,
                 std::vector<double> const & avg_cov_by_readlen,
                 bool const is_copy_reference)
{
  long const NUM_SAMPLES = sams.size();

  if (NUM_SAMPLES >= 500)
  {
    ++Options::instance()->minimum_extract_variant_support;
    Options::instance()->minimum_extract_score_over_homref += 3;
  }

  if (NUM_SAMPLES >= 10000)
  {
    ++Options::instance()->minimum_extract_variant_support;
    Options::instance()->minimum_extract_score_over_homref += 3;
  }

  //long const NUM_REGIONS = regions.size();
  //gyper::Options & opts = *(gyper::Options::instance());
  //long const NUM_THREADS = opts.threads;

  /*
  ** First make everything in genotype thread safe before enabling this
  if (NUM_REGIONS > NUM_SAMPLES && NUM_THREADS > NUM_SAMPLES)
  {
    // Genotype regions in parallel
    paw::Station region_station(NUM_THREADS);
    opts.threads = 1;
    opts.max_files_open /= NUM_THREADS;

    for (long r = 0; r < static_cast<long>(regions.size()) - 1l; ++r)
    {
      region_station.add_work(genotype, ref_path, sams, regions[r], output_path, avg_cov_by_readlen);
    }

    assert(regions.size() > 0);
    region_station.add_to_thread(NUM_THREADS - 1,
                                 genotype,
                                 ref_path,
                                 sams,
                                 regions[regions.size() - 1l],
                                 output_path,
                                 avg_cov_by_readlen);

  }
  else
  */
  {
    // Genotype regions serially
    for (auto const & region : regions)
    {
      genotype(ref_path,
               sams,
               region,
               output_path,
               avg_cov_by_readlen,
               is_copy_reference);
    }
  }
}


} // namespace gyper
