#pragma once

#include <string>
#include <thread>
#include <vector>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/sinks/text_file_backend.hpp>
#include <boost/log/utility/setup/file.hpp>

#include <graphtyper/constants.hpp>


namespace gyper
{

class Options
{
public:
  // New options
  std::string log;
  std::string output_dir = "results";
  int threads = std::thread::hardware_concurrency();
  bool verbose{false};
  bool vverbose{false};
  bool no_cleanup{false};
  bool no_asterisks{false};
  bool no_decompose{false};
  bool no_bamshrink{false};
  bool no_variant_overlapping{false};
  long ploidy{2};

  /****
   * FILTERING OPTIONS
   */
  bool filter_on_mapq{true};
  bool filter_on_proper_pairs{true};

  /*******************
   * GENERAL OPTIONS *
   *******************/
  std::vector<std::string> regions = {"."}; // "." means the entire SAM file is read.
  std::string stats = ""; // Filename for statistics file

  /************************
   * CONSTRUCTOR OPTIONS *
   ************************/
  std::string vcf = "";
  bool add_all_variants = false;

  /********************
   * INDEXING OPTIONS *
   ********************/
  uint64_t max_index_labels{32};

  /*******************
   * CALLING OPTIONS *
   *******************/
//  std::size_t minimum_variant_support{5}; // 5 default, 3 in eagle calling
//  double minimum_variant_support_ratio{0.25}; // 0.3 default, 0.2 in eagle calling
//  std::size_t certain_variant_support{14};
//  double certain_variant_support_ratio{0.49};
  bool hq_reads{false};
  long max_files_open{1000l}; // Maximum amount of SAM/BAM/CRAM files can be opened at the same time
  long soft_cap_of_variants_in_100_bp_window{22};
  bool get_sample_names_from_filename{false};
  bool output_all_variants{false};
  bool is_one_genotype_per_haplotype{false};
  std::string variant_suffix_id = "";
//  bool is_perfect_alignments_only{false};

  /********************************
   * HAPLOTYPE EXTRACTION OPTIONS *
   ********************************/
  long max_extracted_haplotypes{64};
  int minimum_extract_variant_support{2};
  int minimum_extract_score_over_homref{6};
  /*******************
   * LOGGING OPTIONS *
   *******************/
  boost::shared_ptr<boost::log::sinks::synchronous_sink<boost::log::sinks::text_file_backend> > sink;

  /***********
   * METHODS *
   ***********/
  static Options * instance(); // Gets the only instance (singleton pattern)
  static const Options * const_instance(); // const version of instance()
  void print();


private:
  Options(); // Prevent construction of new instances
  Options(Options const &); // Prevent copy-construction
  Options & operator=(Options const &); // Prevent copy-assignment

  static Options * _instance;
};

} // namespace gyper
