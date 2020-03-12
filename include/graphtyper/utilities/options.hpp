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
  std::vector<std::string> regions = {"."}; // "." means the entire SAM file is read
  std::string stats{}; // Filename for statistics file

  /*********************
   * BAMSHRINK OPTIONS *
   *********************/
  int bamshrink_max_fraglen{1000};
  int bamshrink_min_matching{55};
  bool bamshrink_is_not_filtering_mapq0{false};
  int bamshrink_min_readlen{75};
  int bamshrink_min_readlen_low_mapq{94};
  int bamshrink_min_unpair_readlen{94};
  long bamshrink_as_filter_threshold{40};
  bool force_use_input_ref_for_cram_reading{false};

  /************************
   * CONSTRUCTOR OPTIONS *
   ************************/
  std::string vcf{};
  std::string prior_vcf{};
  bool add_all_variants{false};

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
  long max_files_open{1000}; // Maximum amount of SAM/BAM/CRAM files can be opened at the same time
  long soft_cap_of_variants_in_100_bp_window{22};
  bool get_sample_names_from_filename{false};
  bool output_all_variants{false};
  bool is_one_genotype_per_haplotype{false};
  std::string variant_suffix_id = "";

  // 7 and 0.26 for >= 1000 samples
  long genotype_aln_min_support{4};
  double genotype_aln_min_support_ratio{0.21};
  long genotype_dis_min_support{8};
  double genotype_dis_min_support_ratio{0.30};


  /********************************
   * HAPLOTYPE EXTRACTION OPTIONS *
   ********************************/
  long max_extracted_haplotypes{64};
  int minimum_extract_variant_support{2};
  int minimum_extract_score_over_homref{9};
  double impurity_threshold{0.15};
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
