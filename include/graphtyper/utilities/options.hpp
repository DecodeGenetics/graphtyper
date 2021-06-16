#pragma once

#include <string>
#include <thread>
#include <vector>

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
  bool no_sample_name_reordering{false};
  bool no_variant_overlapping{false};
  bool normal_and_no_variant_overlapping{false};
  bool is_all_biallelic{false};
  bool is_only_cigar_discovery{false};
  bool is_discovery_only_for_paired_reads{false};
  bool is_sam_merging_allowed{false};
  long ploidy{2};
  bool is_dropping_genotypes{false};
  uint32_t split_var_threshold{SPLIT_VAR_THRESHOLD};
  bool is_segment_calling{false};
  bool force_ignore_segment{false};
  bool uncompressed_sample_names{false};

  /****
   * FILTERING OPTIONS
   */
  bool filter_on_mapq{true};
  bool filter_on_proper_pairs{true};
  bool filter_on_read_bias{true};
  bool filter_on_strand_bias{true};
  bool no_filter_on_begin_pos{false};
  bool no_filter_on_coverage{false};
  int lr_mapq_filter{5};

  /*******************
   * GENERAL OPTIONS *
   *******************/
  std::vector<std::string> regions = {"."}; // "." means the entire SAM file is read
  std::string stats{};                      // Filename for statistics file

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
  long max_index_labels{75};

  /*******************
   * CALLING OPTIONS *
   *******************/
  bool hq_reads{false};
  bool is_csi{false};
  bool force_align_both_orientations{false};
  int sam_flag_filter{3840};
  long max_files_open{864}; // Maximum amount of SAM/BAM/CRAM files can be opened at the same time
  long soft_cap_of_variants_in_100_bp_window{22};
  bool get_sample_names_from_filename{false};
  bool output_all_variants{false};
  bool is_one_genotype_per_haplotype{false};
  bool force_no_filter_bad_alts{false};
  bool force_no_filter_zero_qual{false};
  std::string variant_suffix_id{};
  std::string primer_bedpe{};
  bool is_extra_call_only_iteration{false};

  // 7 and 0.26 for >= 1000 samples
  long genotype_aln_min_support{4};
  double genotype_aln_min_support_ratio{0.21};
  long genotype_dis_min_support{8};
  double genotype_dis_min_support_ratio{0.30};

  // internal vcf options
  long num_alleles_in_batch{250};

  /********************************
   * HAPLOTYPE EXTRACTION OPTIONS *
   ********************************/
  long max_extracted_haplotypes{100};
  int minimum_extract_variant_support{2};
  int minimum_extract_score_over_homref{27};
  double impurity_threshold{0.15};

  /*******************
   * LOGGING OPTIONS *
   *******************/
  //  boost::shared_ptr<boost::log::sinks::synchronous_sink<boost::log::sinks::text_file_backend> > sink;

  /***********
   * METHODS *
   ***********/
  static Options * instance();             // Gets the only instance (singleton pattern)
  static const Options * const_instance(); // const version of instance()
  void print();

private:
  Options();                            // Prevent construction of new instances
  Options(Options const &);             // Prevent copy-construction
  Options & operator=(Options const &); // Prevent copy-assignment

  static Options * _instance;
};

} // namespace gyper

#define GT_DEV 1
