#pragma once

#include <string>
#include <vector>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/sinks/text_file_backend.hpp>
#include <boost/log/utility/setup/file.hpp>

#include <graphtyper/constants.hpp>


namespace gyper
{

enum Verbosity {ERROR = 0, WARNING = 1, INFO = 2};

class Options
{
public:
  /*******************
   * GENERAL OPTIONS *
   *******************/
  Verbosity verbosity = ERROR; // How verbose logging should be
  unsigned threads = 1; // How many threads should be used by Graphtyper (note that RocksDB may also use some additional threads)
  std::vector<std::string> regions = {"."}; // "." means the entire SAM file is read.
  std::string stats = ""; // Filename for statistics file

  /************************
   * CONSTRUCTOR OPTIONS *
   ************************/
  std::string fasta = "";
  std::string vcf = "";
  bool add_all_variants = false;

  /********************
   * INDEXING OPTIONS *
   ********************/
  uint64_t max_index_labels = 32;

  /*******************
   * CALLING OPTIONS *
   *******************/
  std::vector<std::string> sam;
  uint16_t epsilon_0_exponent = 13;
  std::size_t minimum_variant_support = 5;
  double minimum_variant_support_ratio = 0.25;
  std::size_t certain_variant_support = 14;
  double certain_variant_support_ratio = 0.49;
  bool no_new_variants = false;
  bool hq_reads = false;
  std::size_t read_chunk_size = 16384ul; // 2^14
  bool phased_output = false;
  std::size_t soft_cap_of_variants_in_100_bp_window = 11;
  std::size_t soft_cap_of_non_snps_in_100_bp_window = 7;
  std::size_t hard_cap_of_variants_in_100_bp_window = 18;
  std::size_t hard_cap_of_non_snps_in_100_bp_window = 11;
  bool is_segment_calling = false;
  uint16_t max_merge_variant_dist = SPLIT_VAR_THRESHOLD + 1;
  bool get_sample_names_from_filename = false;
  bool output_all_variants = false;
  bool always_query_hamming_distance_one = false;
  bool is_one_genotype_per_haplotype = false;
  std::string variant_suffix_id = "";
  bool is_perfect_alignments_only = false;

  // Insert size options
  uint32_t max_insert_size = 1000;
  uint32_t optimal_insert_size = 300;
  uint32_t max_insert_size_threshold = 150; // Allowed threshold from the optimal insert size to be considered HQ

  // Alignment constraints
  uint32_t MAX_SEED_NUMBER_ALLOWING_MISMATCHES = 64;
  uint32_t MAX_SEED_NUMBER_FOR_WALKING = 256;
  uint32_t MAX_NUM_LOCATIONS_PER_PATH = 256;
  uint32_t MAX_UNIQUE_KMER_POSITIONS = 512;

  /********************************
   * HAPLOTYPE EXTRACTION OPTIONS *
   ********************************/
  std::size_t max_extracted_haplotypes = 42;
  bool skip_breaking_down_extracted_haplotypes = false;

  /*******************
   * LOGGING OPTIONS *
   *******************/
  boost::shared_ptr<boost::log::sinks::synchronous_sink<boost::log::sinks::text_file_backend> > sink;

  /***********
   * METHODS *
   ***********/
  static Options * instance(); // Gets the only instance (singleton pattern)
  void print();


private:
  Options(); // Prevent construction of new instances
  Options(Options const &); // Prevent copy-construction
  Options & operator=(Options const &); // Prevent copy-assignment

  static Options * _instance;
};

} // namespace gyper
