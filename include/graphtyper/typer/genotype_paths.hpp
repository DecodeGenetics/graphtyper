#pragma once

#include <cstdint> // uint8_t, uint16_t
#include <memory>
#include <string> // std::string
#include <vector>

#include <seqan/sequence.h>

#include <graphtyper/constants.hpp> // INSERT_SIZE_WHEN_NOT_PROPER_PAIR

namespace gyper
{

// Forward declarations
class Graph;
class KmerLabel;
class Path;
class VariantCandidate;

struct GenotypePathsDetails
{
  std::string query_name;
  std::string read_group;
  uint32_t score_diff = 0u;
};


class GenotypePaths
{
public:
  std::vector<char> read2;
  std::vector<char> qual2;
  std::vector<Path> paths;
  uint16_t read_length = 0;
  uint16_t flags = 0;
  uint32_t longest_path_length = 0;
  uint32_t original_pos = 0; // 0-based position from global alignment
  uint8_t mapq = 255;
  int32_t ml_insert_size = INSERT_SIZE_WHEN_NOT_PROPER_PAIR;

#ifndef NDEBUG
  std::unique_ptr<GenotypePathsDetails> details; // Only used when statistics are kept
#endif // NDEBUG

  /****************
   * CONSTRUCTORS *
   ****************/
  GenotypePaths();
  explicit GenotypePaths(int16_t _flags, std::size_t _read_length);
  explicit GenotypePaths(GenotypePaths const & b);
  explicit GenotypePaths(GenotypePaths && b);
  GenotypePaths & operator=(GenotypePaths const &);
  GenotypePaths & operator=(GenotypePaths &&);
  ~GenotypePaths() = default;

  /***********************
   * CLASS MODIFICATIONS *
   ***********************/
  void add_next_kmer_labels(std::vector<KmerLabel> const & ll,
                            uint32_t start_index,
                            uint32_t read_end_index,
                            int mismatches = 0
                            );

  void add_prev_kmer_labels(std::vector<KmerLabel> const & ll,
                            uint32_t const read_start_index,
                            uint32_t const read_end_index,
                            int const mismatches = 0
                            );

  void clear_paths();

  void walk_read_ends(seqan::IupacString const & read,
                      int maximum_mismatches, // -1 if no limit
                      gyper::Graph const & graph);

  void walk_read_starts(seqan::IupacString const & read,
                        int maximum_mismatches, // -1 if no limit
                        gyper::Graph const & graph);

  // Path filtering
  void remove_short_paths();
  void remove_support_from_read_ends();
  void remove_paths_within_variant_node();
  void remove_paths_with_too_many_mismatches();
  void remove_non_ref_paths_when_read_matches_ref();
  void remove_fully_special_paths();

  std::vector<VariantCandidate> find_new_variants() const;

  void update_longest_path_size();

  /*********************
   * CLASS INFORMATION *
   *********************/
  std::size_t longest_path_size() const;
//  std::vector<Path> longest_paths() const;
  bool all_paths_unique() const;
  bool all_paths_fully_aligned() const;
  bool is_purely_reference() const;
  bool check_no_variant_is_missing() const;

#ifndef NDEBUG
  std::string to_string() const;
#endif // NDEBUG

  bool is_proper_pair() const;

};

int compare_pair_of_genotype_paths(GenotypePaths const & geno1, GenotypePaths const & geno2);
int compare_pair_of_genotype_paths(std::pair<GenotypePaths *, GenotypePaths *> const & genos1_ptr,
                                   std::pair<GenotypePaths *, GenotypePaths *> const & genos2_ptr);

int compare_pair_of_genotype_paths(std::pair<GenotypePaths, GenotypePaths> & genos1,
                                   std::pair<GenotypePaths, GenotypePaths> & genos2);

} // namespace gyper
