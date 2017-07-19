#pragma once

#include <cstdint> // uint8_t, uint16_t
#include <vector>

#include <seqan/sequence.h>

#include <graphtyper/graph/graph.hpp>
#include <graphtyper/index/mem_index.hpp>
#include <graphtyper/typer/path.hpp>
#include <graphtyper/typer/variant_candidate.hpp>


namespace gyper
{

class GenotypePaths
{
public:
  std::vector<char> read;
  std::vector<char> qual;
  std::vector<Path> paths;
  uint32_t longest_path_length = 0;
  uint8_t mapq = 255u;
  uint32_t original_pos = 0; /* 0-based position from global alignment */
  bool is_first_in_pair = true;
  bool forward_strand = true;
  bool is_originally_unaligned = false;
  bool is_originally_clipped = false;

  /****************
   * CONSTRUCTORS *
   ****************/
  GenotypePaths();
  GenotypePaths(seqan::IupacString const & read,
                seqan::CharString const & qual,
                uint8_t const mapq = 255u
                );

  /***********************
   * CLASS MODIFICATIONS *
   ***********************/
  void add_next_kmer_labels(std::vector<KmerLabel> const & ll, uint32_t const start_index, uint32_t const read_end_index, std::size_t const mismatches = 0);
  void add_prev_kmer_labels(std::vector<KmerLabel> const & ll, uint32_t const read_start_index, uint32_t const read_end_index, std::size_t const mismatches = 0);

  void clear_paths();

  void walk_read_ends(seqan::IupacString const & read, int maximum_mismatches = -1, gyper::Graph const & graph = gyper::graph);
  void walk_read_starts(seqan::IupacString const & read, int maximum_mismatches = -1, gyper::Graph const & graph = gyper::graph);

  // Path filtering
  void remove_short_paths();
  void remove_paths_with_no_variants();
  void remove_paths_with_too_many_mismatches();

  std::vector<VariantCandidate> find_new_variants() const;

  void update_longest_path_size();

  /*********************
   * CLASS INFORMATION *
   *********************/
  std::size_t longest_path_size() const;
  std::vector<Path> longest_paths() const;
  bool all_paths_unique() const;
  bool all_paths_fully_aligned() const;
  bool is_purely_reference() const;
  bool check_no_variant_is_missing() const;

};

int compare_pair_of_genotype_paths(GenotypePaths const & geno1, GenotypePaths const & geno2);
int compare_pair_of_genotype_paths(std::pair<GenotypePaths, GenotypePaths> const & genos1, std::pair<GenotypePaths, GenotypePaths> const & genos2);

} // namespace gyper
