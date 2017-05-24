#pragma once

#include <cstdint>
#include <mutex>
#include <vector>


namespace gyper
{

class GenotypePaths;
class VariantCandidate;

class ReferenceDepth
{
public:
  uint32_t reference_offset = 0;
  std::vector<uint16_t> depth;
  std::vector<uint32_t> local_depth; // List of indexes to increase the depth of, they might be duplicates

  /****************
   * CONSTRUCTORS *
   ****************/
  ReferenceDepth();

  /***************
   * INFORMATION *
   ***************/
  std::size_t start_pos_to_index(uint32_t const start_pos) const;
  std::size_t end_pos_to_index(uint32_t const end_pos) const;
  std::string print_non_zero_depth(std::size_t const MIN_DEPTH = 1) const;
  uint16_t get_read_depth(VariantCandidate const & var) const;

  /****************
   * MODIFICATION *
   ****************/
  void add_genotype_paths(GenotypePaths const & geno);
  void add_reference_depths_from(ReferenceDepth const & ref_depth);
  void clear();
  void increase_local_depth_by_one(std::size_t const start_pos, std::size_t const end_pos);
  void commit_local_depth();
  void resize_depth();
};

class GlobalReferenceDepth
{
public:
  uint32_t reference_offset = 0;
  std::vector<std::vector<uint16_t> > depths;
  std::vector<std::mutex> mutable reference_depth_mutexes;

  /***************
   * INFORMATION *
   ***************/
  std::size_t start_pos_to_index(uint32_t const start_pos) const;
  std::size_t end_pos_to_index(uint32_t const end_pos, std::size_t const depth_size) const;

  /****************
   * MODIFICATION *
   ****************/
  void set_pn_count(std::size_t const pn_count);
  void add_reference_depths_from(ReferenceDepth const & ref_depth, std::size_t const pn_index);
  uint16_t get_read_depth(VariantCandidate const & var, std::size_t const pn_index) const;
  uint64_t get_total_read_depth_of_samples(VariantCandidate const & var, std::vector<uint32_t> const & pn_indexes) const;
};

extern GlobalReferenceDepth global_reference_depth;

} // namespace gyper
