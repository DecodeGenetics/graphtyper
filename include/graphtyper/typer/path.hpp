#pragma once

#include <bitset> // std::bitset<Size>
#include <cstdint> // uint16_t, uint32_t
#include <vector> // std::vector<Type>

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/graph.hpp>
#include <graphtyper/graph/haplotype.hpp>
#include <graphtyper/index/kmer_label.hpp>

namespace gyper
{

class Path
{
public:
  uint32_t start = 0; /** \brief The start position of this path. Can be a special position and is used for connecting two paths. */
  uint32_t end = 0; /** \brief The end position of this path. Can be a special position and is used for connecting two paths. */
  uint16_t read_start_index = 0;
  uint16_t read_end_index = 0;
  std::vector<uint32_t> var_order;
  std::vector<std::bitset<MAX_NUMBER_OF_HAPLOTYPES> > nums;
  std::size_t mismatches;

  /*********************
   * PATH CONSTRUCTORS *
   *********************/
  Path() noexcept;
  Path(KmerLabel const & l, uint16_t const read_start_index, uint16_t const read_end_index, uint16_t const mismatches = 0) noexcept;
  Path(Path const & p1, Path const & p2) noexcept;

  /**********************
   * PATH MODIFICATIONS *
   **********************/
  void merge_with_current(KmerLabel const & l);

  /********************
   * PATH INFORMATION *
   ********************/
  uint32_t start_pos() const;
  uint32_t end_pos() const;
  uint32_t start_correct_pos() const;
  uint32_t end_correct_pos() const;
  uint32_t start_ref_reach_pos() const;
  uint32_t end_ref_reach_pos() const;
  uint32_t size() const;
  uint32_t get_read_end_index(uint32_t const read_length) const;
};

} // namespace gyper
