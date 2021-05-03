#pragma once

#include <bitset> // std::bitset<Size>
#include <cstdint> // uint16_t, uint32_t
#include <set>
#include <vector> // std::vector<Type>

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/haplotype.hpp>
#include <graphtyper/index/kmer_label.hpp>

namespace gyper
{

class Graph;

class Path
{
public:
  /** \brief The start position of this path. Can be a special position and is used for connecting
   * two paths. */
  uint32_t start = 0;

  /** \brief The end position of this path. Can be a special position and is used for connecting
   * two paths. */
  uint32_t end = 0;

  /** \brief 0-based start index of the part of the read which aligned. */
  uint16_t read_start_index = 0;

  /** \brief 0-based end index of the part of the read which aligned.
   *
   * The reported value is included.
   */
  uint16_t read_end_index = 0;
  std::vector<uint32_t> var_order;
  std::vector<std::set<uint64_t>> nums = std::vector<std::set<uint64_t>>(0);
  uint16_t mismatches;

  /*********************
   * PATH CONSTRUCTORS *
   *********************/
  Path() = default;
  Path(Graph const & graph,
       KmerLabel const & l,
       uint16_t read_start_index,
       uint16_t read_end_index,
       uint16_t mismatches = 0) noexcept;

  Path(Path const & p1, Path const & p2) noexcept;
  Path(Path const &)                     noexcept = default;
  Path(Path &&)                          noexcept = default;
  Path & operator=(Path const &)         noexcept = default;
  Path & operator=(Path &&)              noexcept = default;

  /**********************
   * PATH MODIFICATIONS *
   **********************/
  void erase_var_order(long index);
  void erase_ref_support(long index);
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
  uint32_t get_read_end_index(uint32_t read_length) const;
  bool is_reference() const;
  bool is_purely_reference() const;
  bool is_empty() const;
};

} // namespace gyper
