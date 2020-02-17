#pragma once

#include <cstdint> // intX_t, uintX_t
#include <map> // std::map
#include <string> // std::string
#include <vector> // std::vector

#include <graphtyper/graph/read_strand.hpp>

namespace gyper
{

class VarStats
{
public:
  /** Clipped reads */
  uint32_t clipped_reads = 0u;

  /** Strand bias per allele */
  std::vector<ReadStrand> read_strand{};

  /** MapQ statistics */
  uint64_t mapq_squared = 0u;

  /**
   * CONSTRUCTORS
   */
  VarStats(uint16_t allele_count) noexcept;

  /**
   * MODIFIERS
   */
  void add_mapq(uint8_t const new_mapq);

  /**
   * CLASS INFORMATION
   */
  std::string get_forward_strand_bias() const;
  std::string get_reverse_strand_bias() const;

  // Read pair specfici biases
  std::string get_r1_forward_strand_bias() const;
  std::string get_r2_forward_strand_bias() const;
  std::string get_r1_reverse_strand_bias() const;
  std::string get_r2_reverse_strand_bias() const;


};

template <class T>
std::string join_strand_bias(std::vector<T> const & bias);
std::string join_strand_bias(std::vector<uint32_t> const & r1bias, std::vector<uint32_t> const & r2bias);

std::vector<std::string> split_bias_to_strings(std::string const & bias);
std::vector<uint32_t> split_bias_to_numbers(std::string const & bias);
std::vector<uint32_t> get_strand_bias(std::map<std::string, std::string> const & infos, std::string const & bias);
uint32_t get_accumulated_strand_bias(std::map<std::string, std::string> const & infos, std::string const & bias);
std::vector<uint16_t> get_list_of_uncalled_alleles(std::string const & ac);

} // namespace gyper
