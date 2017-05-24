#pragma once

#include <cstdint> // intX_t, uintX_t
#include <map> // std::map
#include <string> // std::string
#include <vector> // std::vector


namespace gyper
{

class VarStats
{
public:
  /** Mapping quality (MQ) */
  uint64_t mapq_root_total = 0u;
  uint32_t mapq_count = 0u;
  uint32_t mapq_zero_count = 0u;

  /** Clipped reads */
  uint32_t clipped_reads = 0u;

  /** Realignment statistics */
  uint32_t unaligned_reads = 0u;
  std::vector<uint32_t> realignment_distance;
  std::vector<uint32_t> realignment_count;

  /** Strand bias per allele */
  std::vector<uint32_t> r1_strand_forward;
  std::vector<uint32_t> r1_strand_reverse;
  std::vector<uint32_t> r2_strand_forward;
  std::vector<uint32_t> r2_strand_reverse;

  /** Graph complexity */
  uint8_t graph_complexity = 0u;

  /**
   * CONSTRUCTORS
   */
  VarStats(uint16_t allele_count) noexcept;

  /**
   * MODIFIERS
   */
  void add_mapq(uint8_t const new_mapq);
  void add_realignment_distance(uint8_t const allele_id, uint32_t const original_pos, uint32_t const new_pos);

  /**
   * CLASS INFORMATION
   */
  uint8_t get_rms_mapq() const;
  std::string get_realignment_count() const;
  std::string get_realignment_distance() const;
  std::string get_forward_strand_bias() const;
  std::string get_reverse_strand_bias() const;
  std::string get_unaligned_count() const;

  // Read pair specfici biases
  std::string get_r1_forward_strand_bias() const;
  std::string get_r2_forward_strand_bias() const;
  std::string get_r1_reverse_strand_bias() const;
  std::string get_r2_reverse_strand_bias() const;
};

std::string join_strand_bias(std::vector<uint32_t> const & bias);
std::string join_strand_bias(std::vector<uint32_t> const & r1bias, std::vector<uint32_t> const & r2bias);
std::vector<std::string> split_bias_to_strings(std::string const & bias);
std::vector<uint32_t> split_bias_to_numbers(std::string const & bias);
std::vector<uint32_t> get_strand_bias(std::map<std::string, std::string> const & infos, std::string const & bias);
uint32_t get_accumulated_strand_bias(std::map<std::string, std::string> const & infos, std::string const & bias);
std::vector<uint16_t> get_list_of_uncalled_alleles(std::string const & ac);

} // namespace gyper
