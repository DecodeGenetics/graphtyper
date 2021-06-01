#pragma once

#include <cstdint> // intX_t, uintX_t
#include <map>     // std::map
#include <string>  // std::string
#include <vector>  // std::vector

#include <cereal/access.hpp>
#include <cereal/types/utility.hpp>

#include <graphtyper/graph/read_strand.hpp>

namespace gyper
{
struct VarStatsPerAllele
{
  uint64_t clipped_bp{0u};
  uint64_t mapq_squared{0u};
  uint32_t score_diff{0u};
  uint32_t mismatches{0u};
  uint64_t qd_qual{0u};
  uint64_t qd_depth{0u};
  uint64_t total_depth{0u};
  uint32_t ac{0};
  uint32_t pass_ac{0};
  uint16_t maximum_alt_support{0u};
  double maximum_alt_support_ratio{0u};
  std::pair<uint32_t, uint32_t> het_multi_allele_depth{0u, 0u};
  std::pair<uint32_t, uint32_t> hom_multi_allele_depth{0u, 0u};

  template <class Archive>
  void serialize(Archive & ar, unsigned int /*version*/)
  {
    ar & clipped_bp;
    ar & mapq_squared;
    ar & score_diff;
    ar & mismatches;
    ar & qd_qual;
    ar & qd_depth;
    ar & total_depth;
    ar & ac;
    ar & pass_ac;
    ar & maximum_alt_support;
    ar & maximum_alt_support_ratio;
    ar & het_multi_allele_depth;
    ar & hom_multi_allele_depth;
  }
};

class VarStats
{
  friend class cereal::access;

public:
  /** Stats per allele (excluding ReadStrand) */
  std::vector<VarStatsPerAllele> per_allele{};

  /** Strand bias per allele */
  std::vector<ReadStrand> read_strand{};

  /** Clipped reads */
  uint32_t clipped_reads{0u};

  /** MapQ statistics */
  uint64_t mapq_squared{0u};

  // Other statistics for INFO field
  uint32_t n_genotyped{0};
  uint32_t n_calls{0};
  uint32_t n_passed_calls{0};
  uint32_t n_ref_ref{0};
  uint32_t n_ref_alt{0};
  uint32_t n_alt_alt{0};
  uint8_t n_max_alt_proper_pairs{0};
  uint64_t seqdepth{0};

  std::pair<uint32_t, uint32_t> het_allele_depth = {0ul, 0ul}; // First is the first call, second is the second call
  std::pair<uint32_t, uint32_t> hom_allele_depth = {0ul,
                                                    0ul}; // First is the called allele, second is not the called one

  /**
   * CONSTRUCTORS
   */
  VarStats() = default;
  VarStats(std::size_t allele_count) noexcept;

  /**
   * CLASS INFORMATION
   */

  void write_stats(std::map<std::string, std::string> & infos) const;
  void write_per_allele_stats(std::map<std::string, std::string> & infos) const;
  void write_read_strand_stats(std::map<std::string, std::string> & infos) const;

  /**
   * CLASS MODIFIERS
   */
  void add_stats(VarStats const & stats);
  void read_stats(std::map<std::string, std::string> const & infos);

private:
  template <class Archive>
  void serialize(Archive & ar, unsigned int version);
};

template <class T>
std::string join_strand_bias(std::vector<T> const & bias);
std::string join_strand_bias(std::vector<uint32_t> const & r1bias, std::vector<uint32_t> const & r2bias);

std::vector<std::string> split_bias_to_strings(std::string const & bias);
std::vector<uint32_t> split_bias_to_numbers(std::string const & bias);
std::vector<uint32_t> get_strand_bias(std::map<std::string, std::string> const & infos, std::string const & bias);
long get_accumulated_strand_bias(std::map<std::string, std::string> const & infos, std::string const & bias);
long get_accumulated_alt_strand_bias(std::map<std::string, std::string> const & infos, std::string const & bias);
std::vector<uint16_t> get_list_of_uncalled_alleles(std::string const & ac);

} // namespace gyper
