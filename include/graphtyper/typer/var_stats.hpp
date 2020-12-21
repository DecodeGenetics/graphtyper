#pragma once

#include <cstdint> // intX_t, uintX_t
#include <map> // std::map
#include <string> // std::string
#include <vector> // std::vector

#include <boost/serialization/access.hpp>


#include <graphtyper/graph/read_strand.hpp>

namespace gyper
{

struct VarStatsPerAllele
{
  uint64_t clipped_bp{0u};
  uint64_t mapq_squared{0u};
  uint32_t score_diff{0u};
  uint32_t mismatches{0u};

  template <class Archive>
  void
  serialize(Archive & ar, unsigned int /*version*/)
  {
    ar & clipped_bp;
    ar & mapq_squared;
    ar & score_diff;
    ar & mismatches;
  }


};


class VarStats
{
  friend class boost::serialization::access;

public:
  /** Stats per allele (excluding ReadStrand) */
  std::vector<VarStatsPerAllele> per_allele{};

  /** Strand bias per allele */
  std::vector<ReadStrand> read_strand{};

  /** Clipped reads */
  uint32_t clipped_reads{0u};

  /** MapQ statistics */
  uint64_t mapq_squared{0u};

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
