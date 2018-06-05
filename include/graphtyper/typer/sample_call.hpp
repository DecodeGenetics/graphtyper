#pragma once

#include <cstdint> // uint8_t, uint16_t
#include <vector> // std::vector

namespace gyper
{

struct AlleleCoverage;

class SampleCall
{
public:
  SampleCall() noexcept;
  SampleCall(std::vector<uint8_t> const & phred,
             std::vector<uint16_t> const & coverage,
             uint8_t ambiguous_depth,
             uint8_t ambiguous_depth_alt,
             uint8_t alt_proper_pair_depth
    ) noexcept;

  // If R is the number of alleles, then phred should be of size R * (R + 1) / 2 and coverage of size R.
  std::vector<uint8_t> phred; // GT and PL
  std::vector<uint16_t> coverage; // AD and DP
  uint16_t ref_total_depth = 0u; // Total reads that support the reference allele
  uint16_t alt_total_depth = 0u; // Total reads that support alternative allele(s)
  uint8_t ambiguous_depth = 0u; //
  uint8_t alt_proper_pair_depth = 0u; // Total reads in proper pairs that support the alternative allele(s)

  std::pair<uint16_t, uint16_t> get_gt_call() const;
  uint8_t get_gq() const;

  /**
   * MODIFIERS
   */
  void change_to_ref_vs_all();

};

/**
 * \brief Create a biallelic call from a multi-allelic call
 * \param old_call Old multi-allelic call
 * \param aa index of alternative allele
 * \return Biallelic call for alternative allele at index aa
 */
SampleCall
make_bi_allelic_call(SampleCall const & old_call, long aa);

} // namespace gyper
