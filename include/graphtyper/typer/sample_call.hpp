#pragma once

#include <cstdint> // uint8_t, uint16_t
#include <vector> // std::vector

#include <graphtyper/utilities/graph_help_functions.hpp>

namespace gyper
{

class SampleCall
{
public:
  SampleCall() noexcept;
  SampleCall(std::vector<uint8_t> const & phred, std::vector<uint16_t> const & coverage/*, uint8_t const & ambiguous_depth*/) noexcept;

  // If R is the number of alleles, then phred should be of size R * (R + 1) / 2 and coverage of size R.
  std::vector<uint8_t> phred; // GT and PL
  std::vector<uint16_t> coverage; // AD and DP
  // uint8_t ambiguous_depth = 0;

  std::pair<uint16_t, uint16_t> get_gt_call() const;
  uint8_t get_gq() const;

};

} // namespace gyper
