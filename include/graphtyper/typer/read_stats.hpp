#pragma once

#include <cstdint>
#include <vector>

class ReadStats
{
public:
  ReadStats(std::size_t const num_samples) :
    insert_sizes(num_samples),
    mismatches(num_samples),
    ref_read_abs_pos(num_samples),
    alt_read_abs_pos(num_samples),
    other_read_abs_pos(num_samples)
  {
  }

  std::vector<std::vector<uint32_t>> insert_sizes;
  std::vector<std::vector<uint8_t>> mismatches;
  std::vector<std::vector<std::pair<uint32_t, uint32_t>>> ref_read_abs_pos;
  std::vector<std::vector<std::pair<uint32_t, uint32_t>>> alt_read_abs_pos;
  std::vector<std::vector<std::pair<uint32_t, uint32_t>>> other_read_abs_pos;
};
