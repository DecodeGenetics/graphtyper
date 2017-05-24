#pragma once

#include <vector>

namespace gyper
{

class AlleleCover
{
public:
  uint16_t marker_size;
  uint16_t fully_covered_count = 0; // The number reads which fully cover this variant
  std::vector<int16_t> partial_covers; // Each entry represents how a single read

  AlleleCover(uint16_t const _marker_size = 0)
    : marker_size(_marker_size)
  {}
};

} // namespace gyper
