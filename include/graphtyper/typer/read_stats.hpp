#pragma once

#include <cstdint>
#include <vector>

struct Mismatch
{
  uint32_t pos = 0;
  char new_base = 'N';
}

class ReadStats
{
public:
  std::vector<int32_t> insert_sizes;
  std::vector<Mismatch> mismatches;
};
