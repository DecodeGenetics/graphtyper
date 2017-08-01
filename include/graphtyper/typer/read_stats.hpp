#pragma once

#include <cstdint>
#include <vector>

struct Mismatch
{
  uint32_t pos = 0;
  char new_base = 'N';
};

class ReadStats
{
public:
  ReadStats(std::size_t const num_samples)
  {
    insert_sizes.resize(num_samples);
    mismatches.resize(num_samples);
  }

  std::vector<std::vector<uint32_t> > insert_sizes;
  std::vector<std::vector<Mismatch> > mismatches;
};
