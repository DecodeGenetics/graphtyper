#pragma once

#include <cstdint>
#include <vector>

#include <parallel_hashmap/phmap.h>

#include <graphtyper/index/kmer_label.hpp>


namespace gyper
{

class PHIndex
{
public:
  using PHtype = phmap::parallel_flat_hash_map<uint64_t, std::vector<KmerLabel> >;

  PHtype hamming0;

  PHIndex() = default;

  void put(uint64_t const key, KmerLabel && label);
  void put(uint64_t const key, std::vector<KmerLabel> && labels);

  bool check() const;
  std::vector<KmerLabel> get(uint64_t const key) const;
  std::vector<KmerLabel> get(std::vector<uint64_t> const & keys) const;
  std::vector<std::vector<KmerLabel> > multi_get(std::vector<std::vector<uint64_t> > const & keys) const;
};

} // namespace gyper
