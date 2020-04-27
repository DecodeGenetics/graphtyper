#pragma once

#include <cstdint>
#include <unordered_map>
#include <vector>

#include <parallel_hashmap/phmap.h>

#include <graphtyper/index/kmer_label.hpp>


namespace gyper
{

class PHIndex
{
public:
  using PHtype = phmap::flat_hash_map<uint64_t, std::vector<KmerLabel> >;
  //using PHtype = phmap::node_hash_map<uint64_t, std::vector<KmerLabel> >;
  //using PHtype = std::unordered_map<uint64_t, std::vector<KmerLabel> >;

  PHtype hamming0;

  PHIndex() = default;
  PHIndex(PHIndex const &) = delete; // No copy
  PHIndex(PHIndex &&) = default;
  PHIndex & operator=(PHIndex const &) = delete; // No copy
  PHIndex & operator=(PHIndex &&) = default;
  ~PHIndex() = default;

  void put(uint64_t const key, KmerLabel && label);
  void put(uint64_t const key, std::vector<KmerLabel> && labels);

  bool check() const;
  std::vector<KmerLabel> get(uint64_t const key) const;
  std::vector<KmerLabel> get(std::vector<uint64_t> const & keys) const;
  std::vector<std::vector<KmerLabel> > multi_get(std::vector<std::vector<uint64_t> > const & keys) const;
};

} // namespace gyper
