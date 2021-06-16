#pragma once

#include <cstdint> // uint64_t
#include <string>  // std::string
#include <vector>  // std::vector

#include <graphtyper/index/kmer_label.hpp> // gyper::KmerLabel

namespace gyper
{
class Graph;

template <typename HashTable>
class Index
{
public:
  std::size_t const MAX_BUFFER = 100000000; // One-hundred million or in worst case 1800 MB if we disregard overhead
                                            // (which is probably something like 2-5x)
  std::unordered_map<uint64_t, std::vector<KmerLabel>> buffer_map;
  bool opened = false;
  HashTable hamming0;

  Index();
  Index(Index const & cp_index); // = default;
  Index(Index && mv_index);      // = default;
  Index(std::string const & f, bool const clear_first = false, bool const read_only = false);
  ~Index();

  Index & operator=(Index &&);

  bool exists(uint64_t const key) const;
  std::vector<KmerLabel> get(uint64_t const key) const;
  std::vector<std::vector<KmerLabel>> multi_get(std::vector<std::vector<uint64_t>> const & multi_keys) const;

  void open(std::string const & f, bool clear_first = false, bool read_only = false);
  void close();
  void put(uint64_t const key, KmerLabel && value);
  void put(uint64_t const key, std::vector<KmerLabel> && values);
  bool check();
  void commit();
  std::size_t size();
  void clear();
  void construct(bool const read_only = false);
};

} // namespace gyper
