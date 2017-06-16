#pragma once

#include <vector> // std::vector
#include <unordered_map> // std::unordered_map

#include <google/dense_hash_map> // google::dense_hash_map

#include <graphtyper/index/kmer_label.hpp> // gyper::KmerLabel


namespace gyper
{

class Graph;

class MemIndex
{
public:
  uint64_t empty_key = 0ul;
  google::dense_hash_map<uint64_t, std::vector<KmerLabel> > hamming0;
  std::unordered_map<uint64_t, uint64_t> hamming1;

  MemIndex();
  void load();
  void generate_hamming1_hash_map();
  std::vector<std::vector<KmerLabel> > multi_get(std::vector<std::vector<uint64_t> > const & keys);
  std::vector<std::vector<KmerLabel> > multi_get_hamming1(std::vector<std::vector<uint64_t> > const & keys);
};


MemIndex load_secondary_mem_index(std::string const & secondary_index_path, Graph & secondary_graph);

extern MemIndex mem_index;

} // namespace gyper
