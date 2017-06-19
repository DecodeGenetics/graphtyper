#include <array> // std::array
#include <cstdint> // uint64_t
#include <vector> // std::vector
#include <unordered_map> // std::unordered_map

#include <graphtyper/graph/graph.hpp> // gyper::Graph
#include <graphtyper/index/indexer.hpp>
#include <graphtyper/index/rocksdb.hpp> // gyper::index
#include <graphtyper/index/mem_index.hpp> // gyper::MemIndex
#include <graphtyper/utilities/type_conversions.hpp>

#include <google/dense_hash_map> // google::dense_hash_map


namespace gyper
{

MemIndex::MemIndex()
{}


void
MemIndex::load()
{
  assert(index.hamming0.db); // Index is open
  assert(index.opened);
  this->hamming0 = google::dense_hash_map<uint64_t, std::vector<KmerLabel> >();

  for (uint64_t key = 0; key < 0xFFFFFFFFFFFFFFFFull; ++key)
  {
    if (!index.exists(key))
    {
      empty_key = key; // Empty key is a class member
      break;
    }
  }

  this->hamming0.set_empty_key(empty_key);
  std::size_t num_keys = 0;
  rocksdb::Iterator* it = index.hamming0.db->NewIterator(rocksdb::ReadOptions());
  assert (it);

  for (it->SeekToFirst(); it->Valid(); it->Next())
  {
    uint64_t const key = key_to_uint64_t(it->key().ToString());
    ++num_keys;
    this->hamming0[key] = value_to_labels(it->value().ToString());
  }

  assert(it->status().ok()); // Check for any errors
  delete it;
}


void
MemIndex::generate_hamming1_hash_map()
{
  hamming1.clear();
  std::vector<uint64_t> duplicate_keys;

  for (auto it = hamming0.begin(); it != hamming0.end(); ++it)
  {
    std::array<uint64_t, 96> hamming1_keys = to_uint64_vec_hamming_distance_1(it->first);

    for (auto const & hamming1_key : hamming1_keys)
    {
      auto find_it = hamming1.find(hamming1_key);

      if (find_it == hamming1.end())
        hamming1[hamming1_key] = it->first;
      else
        duplicate_keys.push_back(hamming1_key);
    }
  }

  // Delete all duplicated keys
  for (auto const & d_key : duplicate_keys)
    hamming1.erase(d_key);
}


std::vector<std::vector<KmerLabel> >
MemIndex::multi_get(std::vector<std::vector<uint64_t> > const & keys) const
{
  std::vector<std::vector<KmerLabel> > labels(keys.size());

  for (std::size_t i = 0; i < keys.size(); ++i)
  {
    for (std::size_t j = 0; j < keys[i].size(); ++j)
    {
      if (keys[i][j] != empty_key)
      {
        auto find_it = hamming0.find(keys[i][j]);

        if (find_it != hamming0.end())
          std::copy(find_it->second.begin(), find_it->second.end(), std::back_inserter(labels[i]));
      }
    }
  }

  return labels;
}


std::vector<std::vector<KmerLabel> >
MemIndex::multi_get_hamming1(std::vector<std::vector<uint64_t> > const & keys) const
{
  std::vector<std::vector<KmerLabel> > labels(keys.size());

  for (std::size_t i = 0; i < keys.size(); ++i)
  {
    for (std::size_t j = 0; j < keys[i].size(); ++j)
    {
      auto find_it = hamming1.find(keys[i][j]);

      if (find_it != hamming1.end())
      {
        auto hamming0_find_it = hamming0.find(find_it->second);
        assert(hamming0_find_it != hamming0.end());
        std::copy(hamming0_find_it->second.begin(), hamming0_find_it->second.end(), std::back_inserter(labels[i]));
      }
    }
  }

  assert(keys.size() == labels.size());
  return labels;
}


MemIndex
load_secondary_mem_index(std::string const & secondary_index_path, Graph & secondary_graph)
{
  // Swap graphs
  std::swap(graph, secondary_graph);

  Index<RocksDB> secondary_index = load_secondary_index(secondary_index_path);
  std::swap(index, secondary_index);
  MemIndex secondary_mem_index;
  secondary_mem_index.load();
  std::swap(index, secondary_index);
  secondary_index.close();

  // Swap graphs back to the way they were
  std::swap(graph, secondary_graph);

  return secondary_mem_index;
}


MemIndex mem_index;

}
