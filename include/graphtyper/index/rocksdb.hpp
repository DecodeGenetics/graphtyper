#pragma once

#include <string> // std::string
#include <unordered_map>
#include <vector> // std::vector

#include <graphtyper/index/index.hpp>

#include <rocksdb/db.h>
#include <rocksdb/env.h>
#include <rocksdb/slice.h>
#include <rocksdb/options.h>
#include <rocksdb/statistics.h>


namespace gyper
{

std::vector<KmerLabel> value_to_labels(std::string const & value);
uint64_t key_to_uint64_t(std::string const & key_str);


class RocksDB
{
public:
  rocksdb::DB * db = nullptr;
  rocksdb::Status s;
  rocksdb::Options options;
  bool destroy_db_after_use = true;
  bool opened = false;
  std::string filename;

  /** CONSTRUCTORS */
  RocksDB() = default;
  RocksDB & operator=(RocksDB const &) = default;
  RocksDB & operator=(RocksDB &&) = default;
  RocksDB(RocksDB const &) = default;
};


extern Index<RocksDB> index; // Global index

} // namespace gyper
