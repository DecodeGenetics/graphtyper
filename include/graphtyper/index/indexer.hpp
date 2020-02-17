#pragma once

#include <cstdint>
#include <unordered_map>
#include <list>
#include <memory>
#include <fstream>
#include <string>

#include <boost/archive/binary_iarchive.hpp>

#include <graphtyper/constants.hpp>
#include <graphtyper/index/kmer_label.hpp>
#include <graphtyper/index/index_entry.hpp>
#include <graphtyper/index/rocksdb.hpp>


namespace gyper
{

using TEntrySublist = std::deque<IndexEntry>;
using TEntryList = std::deque<TEntrySublist>;

void index_graph(std::string const & index_path);
void index_graph(std::string const & graph_path, std::string const & index_path);
void load_index(std::string const & index_path);
Index<RocksDB> load_secondary_index(std::string const & index_path);

} // namespace gyper
