#pragma once

#include <deque>

#include <graphtyper/index/index_entry.hpp>
#include <graphtyper/index/ph_index.hpp>


namespace gyper
{

class PHIndex;
class Graph;

using TEntrySublist = std::deque<IndexEntry>;
using TEntryList = std::deque<TEntrySublist>;

PHIndex index_graph(Graph const & graph);
PHIndex index_graph(std::string const & graph_path);

} // namespace gyper
