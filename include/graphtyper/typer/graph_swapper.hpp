#pragma once

#include <graphtyper/graph/graph.hpp>
#include <graphtyper/index/mem_index.hpp>


namespace gyper
{

void swap_graph_and_index(Graph && secondary_graph, MemIndex && secondary_mem_index);

} // namespace gyper
