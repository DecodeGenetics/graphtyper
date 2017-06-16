#include <utility> // std::swap

#include <graphtyper/graph/graph.hpp> // gyper::Graph
#include <graphtyper/typer/graph_swapper.hpp>

namespace gyper
{

void
swap_graph_and_index(Graph & secondary_graph, MemIndex & secondary_mem_index)
{
  std::swap(graph, secondary_graph);
  std::swap(mem_index, secondary_mem_index);
}

} // namespace gyper
