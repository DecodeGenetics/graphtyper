#pragma once

#include <string>

#include <graphtyper/graph/graph.hpp>

namespace gyper
{
void save_graph(std::string const & graph_path);
void load_graph(std::string const & graph_path);
Graph load_secondary_graph(std::string const & graph_path);

} // namespace gyper
