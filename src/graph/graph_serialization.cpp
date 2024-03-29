#include <fstream>
#include <string>

#include <cereal/archives/binary.hpp>

#include <graphtyper/graph/absolute_position.hpp>
#include <graphtyper/graph/graph.hpp>
#include <graphtyper/graph/graph_serialization.hpp>
#include <graphtyper/utilities/logging.hpp>

namespace gyper
{
void save_graph(std::string const & graph_path)
{
  std::ofstream ofs(graph_path.c_str(), std::ios::binary);

  if (!ofs.is_open())
  {
    print_log(log_severity::error, "[graphtyper::constructor] Could not save graph at '", graph_path, "'");
    std::exit(1);
  }

  cereal::BinaryOutputArchive oa(ofs);
  oa << graph;
}

void load_graph(std::string const & graph_path)
{
  gyper::graph.clear();
  gyper::graph = Graph();
  std::ifstream ifs(graph_path.c_str(), std::ios::binary);

  if (!ifs.is_open())
  {
    print_log(log_severity::error, "[graphtyper::constructor] Could not load graph at '", graph_path, "'");
    std::exit(1);
  }

  cereal::BinaryInputArchive ia(ifs);
  ia >> graph;
  assert(graph.size() > 0u);

  // Create a reference genome each time the graph is loaded
  graph.generate_reference_genome();
  absolute_pos.calculate_offsets(graph.contigs);
}

Graph load_secondary_graph(std::string const & graph_path)
{
  Graph second_graph = Graph();
  std::ifstream ifs(graph_path.c_str(), std::ios::binary);

  if (!ifs.is_open())
  {
    print_log(log_severity::error, "[graphtyper::constructor] Could not load graph at '", graph_path, "'");
    std::exit(1);
  }

  cereal::BinaryInputArchive ia(ifs);
  ia >> second_graph;
  assert(second_graph.size() > 0u);

  // Create a reference genome each time the graph is loaded
  second_graph.generate_reference_genome();

  return second_graph;
}

} // namespace gyper
