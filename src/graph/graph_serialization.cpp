#include <fstream>
#include <string>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/log/trivial.hpp>

#include <graphtyper/graph/graph.hpp>
#include <graphtyper/graph/graph_serialization.hpp>


namespace gyper
{

void
save_graph(std::string const & graph_path)
{
  std::ofstream ofs(graph_path.c_str(), std::ios::binary);

  if (!ofs.is_open())
  {
    BOOST_LOG_TRIVIAL(fatal) << "[graphtyper::constructor] Could not save graph at '" << graph_path << "'";
    std::exit(1);
  }

  boost::archive::binary_oarchive oa(ofs);
  oa << graph;
}


void
load_graph(std::string const & graph_path)
{
  graph.clear();
  graph = Graph();
  std::ifstream ifs(graph_path.c_str(), std::ios::binary);

  if (!ifs.is_open())
  {
    BOOST_LOG_TRIVIAL(fatal) << "[graphtyper::constructor] Could not load graph at '" << graph_path << "'";
    std::exit(1);
  }

  boost::archive::binary_iarchive ia(ifs);
  ia >> graph;
  assert(graph.size() > 0u);

  // Create a reference genome each time the graph is loaded
  graph.generate_reference_genome();
}


Graph
load_secondary_graph(std::string const & graph_path)
{
  Graph second_graph = Graph();
  std::ifstream ifs(graph_path.c_str(), std::ios::binary);

  if (!ifs.is_open())
  {
    BOOST_LOG_TRIVIAL(fatal) << "[graphtyper::constructor] Could not load graph at '" << graph_path << "'";
    std::exit(1);
  }

  boost::archive::binary_iarchive ia(ifs);
  ia >> second_graph;
  assert(second_graph.size() > 0u);

  // Create a reference genome each time the graph is loaded
  second_graph.generate_reference_genome();

  return second_graph;
}


} // namespace gyper
