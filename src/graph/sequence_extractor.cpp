#include <iostream>
#include <sstream>

#include <graphtyper/graph/graph.hpp>
#include <graphtyper/graph/graph_serialization.hpp> // gyper::load_graph
#include <graphtyper/graph/sequence_extractor.hpp>
#include <graphtyper/utilities/io.hpp>

namespace gyper
{
void extract_to_fasta(std::string const & file_name, std::string const & graph_path, uint64_t n, uint64_t b, uint64_t e)
{
  load_graph(graph_path);

  std::stringstream ss;
  ss << ">reference\n";
  std::vector<char> ref = graph.get_all_ref();
  ss << std::string(ref.begin(), ref.end()) << std::endl;
  write_to_file(ss.str(), file_name);
  ss.str(std::string());

  for (unsigned i = 0; i < n; ++i)
  {
    std::vector<char> c = graph.walk_random_path(b, e);
    ss << "> " << i << std::endl;
    ss << std::string(c.begin(), c.end()) << std::endl;
    append_to_file(ss.str(), file_name);
    ss.str(std::string());
  }
}

} // namespace gyper
