#include <catch.hpp>

#include <stdio.h>
#include <climits>
#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/graph.hpp>
#include <graphtyper/graph/constructor.hpp>
#include <graphtyper/graph/graph_serialization.hpp>
#include <graphtyper/graph/label.hpp>
#include <graphtyper/utilities/type_conversions.hpp>


namespace
{

bool inline
is_file(std::string filename)
{
  struct stat sb;

  if (stat(filename.c_str(), &sb) == 0 && (S_ISREG(sb.st_mode) || S_ISLNK(sb.st_mode)))
  {
    return true;
  }

  return false;
}


} // anon namespace


void
create_graph(std::string fasta, std::string vcf, std::vector<std::string> regions, bool use_absolute_positions = true)
{
  std::stringstream vcf_path;
  vcf_path << gyper_SOURCE_DIRECTORY << vcf;
  std::stringstream reference_path;
  reference_path << gyper_SOURCE_DIRECTORY << fasta;

  gyper::construct_graph(reference_path.str().c_str(), vcf_path.str().c_str(), regions, use_absolute_positions);
  std::stringstream graph_path;
  graph_path << gyper_SOURCE_DIRECTORY << "/test/data/graphs";

  std::string graph_directory(graph_path.str());

  // Check if directory exists, and of not, create it
  struct stat st;

  if (stat(graph_directory.c_str(), &st) == -1)
  {
    mkdir(graph_directory.c_str(), 0755);
  }

  graph_path << fasta.substr(20, fasta.size() - 3 - 20);

  for (auto region : regions)
  {
    graph_path << "_" << region;
  }

  graph_path << ".grf";
  REQUIRE(gyper::graph.size() != 0);

  {
    std::ofstream ofs(graph_path.str().c_str(), std::ios::binary);
    REQUIRE(ofs.is_open());
    boost::archive::binary_oarchive oa(ofs);
    oa << gyper::graph;
  }

  // test open
  {
    gyper::Graph new_graph;
    std::ifstream ifs(graph_path.str().c_str(), std::ios::binary);
    REQUIRE(ifs.is_open());
    boost::archive::binary_iarchive ia(ifs);
    ia >> new_graph;

    REQUIRE(new_graph.size() == gyper::graph.size());
    REQUIRE(new_graph.get_genomic_region().rID == gyper::graph.get_genomic_region().rID);
    REQUIRE(new_graph.get_genomic_region().chr == gyper::graph.get_genomic_region().chr);
    REQUIRE(new_graph.get_genomic_region().begin == gyper::graph.get_genomic_region().begin);
    REQUIRE(new_graph.get_genomic_region().end == gyper::graph.get_genomic_region().end);
    // REQUIRE(new_graph.actual_poses == gyper::graph.actual_poses);
  }
}


void
create_graph(std::string fasta, std::string vcf, std::string region, bool use_absolute_positions = true)
{
  std::vector<std::string> regions = {region};
  create_graph(fasta, vcf, regions, use_absolute_positions);
}


TEST_CASE("Construct test graph (chr1)")
{
  create_graph("/test/data/reference/index_test.fa", "/test/data/reference/index_test.vcf.gz", "chr1", false);

  REQUIRE(gyper::graph.ref_nodes.size() == 2);
  REQUIRE(gyper::graph.var_nodes.size() == 2);

  std::vector<gyper::RefNode> const & ref_nodes = gyper::graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = gyper::graph.var_nodes;

  SECTION("The nodes should be correctly connected")
  {
    REQUIRE(ref_nodes[0].out_degree() == 2);
    REQUIRE(ref_nodes[0].get_var_index(0) == 0);
    REQUIRE(ref_nodes[0].get_var_index(1) == 1);

    REQUIRE(var_nodes[0].out_degree() == 1);
    REQUIRE(var_nodes[0].get_out_ref_index() == 1);

    REQUIRE(var_nodes[1].out_degree() == 1);
    REQUIRE(var_nodes[1].get_out_ref_index() == 1);

    REQUIRE(ref_nodes[1].out_degree() == 0);
  }

  SECTION("The nodes should have the correct order")
  {
    REQUIRE(ref_nodes[0].get_label().order == 0);
    REQUIRE(var_nodes[0].get_label().order == 36);
    REQUIRE(var_nodes[1].get_label().order == 36);
    REQUIRE(ref_nodes[1].get_label().order == 37);
  }

  SECTION("The nodes should have a label with the correct DNA bases")
  {
    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec("AGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTTT"));
    REQUIRE(var_nodes[0].get_label().dna == gyper::to_vec("C"));
    REQUIRE(var_nodes[1].get_label().dna == gyper::to_vec("G"));
    REQUIRE(ref_nodes[1].get_label().dna == gyper::to_vec("CCCAGGTTTCCCCAGGTTTCCCCTTTGGA"));
  }
}


TEST_CASE("Construct test graph (chr1) but with absolute positions")
{
  using namespace gyper;

  create_graph("/test/data/reference/index_test.fa", "/test/data/reference/index_test.vcf.gz", "chr1", true);

  REQUIRE(graph.ref_nodes.size() == 2);
  REQUIRE(graph.var_nodes.size() == 2);

  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  SECTION("The nodes should be correctly connected")
  {
    REQUIRE(ref_nodes[0].out_degree() == 2);
    REQUIRE(ref_nodes[0].get_var_index(0) == 0);
    REQUIRE(ref_nodes[0].get_var_index(1) == 1);

    REQUIRE(var_nodes[0].out_degree() == 1);
    REQUIRE(var_nodes[0].get_out_ref_index() == 1);

    REQUIRE(var_nodes[1].out_degree() == 1);
    REQUIRE(var_nodes[1].get_out_ref_index() == 1);

    REQUIRE(ref_nodes[1].out_degree() == 0);
  }

  SECTION("The nodes should have the correct order")
  {
    REQUIRE(ref_nodes[0].get_label().order == 0 + 1);
    REQUIRE(var_nodes[0].get_label().order == 36 + 1);
    REQUIRE(var_nodes[1].get_label().order == 36 + 1);
    REQUIRE(ref_nodes[1].get_label().order == 37 + 1);
  }

  SECTION("The nodes should have a label with the correct DNA bases")
  {
    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec("AGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTTT"));
    REQUIRE(var_nodes[0].get_label().dna == gyper::to_vec("C"));
    REQUIRE(var_nodes[1].get_label().dna == gyper::to_vec("G"));
    REQUIRE(ref_nodes[1].get_label().dna == gyper::to_vec("CCCAGGTTTCCCCAGGTTTCCCCTTTGGA"));
  }
}


TEST_CASE("Construct test graph (chr2)")
{
  using namespace gyper;

  create_graph("/test/data/reference/index_test.fa", "/test/data/reference/index_test.vcf.gz", "chr2", false);

  REQUIRE(graph.ref_nodes.size() == 3);
  REQUIRE(graph.var_nodes.size() == 4);

  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  SECTION("The nodes should be correctly connected")
  {
    REQUIRE(ref_nodes[0].out_degree() == 2);
    REQUIRE(ref_nodes[0].get_var_index(0) == 0);
    REQUIRE(ref_nodes[0].get_var_index(1) == 1);

    REQUIRE(var_nodes[0].get_out_ref_index() == 1);
    REQUIRE(var_nodes[1].get_out_ref_index() == 1);

    REQUIRE(ref_nodes[1].get_var_index(0) == 2);
    REQUIRE(ref_nodes[1].get_var_index(1) == 3);

    REQUIRE(var_nodes[2].get_out_ref_index() == 2);
    REQUIRE(var_nodes[3].get_out_ref_index() == 2);

    REQUIRE(ref_nodes[2].out_degree() == 0);
  }

  SECTION("The nodes should have the correct order")
  {
    REQUIRE(ref_nodes[0].get_label().order == 0);
    REQUIRE(var_nodes[0].get_label().order == 1);
    REQUIRE(var_nodes[1].get_label().order == 1);
    REQUIRE(ref_nodes[1].get_label().order == 2);
    REQUIRE(var_nodes[2].get_label().order == 2);
    REQUIRE(var_nodes[3].get_label().order == 2);
    REQUIRE(ref_nodes[2].get_label().order == 3);
  }

  SECTION("The nodes should have a label with the correct DNA bases")
  {
    // CCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTTTGGACCC
    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec("C"));
    REQUIRE(var_nodes[0].get_label().dna == gyper::to_vec("C"));
    REQUIRE(var_nodes[1].get_label().dna == gyper::to_vec("A"));
    REQUIRE(ref_nodes[1].get_label().dna == gyper::to_vec(""));
    REQUIRE(var_nodes[2].get_label().dna == gyper::to_vec("C"));
    REQUIRE(var_nodes[3].get_label().dna == gyper::to_vec("A"));
    REQUIRE(ref_nodes[2].get_label().dna == gyper::to_vec("CAGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTTTGGACCC"));
  }
}


TEST_CASE("Construct test graph (chr3)")
{
  using namespace gyper;

  // AAAACAAAATAAAACAAAATAAAAGAAAACAAAATAAAACAAAATAAAAGAAAACATTATAAAACA
  // chr3 31 rs4 A G,GA

  create_graph("/test/data/reference/index_test.fa", "/test/data/reference/index_test.vcf.gz", "chr3", false);

  REQUIRE(graph.ref_nodes.size() == 2);
  REQUIRE(graph.ref_nodes[0].get_label().dna == gyper::to_vec("AAAACAAAATAAAACAAAATAAAAGAAAAC"));
  REQUIRE(graph.ref_nodes[1].get_label().dna == gyper::to_vec("AAAACAAAATAAAAGAAAACATTATAAAACA"));
  REQUIRE(graph.var_nodes.size() == 3);
  REQUIRE(graph.var_nodes[0].get_label().dna == gyper::to_vec("AAAAT"));
  REQUIRE(graph.var_nodes[1].get_label().dna == gyper::to_vec("GAAAT"));
  REQUIRE(graph.var_nodes[2].get_label().dna == gyper::to_vec("GAAAAT"));

  REQUIRE(graph.actual_poses.size() == 1);
  REQUIRE(graph.actual_poses[0] == 35);
  REQUIRE(graph.ref_reach_poses.size() == 1);
  REQUIRE(graph.ref_reach_poses[0] == 34);
  REQUIRE(std::distance(graph.ref_reach_to_special_pos.begin(), graph.ref_reach_to_special_pos.end()) == 1);
  REQUIRE(graph.ref_reach_to_special_pos.count(34) == 1);
}
