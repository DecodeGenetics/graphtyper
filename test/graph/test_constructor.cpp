#include <catch.hpp>

#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>

#include <cereal/archives/binary.hpp>

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/graph.hpp>
#include <graphtyper/graph/constructor.hpp>
#include <graphtyper/graph/graph_serialization.hpp>
#include <graphtyper/utilities/type_conversions.hpp>
#include <graphtyper/utilities/options.hpp>
#include <graphtyper/constants.hpp>

#include "../help_functions.hpp" // create_test_graph


TEST_CASE("Construct test graph (chr1) (not with absolute positions)")
{
  gyper::print_log(gyper::log_severity::debug, __HERE__, " TEST CASE: Construct test graph (chr1)");

  create_test_graph("/test/data/reference/index_test.fa", "/test/data/reference/index_test.vcf.gz", "chr1", false);

  REQUIRE(gyper::graph.ref_nodes.size() == 2);
  REQUIRE(gyper::graph.var_nodes.size() == 2);

  std::vector<gyper::RefNode> const & ref_nodes = gyper::graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = gyper::graph.var_nodes;

  // The nodes should be correctly connected
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

  // The nodes should have the correct order
  {
    REQUIRE(ref_nodes[0].get_label().order == 1);
    REQUIRE(var_nodes[0].get_label().order == 37);
    REQUIRE(var_nodes[1].get_label().order == 37);
    REQUIRE(ref_nodes[1].get_label().order == 38);
  }

  // The nodes should have a label with the correct DNA bases
  {
    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec("AGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTTT"));
    REQUIRE(var_nodes[0].get_label().dna == gyper::to_vec("C"));
    REQUIRE(var_nodes[1].get_label().dna == gyper::to_vec("G"));
    REQUIRE(ref_nodes[1].get_label().dna == gyper::to_vec("CCCAGGTTTCCCCAGGTTTCCCCTTTGGA"));
  }
}


TEST_CASE("Construct test graph (chr1)")
{
  using namespace gyper;

  gyper::print_log(gyper::log_severity::debug, __HERE__, " TEST CASE: Construct test graph (chr1) but with absolute positions.");
  create_test_graph("/test/data/reference/index_test.fa", "/test/data/reference/index_test.vcf.gz", "chr1", true);

  REQUIRE(graph.ref_nodes.size() == 2);
  REQUIRE(graph.var_nodes.size() == 2);

  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  // SECTION("The nodes should be correctly connected")
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

  // SECTION("The nodes should have the correct order")
  {
    REQUIRE(ref_nodes[0].get_label().order == 1);
    REQUIRE(var_nodes[0].get_label().order == 37);
    REQUIRE(var_nodes[1].get_label().order == 37);
    REQUIRE(ref_nodes[1].get_label().order == 38);
  }

  // SECTION("The nodes should have a label with the correct DNA bases")
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

  Options::instance()->add_all_variants = false;
  gyper::print_log(gyper::log_severity::debug, __HERE__, " Construct test graph (chr2).");
  create_test_graph("/test/data/reference/index_test.fa", "/test/data/reference/index_test.vcf.gz", "chr2", true);

  REQUIRE(graph.ref_nodes.size() == 3);
  REQUIRE(graph.var_nodes.size() == 4);

  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  // SECTION("The nodes should be correctly connected")
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

  // SECTION("The nodes should have the correct order")
  {
    REQUIRE(ref_nodes[0].get_label().order == 1);
    REQUIRE(var_nodes[0].get_label().order == 2);
    REQUIRE(var_nodes[1].get_label().order == 2);
    REQUIRE(ref_nodes[1].get_label().order == 3);
    REQUIRE(var_nodes[2].get_label().order == 3);
    REQUIRE(var_nodes[3].get_label().order == 3);
    REQUIRE(ref_nodes[2].get_label().order == 4);
  }

  // SECTION("The nodes should have a label with the correct DNA bases")
  {
    // CCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTTTGGACCC
    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec("C"));
    REQUIRE(var_nodes[0].get_label().dna == gyper::to_vec("C"));
    REQUIRE(var_nodes[1].get_label().dna == gyper::to_vec("A"));
    REQUIRE(ref_nodes[1].get_label().dna == gyper::to_vec(""));
    REQUIRE(var_nodes[2].get_label().dna == gyper::to_vec("C"));
    REQUIRE(var_nodes[3].get_label().dna == gyper::to_vec("A"));
    REQUIRE(ref_nodes[2].get_label().dna ==
            gyper::to_vec("CAGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTTTGGACCC"));
  }
}


TEST_CASE("Construct test graph (chr3)")
{
  using namespace gyper;

  // AAAACAAAATAAAACAAAATAAAAGAAAACAAAATAAAACAAAATAAAAGAAAACATTATAAAACA
  // chr3 31 rs4 A G,GA

  gyper::print_log(gyper::log_severity::debug, "[", __HERE__, "] TEST CASE: Construct test graph (chr3).");
  create_test_graph("/test/data/reference/index_test.fa", "/test/data/reference/index_test.vcf.gz", "chr3", true);

  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  // SECTION("Nodes are correctly connected")
  {
    REQUIRE(ref_nodes[0].out_degree() == 3);
    REQUIRE(ref_nodes[0].get_var_index(0) == 0);
    REQUIRE(ref_nodes[0].get_var_index(1) == 1);
    REQUIRE(ref_nodes[0].get_var_index(2) == 2);

    REQUIRE(var_nodes[0].get_out_ref_index() == 1);
    REQUIRE(var_nodes[1].get_out_ref_index() == 1);
    REQUIRE(var_nodes[2].get_out_ref_index() == 1);
  }

  // SECTION("Nodes have the correct order")
  {
    REQUIRE(ref_nodes[0].get_label().order == 1);
    REQUIRE(var_nodes[0].get_label().order == 31);
    REQUIRE(var_nodes[1].get_label().order == 31);
    REQUIRE(var_nodes[2].get_label().order == 31);
    REQUIRE(ref_nodes[1].get_label().order == 32);
  }

  // SECTION("Nodes have the correct bases")
  {
    REQUIRE(graph.ref_nodes.size() == 2);
    REQUIRE(graph.ref_nodes[0].get_label().dna == gyper::to_vec("AAAACAAAATAAAACAAAATAAAAGAAAAC"));
    REQUIRE(graph.ref_nodes[1].get_label().dna == gyper::to_vec("AAATAAAACAAAATAAAAGAAAACATTATAAAACA"));
    REQUIRE(graph.var_nodes.size() == 3);
    REQUIRE(graph.var_nodes[0].get_label().dna == gyper::to_vec("A"));
    REQUIRE(graph.var_nodes[1].get_label().dna == gyper::to_vec("G"));
    REQUIRE(graph.var_nodes[2].get_label().dna == gyper::to_vec("GA"));

    REQUIRE(graph.actual_poses.size() == 1);
    REQUIRE(graph.actual_poses[0] == 32);
    REQUIRE(graph.ref_reach_poses.size() == 1);
    REQUIRE(graph.ref_reach_poses[0] == 31);
    REQUIRE(std::distance(graph.ref_reach_to_special_pos.begin(), graph.ref_reach_to_special_pos.end()) == 1);
    REQUIRE(graph.ref_reach_to_special_pos.count(31) == 1);
  }
}


TEST_CASE("Construct test graph (chr8) in a region that fully overlaps only a second indel")
{
  using namespace gyper;
  Options::instance()->add_all_variants = false;

  // TGCAAATCTCATATATATATATATATATATATATATATATATATATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCCAA
  //chr8 31 ATATATATATATATATTTTTTTTTTTT,A
  //chr8 39 ATATATATTTTTTTTTTT,A

  gyper::print_log(gyper::log_severity::debug, __HERE__, " TEST CASE: Construct test graph (chr8) in a region "
                          , "that fully overlaps only a second indel");

  create_test_graph("/test/data/reference/index_test.fa", "/test/data/reference/index_test.vcf.gz", "chr8:1-56", true);

  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  {
    REQUIRE(ref_nodes.size() == 2);
    REQUIRE(var_nodes.size() == 2);
  }

  // Nodes are correctly connected
  {
    REQUIRE(ref_nodes[0].out_degree() == 2);
    REQUIRE(ref_nodes[0].get_var_index(0) == 0);
    REQUIRE(ref_nodes[0].get_var_index(1) == 1);
    REQUIRE(ref_nodes[1].out_degree() == 0);

    REQUIRE(var_nodes[0].get_out_ref_index() == 1);
    REQUIRE(var_nodes[1].get_out_ref_index() == 1);
  }

  // Nodes have the correct order
  {
    REQUIRE(ref_nodes[0].get_label().order == 1);
    REQUIRE(var_nodes[0].get_label().order == 39);
    REQUIRE(var_nodes[1].get_label().order == 39);
    REQUIRE(ref_nodes[1].get_label().order == 57);
  }

  // Nodes have the correct bases
  {
    REQUIRE(graph.ref_nodes.size() == 2);
    REQUIRE(graph.ref_nodes[0].get_label().dna == gyper::to_vec("TGCAAATCTCATATATATATATATATATATATATATAT"));
    REQUIRE(graph.ref_nodes[1].get_label().dna == gyper::to_vec(""));
    REQUIRE(graph.var_nodes.size() == 2);
    REQUIRE(graph.var_nodes[0].get_label().dna == gyper::to_vec("ATATATATTTTTTTTTTT"));
    REQUIRE(graph.var_nodes[1].get_label().dna == gyper::to_vec("A"));

    REQUIRE(graph.actual_poses.size() == 0);
  }
}


/*
TEST_CASE("Construct test graph with SV deletion (chr5)")
{
  using namespace gyper;
  using gyper::to_vec;

  // 70A 70C 70G 70T
  // chr5 70 A70C SVTYPE=DEL,SVSIZE=70

  create_test_graph("/test/data/reference/index_test.fa", "/test/data/reference/index_test.vcf.gz", "chr5", false);

  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  REQUIRE(ref_nodes.size() == 2);
  REQUIRE(var_nodes.size() == 2);

  SECTION("Nodes are correctly connected")
  {
    REQUIRE(ref_nodes[0].out_degree() == 2);
    REQUIRE(ref_nodes[0].get_var_index(0) == 0);
    REQUIRE(ref_nodes[0].get_var_index(1) == 1);

    REQUIRE(var_nodes[0].get_out_ref_index() == 1);
    REQUIRE(var_nodes[1].get_out_ref_index() == 1);
  }

  SECTION("The nodes have the correct bases")
  {
    REQUIRE(ref_nodes[0].get_label().dna == to_vec("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    REQUIRE(var_nodes[0].get_label().dna == to_vec("A"));
    REQUIRE(var_nodes[1].get_label().dna == to_vec("AGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT<SV:0000000<"));
    REQUIRE(ref_nodes[1].get_label().dna == to_vec("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
  }

  SECTION("The nodes should have the correct order")
  {
    REQUIRE(ref_nodes[0].get_label().order == 0);
    REQUIRE(var_nodes[0].get_label().order == 69);
    REQUIRE(var_nodes[1].get_label().order == 69);
    REQUIRE(ref_nodes[1].get_label().order == 70);
  }
}


TEST_CASE("Construct test graph with SV duplication (chr6)")
{
  using namespace gyper;

  create_test_graph("/test/data/reference/index_test.fa", "/test/data/reference/index_test.vcf.gz", "chr6", false);

  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  REQUIRE(ref_nodes.size() == 3);
  REQUIRE(var_nodes.size() == 4);

  SECTION("Nodes are correctly connected")
  {
    REQUIRE(ref_nodes[0].out_degree() == 2);
    REQUIRE(ref_nodes[0].get_var_index(0) == 0);
    REQUIRE(ref_nodes[0].get_var_index(1) == 1);

    REQUIRE(var_nodes[0].get_out_ref_index() == 1);
    REQUIRE(var_nodes[1].get_out_ref_index() == 1);

    REQUIRE(ref_nodes[1].out_degree() == 2);
    REQUIRE(ref_nodes[2].out_degree() == 0);
  }

  SECTION("The nodes have the correct bases")
  {
    REQUIRE(ref_nodes[0].get_label().dna == to_vec("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCC"));
    REQUIRE(var_nodes[0].get_label().dna == to_vec("C"));
    REQUIRE(var_nodes[1].get_label().dna == to_vec("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTT<SV:0000000>AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCC"));
    std::vector<char> const & ref_dna = ref_nodes[1].get_label().dna;
    REQUIRE(ref_dna.size() >= 66);
    REQUIRE(std::vector<char>(ref_dna.begin(), ref_dna.begin() + 66) == to_vec("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCG"));
  }

  SECTION("The nodes should have the correct order")
  {
    REQUIRE(ref_nodes[0].get_label().order == 0);
    REQUIRE(var_nodes[0].get_label().order == 74);
    REQUIRE(var_nodes[1].get_label().order == 74);
    REQUIRE(ref_nodes[1].get_label().order == 75);
  }
}


TEST_CASE("Construct test graph with SV duplication (chr7)")
{
  using namespace gyper;

  create_test_graph("/test/data/reference/index_test.fa", "/test/data/reference/index_test.vcf.gz", "chr7", false);

  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  REQUIRE(ref_nodes.size() == 2);
  REQUIRE(var_nodes.size() == 3);

  SECTION("Nodes are correctly connected")
  {
    REQUIRE(ref_nodes[0].out_degree() == 3);
    REQUIRE(ref_nodes[0].get_var_index(0) == 0);
    REQUIRE(ref_nodes[0].get_var_index(1) == 1);
    REQUIRE(ref_nodes[0].get_var_index(2) == 2);

    REQUIRE(var_nodes[0].get_out_ref_index() == 1);
    REQUIRE(var_nodes[1].get_out_ref_index() == 1);
  }

  SECTION("The nodes have the correct bases")
  {
    REQUIRE(ref_nodes[0].get_label().dna == to_vec("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    REQUIRE(var_nodes[0].get_label().dna == to_vec("A"));
    //REQUIRE(var_nodes[1].get_label().dna == to_vec("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTT<SV:0000000>AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCC"));
    std::vector<char> const & ref_dna = ref_nodes[1].get_label().dna;
    REQUIRE(ref_dna.size() >= 71);
    REQUIRE(std::vector<char>(ref_dna.begin(), ref_dna.begin() + 71) == to_vec("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCG"));
  }

  //SECTION("The nodes should have the correct order")
  //{
  //  REQUIRE(ref_nodes[0].get_label().order == 0);
  //  REQUIRE(var_nodes[0].get_label().order == 74);
  //  REQUIRE(var_nodes[1].get_label().order == 74);
  //  REQUIRE(ref_nodes[1].get_label().order == 75);
  //}
}
*/


/*
TEST_CASE("Construct test graph (chromosome 1, without chr in front)")
{
  using namespace gyper;

  // AAAACAAAATAAAACAAAATAAAAGAAAACAAAATAAAACAAAATAAAAGAAAACATTATAAAACA

  create_test_graph("/test/data/reference/index_test_b37.fa", "/test/data/reference/index_test_b37.vcf.gz", "1", true);

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
    REQUIRE(ref_nodes[0].get_label().order == 0 + 1);
    REQUIRE(var_nodes[0].get_label().order == 36 + 1);
    REQUIRE(var_nodes[1].get_label().order == 36 + 1);
    REQUIRE(ref_nodes[1].get_label().order == 37 + 1);
  }

  SECTION("The nodes should have a label with the correct DNA bases")
  {
    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec("AGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTTT"));
    REQUIRE(var_nodes[0].get_label().dna == gyper::to_vec("C"));
    REQUIRE(var_nodes[1].get_label().dna == gyper::to_vec("A"));
    REQUIRE(ref_nodes[1].get_label().dna == gyper::to_vec("CCCAGGTTTCCCCAGGTTTCCCCTTTGGA"));
  }
}
*/

TEST_CASE("Construct test graph with anti events (chr9)")
{
  using namespace gyper;
  Options::instance()->add_all_variants = true;

  gyper::print_log(gyper::log_severity::debug, __HERE__, " TEST CASE: Construct test graph (chr9).");
  create_test_graph("/test/data/reference/index_test.fa", "/test/data/reference/index_test.vcf.gz", "chr9", true);

  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  {
    REQUIRE(ref_nodes.size() == 3);
    REQUIRE(var_nodes.size() == 4);
  }

  {
    REQUIRE(var_nodes[0].events.size() == 1);
    REQUIRE(var_nodes[0].events.count(-1) == 1);
    REQUIRE(var_nodes[1].events.size() == 1);
    REQUIRE(var_nodes[1].events.count(1) == 1);
    REQUIRE(var_nodes[2].events.size() == 1);
    REQUIRE(var_nodes[2].events.count(-2) == 1);
    REQUIRE(var_nodes[3].events.size() == 1);
    REQUIRE(var_nodes[3].events.count(2) == 1);
  }

  {
    REQUIRE(var_nodes[0].anti_events.size() == 0);
    REQUIRE(var_nodes[1].anti_events.size() == 1);
    REQUIRE(var_nodes[1].anti_events.count(2) == 1);
    REQUIRE(var_nodes[2].anti_events.size() == 0);
    REQUIRE(var_nodes[3].anti_events.size() == 0);
  }

  Options::instance()->add_all_variants = false;
}


TEST_CASE("Construct test graph with anti events (chr10)")
{
  using namespace gyper;
  Options::instance()->add_all_variants = true;

  gyper::print_log(gyper::log_severity::debug, __HERE__, " TEST CASE: Construct test graph (chr10).");
  create_test_graph("/test/data/reference/index_test.fa", "/test/data/reference/index_test.vcf.gz", "chr10", true);

  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  {
    REQUIRE(ref_nodes.size() == 3);
    REQUIRE(var_nodes.size() == 4);
  }

  {
    REQUIRE(var_nodes[0].events.size() == 1);
    REQUIRE(var_nodes[0].events.count(-1) == 1);
    REQUIRE(var_nodes[1].events.size() == 1);
    REQUIRE(var_nodes[1].events.count(1) == 1);
    REQUIRE(var_nodes[2].events.size() == 1);
    REQUIRE(var_nodes[2].events.count(-2) == 1);
    REQUIRE(var_nodes[3].events.size() == 1);
    REQUIRE(var_nodes[3].events.count(2) == 1);
  }

  {
    REQUIRE(var_nodes[0].anti_events.size() == 1);
    REQUIRE(var_nodes[0].anti_events.count(2) == 1);
    REQUIRE(var_nodes[1].anti_events.size() == 0);
    REQUIRE(var_nodes[2].anti_events.size() == 0);
    REQUIRE(var_nodes[3].anti_events.size() == 0);
  }

  Options::instance()->add_all_variants = false;
}


TEST_CASE("Construct test graph with anti events (chr11)")
{
  using namespace gyper;
  Options::instance()->add_all_variants = true;

  gyper::print_log(gyper::log_severity::debug, __HERE__, " TEST CASE: Construct test graph (chr11).");
  create_test_graph("/test/data/reference/index_test.fa", "/test/data/reference/index_test.vcf.gz", "chr11", true);

  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  {
    REQUIRE(ref_nodes.size() == 2);
    REQUIRE(var_nodes.size() == 6);
  }

  Options::instance()->add_all_variants = false;
}
