#include <climits>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>

#include <graphtyper/graph/absolute_position.hpp>
#include <graphtyper/graph/graph.hpp>
#include <graphtyper/graph/label.hpp>
#include <graphtyper/graph/var_record.hpp>
#include <graphtyper/utilities/logging.hpp>
#include <graphtyper/utilities/options.hpp>
#include <graphtyper/utilities/type_conversions.hpp>

#include <catch2/catch.hpp>

std::vector<std::vector<char>> get_var_dna(std::vector<gyper::VarNode> const & var_nodes)
{
  std::vector<std::vector<char>> var_dna;

  for (auto const & v : var_nodes)
    var_dna.push_back(v.get_label().dna);

  return var_dna;
}

TEST_CASE("Get the reference sequence of a graph")
{
  using namespace gyper;

  gyper::print_log(gyper::log_severity::debug, "TEST_CASE: Get the reference sequence of a graph.");

  std::vector<char> reference_sequence;
  char testdata[] = "SGTACGEEF";
  reference_sequence.insert(reference_sequence.end(), testdata, testdata + 9);
  std::vector<gyper::VarRecord> records;

  {
    gyper::VarRecord record;
    record.pos = 1;
    record.ref = {'G', 'T', 'A', 'C', 'G'};
    record.alts = {{'G'}};
    records.push_back(record);

    record.pos = 1;
    record.ref = {'G'};
    record.alts = {{'a'}, {'b'}};
    records.push_back(record);

    record.pos = 2;
    record.ref = {'T'};
    record.alts = {{'c'}};
    records.push_back(record);

    record.pos = 4;
    record.ref = {'C'};
    record.alts = {{'d'}};
    records.push_back(record);

    record.pos = 5;
    record.ref = {'G', 'E', 'E'};
    record.alts = {{'G', 'e'}};
    records.push_back(record);
  }

  graph = gyper::Graph(); // use_absolute_positions
  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());

  // Get the entire reference
  {
    REQUIRE(graph.get_all_ref() == gyper::to_vec("SGTACGEEF"));
    REQUIRE(graph.get_ref(0, 999) == gyper::to_vec("SGTACGEEF"));
  }

  // Get the partial reference
  {
    REQUIRE(graph.get_ref(1, 3) == gyper::to_vec("SG"));
    REQUIRE(graph.get_ref(2, 3) == gyper::to_vec("G"));
    REQUIRE(graph.get_ref(2, 9999) == gyper::to_vec("GTACGEEF"));
    REQUIRE(graph.get_ref(5, 8) == gyper::to_vec("CGE"));
    REQUIRE(graph.get_ref(9, 13) == gyper::to_vec("F"));
  }

  // Out of bounds reference
  {
    REQUIRE(graph.get_ref(10, 13) == gyper::to_vec(""));
    REQUIRE(graph.get_ref(9999, 1999992) == gyper::to_vec(""));
  }

  // To is smaller
  {
    REQUIRE(graph.get_ref(4, 1) == gyper::to_vec(""));
  }
}

TEST_CASE("Graph with a reference only.")
{
  using namespace gyper;

  gyper::print_log(gyper::log_severity::debug, "TEST_CASE: Graph with a reference only.");

  std::vector<char> reference_sequence;
  char testdata[] = "ACCGGGAAAA";
  reference_sequence.insert(reference_sequence.end(), testdata, testdata + 10);

  graph = gyper::Graph();
  graph.add_genomic_region(std::move(reference_sequence), std::vector<gyper::VarRecord>(0), gyper::GenomicRegion());
  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  // The graph has the correct size
  {
    REQUIRE(ref_nodes.size() == 1);
    REQUIRE(var_nodes.size() == 0);
  }

  // The only node has the correct order and dna
  {
    REQUIRE(ref_nodes[0].out_degree() == 0);
    REQUIRE(ref_nodes[0].get_label().order == 1);
    REQUIRE(ref_nodes[0].get_label().dna.size() == 10);
    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec("ACCGGGAAAA"));
  }

  // The graph outputs the correct reference sequence
  {
    REQUIRE(graph.get_all_ref() == gyper::to_vec("ACCGGGAAAA"));
  }
}

TEST_CASE("Graph with two variant records.")
{
  using namespace gyper;

  gyper::print_log(gyper::log_severity::debug, "TEST_CASE: Graph with two variant records.");

  std::vector<char> reference_sequence;
  char testdata[] = "ACCGGGAAAA";
  reference_sequence.insert(reference_sequence.end(), testdata, testdata + 10);
  std::vector<gyper::VarRecord> records;

  {
    gyper::VarRecord record;
    record.pos = 3;
    record.ref = {'G'};
    record.alts = {{'G', 'T'}};

    records.push_back(record);

    record.pos = 6;
    record.ref = {'A'};
    record.alts = {{'A', 'T'}, {'G'}};

    records.push_back(record);
  }

  graph = gyper::Graph();
  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());
  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  // The graph should have the correct size
  {
    REQUIRE(ref_nodes.size() == 3);
    REQUIRE(var_nodes.size() == 5);
    REQUIRE(graph.size() == 8);
  }

  // The nodes should be correctly connected
  {
    REQUIRE(ref_nodes[0].out_degree() == 2);
    REQUIRE(ref_nodes[0].get_var_index(0) == 0);
    REQUIRE(ref_nodes[0].get_var_index(1) == 1);

    REQUIRE(var_nodes[0].out_degree() == 1);
    REQUIRE(var_nodes[0].get_out_ref_index() == 1);

    REQUIRE(var_nodes[1].out_degree() == 1);
    REQUIRE(var_nodes[1].get_out_ref_index() == 1);

    REQUIRE(ref_nodes[1].out_degree() == 3);
    REQUIRE(ref_nodes[1].get_var_index(0) == 2);
    REQUIRE(ref_nodes[1].get_var_index(1) == 3);
    REQUIRE(ref_nodes[1].get_var_index(2) == 4);

    REQUIRE(var_nodes[2].out_degree() == 1);
    REQUIRE(var_nodes[2].get_out_ref_index() == 2);

    REQUIRE(var_nodes[3].out_degree() == 1);
    REQUIRE(var_nodes[3].get_out_ref_index() == 2);

    REQUIRE(var_nodes[4].out_degree() == 1);
    REQUIRE(var_nodes[4].get_out_ref_index() == 2);

    REQUIRE(ref_nodes[2].out_degree() == 0);
  }

  // The nodes should have the correct order
  {
    REQUIRE(ref_nodes[0].get_label().order == 0 + 1);
    REQUIRE(var_nodes[0].get_label().order == 3 + 1);
    REQUIRE(var_nodes[1].get_label().order == 3 + 1);
    REQUIRE(ref_nodes[1].get_label().order == 4 + 1);
    REQUIRE(var_nodes[2].get_label().order == 6 + 1);
    REQUIRE(var_nodes[3].get_label().order == 6 + 1);
    REQUIRE(var_nodes[4].get_label().order == 6 + 1);
    REQUIRE(ref_nodes[2].get_label().order == 7 + 1);
  }

  // The nodes should have a label with the correct DNA bases
  {
    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec("ACC"));
    REQUIRE(var_nodes[0].get_label().dna == gyper::to_vec("G"));
    REQUIRE(var_nodes[1].get_label().dna == gyper::to_vec("GT"));
    REQUIRE(ref_nodes[1].get_label().dna == gyper::to_vec("GG"));
    REQUIRE(var_nodes[2].get_label().dna == gyper::to_vec("A"));
    REQUIRE(var_nodes[3].get_label().dna == gyper::to_vec("AT"));
    REQUIRE(var_nodes[4].get_label().dna == gyper::to_vec("G"));
    REQUIRE(ref_nodes[2].get_label().dna == gyper::to_vec("AAA"));
  }
}

TEST_CASE("Graph can start with a variant record.")
{
  using namespace gyper;
  std::vector<char> reference_sequence;
  char testdata[] = "ACCGGGAAAA";
  reference_sequence.insert(reference_sequence.end(), testdata, testdata + 10);
  std::vector<gyper::VarRecord> records;

  {
    gyper::VarRecord record;
    record.pos = 0;
    record.ref = {'A'};
    record.alts = {{'C'}};

    records.push_back(record);

    record.pos = 6;
    record.ref = {'A'};
    record.alts = {{'A', 'T'}, {'G'}};

    records.push_back(record);
  }

  graph = gyper::Graph();
  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());
  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  // The graph should have the correct size
  {
    REQUIRE(ref_nodes.size() == 3);
    REQUIRE(var_nodes.size() == 5);
    REQUIRE(graph.size() == 8);
  }

  // The nodes should be correctly connected
  {
    REQUIRE(ref_nodes[0].out_degree() == 2);
    REQUIRE(var_nodes[0].out_degree() == 1);
    REQUIRE(var_nodes[1].out_degree() == 1);
    REQUIRE(ref_nodes[1].out_degree() == 3);
    REQUIRE(var_nodes[2].out_degree() == 1);
    REQUIRE(var_nodes[3].out_degree() == 1);
    REQUIRE(var_nodes[4].out_degree() == 1);
    REQUIRE(ref_nodes[2].out_degree() == 0);
  }

  // The nodes should have the correct order
  {
    REQUIRE(ref_nodes[0].get_label().order == 0 + 1);
    REQUIRE(var_nodes[0].get_label().order == 0 + 1);
    REQUIRE(var_nodes[1].get_label().order == 0 + 1);
    REQUIRE(ref_nodes[1].get_label().order == 1 + 1);
    REQUIRE(var_nodes[2].get_label().order == 6 + 1);
    REQUIRE(var_nodes[3].get_label().order == 6 + 1);
    REQUIRE(var_nodes[4].get_label().order == 6 + 1);
    REQUIRE(ref_nodes[2].get_label().order == 7 + 1);
  }

  // The nodes should have a label with the correct DNA bases
  {
    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec(""));
    REQUIRE(var_nodes[0].get_label().dna == gyper::to_vec("A"));
    REQUIRE(var_nodes[1].get_label().dna == gyper::to_vec("C"));
    REQUIRE(ref_nodes[1].get_label().dna == gyper::to_vec("CCGGG"));
    REQUIRE(var_nodes[2].get_label().dna == gyper::to_vec("A"));
    REQUIRE(var_nodes[3].get_label().dna == gyper::to_vec("AT"));
    REQUIRE(var_nodes[4].get_label().dna == gyper::to_vec("G"));
    REQUIRE(ref_nodes[2].get_label().dna == gyper::to_vec("AAA"));
  }
}

TEST_CASE("The reference can contain Ns. Note however that variants cannot overlap the N.")
{
  using namespace gyper;

  gyper::print_log(gyper::log_severity::debug,
                   "[",
                   __HERE__,
                   "] TEST_CASE: The reference can contain Ns. ",
                   "Note however that variants cannot overlap the N.");

  std::vector<char> reference_sequence;
  char testdata[] = "ACCGNGAAAA";
  reference_sequence.insert(reference_sequence.end(), testdata, testdata + 10);
  std::vector<gyper::VarRecord> records;

  {
    gyper::VarRecord record;
    record.pos = 3;
    record.ref = {'G'};
    record.alts = {{'G', 'T'}};

    records.push_back(record);

    record.pos = 6;
    record.ref = {'A'};
    record.alts = {{'A', 'T'}, {'G'}};

    records.push_back(record);
  }

  graph = gyper::Graph();
  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());

  REQUIRE(graph.ref_nodes.size() == 3);
  REQUIRE(graph.var_nodes.size() == 5);
  REQUIRE(graph.size() == 8);

  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  // "The nodes should be correctly connected")
  {
    REQUIRE(ref_nodes[0].out_degree() == 2);
    REQUIRE(ref_nodes[0].get_var_index(0) == 0);
    REQUIRE(ref_nodes[0].get_var_index(1) == 1);

    REQUIRE(var_nodes[0].out_degree() == 1);
    REQUIRE(var_nodes[0].get_out_ref_index() == 1);

    REQUIRE(var_nodes[1].out_degree() == 1);
    REQUIRE(var_nodes[1].get_out_ref_index() == 1);

    REQUIRE(ref_nodes[1].out_degree() == 3);
    REQUIRE(ref_nodes[1].get_var_index(0) == 2);
    REQUIRE(ref_nodes[1].get_var_index(1) == 3);
    REQUIRE(ref_nodes[1].get_var_index(2) == 4);

    REQUIRE(var_nodes[2].out_degree() == 1);
    REQUIRE(var_nodes[2].get_out_ref_index() == 2);

    REQUIRE(var_nodes[3].out_degree() == 1);
    REQUIRE(var_nodes[3].get_out_ref_index() == 2);

    REQUIRE(var_nodes[4].out_degree() == 1);
    REQUIRE(var_nodes[4].get_out_ref_index() == 2);

    REQUIRE(ref_nodes[2].out_degree() == 0);
  }

  // The nodes should have the correct order
  {
    REQUIRE(ref_nodes[0].get_label().order == 1);
    REQUIRE(var_nodes[0].get_label().order == 4);
    REQUIRE(var_nodes[1].get_label().order == 4);
    REQUIRE(ref_nodes[1].get_label().order == 5);
    REQUIRE(var_nodes[2].get_label().order == 7);
    REQUIRE(var_nodes[3].get_label().order == 7);
    REQUIRE(var_nodes[4].get_label().order == 7);
    REQUIRE(ref_nodes[2].get_label().order == 8);
  }

  // The nodes should have a label with the correct DNA bases
  {
    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec("ACC"));
    REQUIRE(var_nodes[0].get_label().dna == gyper::to_vec("G"));
    REQUIRE(var_nodes[1].get_label().dna == gyper::to_vec("GT"));
    REQUIRE(ref_nodes[1].get_label().dna == gyper::to_vec("NG"));

    REQUIRE(var_nodes[2].get_label().dna == gyper::to_vec("A"));
    REQUIRE(var_nodes[3].get_label().dna == gyper::to_vec("AT"));
    REQUIRE(var_nodes[4].get_label().dna == gyper::to_vec("G"));

    REQUIRE(ref_nodes[2].get_label().dna == gyper::to_vec("AAA"));
  }
}

TEST_CASE("The reference can start with Ns.")
{
  using namespace gyper;

  gyper::print_log(gyper::log_severity::debug, "[", __HERE__, "] TEST CASE: The reference can start with Ns.");
  std::vector<char> reference_sequence;
  char testdata[] = "NNCGGGAAAA";
  reference_sequence.insert(reference_sequence.end(), testdata, testdata + 10);
  std::vector<gyper::VarRecord> records;

  {
    gyper::VarRecord record;
    record.pos = 3;
    record.ref = {'G'};
    record.alts = {{'G', 'T'}};

    records.push_back(record);

    record.pos = 6;
    record.ref = {'A'};
    record.alts = {{'A', 'T'}, {'G'}};

    records.push_back(record);
  }

  graph = gyper::Graph();
  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());

  REQUIRE(graph.ref_nodes.size() == 3);
  REQUIRE(graph.var_nodes.size() == 5);
  REQUIRE(graph.size() == 8);

  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  // The nodes should be correctly connected
  {
    REQUIRE(ref_nodes[0].out_degree() == 2);
    REQUIRE(ref_nodes[0].get_var_index(0) == 0);

    REQUIRE(var_nodes[0].out_degree() == 1);
    REQUIRE(var_nodes[0].get_out_ref_index() == 1);
    REQUIRE(var_nodes[1].out_degree() == 1);
    REQUIRE(var_nodes[1].get_out_ref_index() == 1);

    REQUIRE(ref_nodes[1].out_degree() == 3);
    REQUIRE(ref_nodes[1].get_var_index(0) == 2);
    REQUIRE(ref_nodes[1].get_var_index(1) == 3);
    REQUIRE(ref_nodes[1].get_var_index(2) == 4);
  }

  // The nodes should have the correct order
  {
    REQUIRE(ref_nodes[0].get_label().order == 1);
    REQUIRE(var_nodes[0].get_label().order == 4);
    REQUIRE(var_nodes[1].get_label().order == 4);
    REQUIRE(ref_nodes[1].get_label().order == 5);
    REQUIRE(var_nodes[2].get_label().order == 7);
    REQUIRE(var_nodes[3].get_label().order == 7);
    REQUIRE(var_nodes[4].get_label().order == 7);
    REQUIRE(ref_nodes[2].get_label().order == 8);
  }

  // The nodes should have a label with the correct DNA bases
  {
    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec("NNC"));
    REQUIRE(var_nodes[0].get_label().dna == gyper::to_vec("G"));
    REQUIRE(var_nodes[1].get_label().dna == gyper::to_vec("GT"));
    REQUIRE(ref_nodes[1].get_label().dna == gyper::to_vec("GG"));
    REQUIRE(var_nodes[2].get_label().dna == gyper::to_vec("A"));
    REQUIRE(var_nodes[3].get_label().dna == gyper::to_vec("AT"));
    REQUIRE(var_nodes[4].get_label().dna == gyper::to_vec("G"));
    REQUIRE(ref_nodes[2].get_label().dna == gyper::to_vec("AAA"));
  }
}

TEST_CASE("We can start at any location of the reference.")
{
  using namespace gyper;

  gyper::print_log(gyper::log_severity::debug,
                   "[",
                   __HERE__,
                   "] TEST CASE: We can start at any location of the reference.");
  std::vector<char> reference_sequence;
  char testdata[] = "CCGGTAAAT";
  reference_sequence.insert(reference_sequence.end(), testdata, testdata + 9);
  std::vector<gyper::VarRecord> records;

  {
    gyper::VarRecord record;
    record.pos = 3;
    record.ref = {'G', 'G'};
    record.alts = {{'G', 'T'}};

    records.push_back(record);

    record.pos = 6;
    record.ref = {'A'};
    record.alts = {{'A', 'T'}, {'G'}};

    records.push_back(record);
  }

  graph = gyper::Graph();
  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion("chr1:2"));

  REQUIRE(graph.ref_nodes.size() == 3);
  REQUIRE(graph.var_nodes.size() == 5);
  REQUIRE(graph.size() == 8);

  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  // The nodes should be correctly connected
  {
    REQUIRE(ref_nodes[0].out_degree() == 2);
    REQUIRE(ref_nodes[0].get_var_index(0) == 0);
    REQUIRE(ref_nodes[0].get_var_index(1) == 1);

    REQUIRE(var_nodes[0].out_degree() == 1);
    REQUIRE(var_nodes[0].get_out_ref_index() == 1);
    REQUIRE(var_nodes[1].out_degree() == 1);
    REQUIRE(var_nodes[1].get_out_ref_index() == 1);

    REQUIRE(ref_nodes[1].out_degree() == 3);
    REQUIRE(ref_nodes[1].get_var_index(0) == 2);
    REQUIRE(ref_nodes[1].get_var_index(1) == 3);
    REQUIRE(ref_nodes[1].get_var_index(2) == 4);

    REQUIRE(var_nodes[2].get_out_ref_index() == 2);
    REQUIRE(var_nodes[3].get_out_ref_index() == 2);
    REQUIRE(var_nodes[4].get_out_ref_index() == 2);

    REQUIRE(ref_nodes[2].out_degree() == 0);
  }

  // The nodes should have the correct order
  {
    REQUIRE(ref_nodes[0].get_label().order == 2);
    REQUIRE(var_nodes[0].get_label().order == 4);
    REQUIRE(var_nodes[1].get_label().order == 4);
    REQUIRE(ref_nodes[1].get_label().order == 6);
    REQUIRE(var_nodes[2].get_label().order == 7);
    REQUIRE(var_nodes[3].get_label().order == 7);
    REQUIRE(ref_nodes[2].get_label().order == 8);
  }

  // The nodes should have a label with the correct DNA bases
  {
    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec("CC"));
    REQUIRE(var_nodes[0].get_label().dna == gyper::to_vec("GG"));
    REQUIRE(var_nodes[1].get_label().dna == gyper::to_vec("GT"));
    REQUIRE(ref_nodes[1].get_label().dna == gyper::to_vec("T"));
    REQUIRE(var_nodes[2].get_label().dna == gyper::to_vec("A"));
    REQUIRE(var_nodes[3].get_label().dna == gyper::to_vec("AT"));
    REQUIRE(ref_nodes[2].get_label().dna == gyper::to_vec("AAT"));
  }
}

TEST_CASE("Variants can overlap", "[graph]")
{
  using namespace gyper;

  gyper::print_log(gyper::log_severity::debug, __HERE__, " TEST CASE: Variants can overlap.");
  std::vector<char> reference_sequence;
  char testdata[] = "ACGGTAA";
  reference_sequence.insert(reference_sequence.end(), testdata, testdata + 7);
  std::vector<gyper::VarRecord> records;

  {
    // ACGGTAA
    // ACTAA
    // ACGATTAA
    gyper::VarRecord record;
    record.pos = 2;
    record.ref = {'G', 'G', 'T'};
    record.alts = {{'T'}};

    records.push_back(record);

    record.pos = 3;
    record.ref = {'G'};
    record.alts = {{'A', 'T'}};

    records.push_back(record);
  }

  graph = gyper::Graph();
  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());
  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  // The graph should have the correct size
  {
    REQUIRE(ref_nodes.size() == 2);
    REQUIRE(var_nodes.size() == 3);
    REQUIRE(graph.size() == 5);
  }

  // The nodes should be correctly connected
  {
    REQUIRE(ref_nodes[0].out_degree() == 3);
    REQUIRE(ref_nodes[0].get_var_index(0) == 0);
    REQUIRE(ref_nodes[0].get_var_index(1) == 1);
    REQUIRE(ref_nodes[0].get_var_index(2) == 2);

    REQUIRE(var_nodes[0].out_degree() == 1);
    REQUIRE(var_nodes[0].get_out_ref_index() == 1);
    REQUIRE(var_nodes[1].out_degree() == 1);
    REQUIRE(var_nodes[1].get_out_ref_index() == 1);
    REQUIRE(var_nodes[1].out_degree() == 1);
    REQUIRE(var_nodes[1].get_out_ref_index() == 1);

    REQUIRE(ref_nodes[1].out_degree() == 0);
  }

  // The nodes should have the correct order
  {
    REQUIRE(ref_nodes[0].get_label().order == 1);
    REQUIRE(var_nodes[0].get_label().order == 3);
    REQUIRE(var_nodes[1].get_label().order == 3);
    REQUIRE(var_nodes[2].get_label().order == 3);
    REQUIRE(ref_nodes[1].get_label().order == 6);
  }

  // The nodes should have a label with the correct DNA bases
  {
    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec("AC"));
    REQUIRE(var_nodes[0].get_label().dna == gyper::to_vec("GGT"));
    REQUIRE(var_nodes[1].get_label().dna == gyper::to_vec("GATT"));
    REQUIRE(var_nodes[2].get_label().dna == gyper::to_vec("T"));
    REQUIRE(ref_nodes[1].get_label().dna == gyper::to_vec("AA"));
  }
}

TEST_CASE("Variants can overlap. Case where the second variant reaches further.")
{
  using namespace gyper;

  gyper::print_log(gyper::log_severity::debug,
                   "[",
                   __HERE__,
                   "] TEST CASE: Variants can overlap. ",
                   "Case where the second variant reaches further.");

  std::vector<char> reference_sequence;
  char testdata[] = "ACGGTAA";
  reference_sequence.insert(reference_sequence.end(), testdata, testdata + 7);
  std::vector<gyper::VarRecord> records;

  {
    // ACGGTAA
    // ACTAA
    // ACGCA
    gyper::VarRecord record;
    record.pos = 2;
    record.ref = {'G', 'G', 'T'};
    record.alts = {{'T'}};

    records.push_back(record);

    record.pos = 3;
    record.ref = {'G', 'T', 'A'};
    record.alts = {{'C'}};

    records.push_back(record);
  }

  graph = gyper::Graph();
  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());
  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  // "The graph should have the correct size"
  {
    REQUIRE(ref_nodes.size() == 2);
    REQUIRE(var_nodes.size() == 3);
    REQUIRE(graph.size() == 5);
  }

  // "The nodes should be correctly connected"
  {
    REQUIRE(ref_nodes[0].out_degree() == 3);
    REQUIRE(ref_nodes[0].get_var_index(0) == 0);
    REQUIRE(ref_nodes[0].get_var_index(1) == 1);
    REQUIRE(ref_nodes[0].get_var_index(2) == 2);

    REQUIRE(var_nodes[0].out_degree() == 1);
    REQUIRE(var_nodes[0].get_out_ref_index() == 1);
    REQUIRE(var_nodes[1].out_degree() == 1);
    REQUIRE(var_nodes[1].get_out_ref_index() == 1);
    REQUIRE(var_nodes[2].out_degree() == 1);
    REQUIRE(var_nodes[2].get_out_ref_index() == 1);

    REQUIRE(ref_nodes[1].out_degree() == 0);
  }

  // "The nodes should have the correct order"
  {
    REQUIRE(ref_nodes[0].get_label().order == 1);
    REQUIRE(var_nodes[0].get_label().order == 3);
    REQUIRE(var_nodes[1].get_label().order == 3);
    REQUIRE(var_nodes[2].get_label().order == 3);
    REQUIRE(ref_nodes[1].get_label().order == 7);
  }

  // "The nodes should have a label with the correct DNA bases"
  {
    // ACGGTAA
    // ACGCA
    // ACTAA
    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec("AC"));
    REQUIRE(var_nodes[0].get_label().dna == gyper::to_vec("GGTA"));
    REQUIRE(var_nodes[1].get_label().dna == gyper::to_vec("GC"));
    REQUIRE(var_nodes[2].get_label().dna == gyper::to_vec("TA"));
    REQUIRE(ref_nodes[1].get_label().dna == gyper::to_vec("A"));
  }
}

/*
TEST_CASE("Two variants next to each other won't overlap")
{
  using namespace gyper;
  Options::instance()->add_all_variants = true;
  std::vector<char> reference_sequence;
  char testdata[] = "ACGCTAA";
  reference_sequence.insert(reference_sequence.end(), testdata, testdata + 7);
  std::vector<gyper::VarRecord> records;

  {
    // ACGCTAA
    // ACTCTAA
    // ACGGTAA
    // ACGTTAA
    gyper::VarRecord record;
    record.pos = 2;
    record.ref = {'G'};
    record.alts = {{'T'}};

    records.push_back(record);

    record.pos = 3;
    record.ref = {'C'};
    record.alts = {{'G'}, {'T'}};

    records.push_back(record);
  }

  graph = gyper::Graph();
  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());
  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  // "The graph should have the correct size"
  {
    REQUIRE(ref_nodes.size() == 2);
    REQUIRE(var_nodes.size() == 6);
    REQUIRE(graph.size() == 8);
  }

  // "The nodes should be correctly connected"
  {
    REQUIRE(ref_nodes[0].out_degree() == 2);
    REQUIRE(ref_nodes[0].get_var_index(0) == 0);
    REQUIRE(ref_nodes[0].get_var_index(1) == 1);

    REQUIRE(var_nodes[0].out_degree() == 1);
    REQUIRE(var_nodes[0].get_out_ref_index() == 1);
    REQUIRE(var_nodes[1].out_degree() == 1);
    REQUIRE(var_nodes[1].get_out_ref_index() == 1);

    REQUIRE(ref_nodes[1].out_degree() == 3);
    REQUIRE(ref_nodes[1].get_var_index(0) == 2);
    REQUIRE(ref_nodes[1].get_var_index(1) == 3);
    REQUIRE(ref_nodes[1].get_var_index(2) == 4);

    REQUIRE(var_nodes[2].get_out_ref_index() == 2);
    REQUIRE(var_nodes[3].get_out_ref_index() == 2);
    REQUIRE(var_nodes[4].get_out_ref_index() == 2);

    REQUIRE(ref_nodes[2].out_degree() == 0);
  }

  // "The nodes should have the correct order"
  {
    REQUIRE(ref_nodes[0].get_label().order == 1);
    REQUIRE(var_nodes[0].get_label().order == 3);
    REQUIRE(var_nodes[1].get_label().order == 3);
    REQUIRE(ref_nodes[1].get_label().order == 4);
    REQUIRE(var_nodes[2].get_label().order == 4);
    REQUIRE(var_nodes[3].get_label().order == 4);
    REQUIRE(var_nodes[4].get_label().order == 4);
    REQUIRE(ref_nodes[2].get_label().order == 5);

  }

  // "The nodes should have a label with the correct DNA bases"
  {
    // ACGCTAA
    // ACTCTAA
    // ACGGTAA
    // ACGTTAA
    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec("AC"));
    REQUIRE(var_nodes[0].get_label().dna == gyper::to_vec("G"));
    REQUIRE(var_nodes[1].get_label().dna == gyper::to_vec("T"));
    REQUIRE(ref_nodes[1].get_label().dna == gyper::to_vec(""));
    REQUIRE(var_nodes[2].get_label().dna == gyper::to_vec("C"));
    REQUIRE(var_nodes[3].get_label().dna == gyper::to_vec("G"));
    REQUIRE(var_nodes[4].get_label().dna == gyper::to_vec("T"));
    REQUIRE(ref_nodes[2].get_label().dna == gyper::to_vec("TAA"));
  }

  Options::instance()->add_all_variants = false;
}
*/

/*
TEST_CASE("Two variants with 3 bp between them")
{
  using namespace gyper;
  Options::instance()->add_all_variants = true;
  std::vector<char> reference_sequence;
  char testdata[] = "ACGCTAA";
  reference_sequence.insert(reference_sequence.end(), testdata, testdata + 7);
  std::vector<gyper::VarRecord> records;

  {
    // ACGCTAA
    // ACTCTAA
    // ACGGTAA
    // ACGTTAA
    gyper::VarRecord record;
    record.pos = 1;
    record.ref = {'C'};
    record.alts = {{'T'}};

    records.push_back(record);

    record.pos = 4;
    record.ref = {'T'};
    record.alts = {{'G'}, {'A'}};

    records.push_back(record);
  }

  graph = gyper::Graph();
  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());
  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  // "The graph should have the correct size"
  {
    REQUIRE(ref_nodes.size() == 2);
    REQUIRE(var_nodes.size() == 6);
    REQUIRE(graph.size() == 8);
  }

  {
    std::vector<std::vector<char> > var_dna = get_var_dna(var_nodes);

    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec("A"));

    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("CGCT")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("CGCG")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("CGCA")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("TGCG")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("TGCA")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("TGCT")) != var_dna.end());

    REQUIRE(ref_nodes[1].get_label().dna == gyper::to_vec("AA"));
  }

  Options::instance()->add_all_variants = false;
}
*/

TEST_CASE(
  "When merging a deletion covering multiple short variants, all combinations of the variants need to be added.",
  "[test2]")
{
  using namespace gyper;

  gyper::print_log(gyper::log_severity::debug,
                   __HERE__,
                   " TEST CASE: ",
                   "When merging a deletion covering multiple short variants, all ",
                   "combinations of the variants need to be added.");

  Options::instance()->add_all_variants = true;
  std::vector<char> reference_sequence;
  char testdata[] = "SSGTAEE";
  reference_sequence.insert(reference_sequence.end(), testdata, testdata + 7);
  std::vector<gyper::VarRecord> records;

  {
    // ACGTGCG reference
    // ACG deletion
    // ACGaGCG snp3a
    // ACGbGCG snp3b
    // ACGTcCG snp4c
    // ACGTdCG snp4d
    // ACGacCG snp3a4c
    // ACGadCG snp3a4d
    // ACGbcCG snp3b4c
    // ACGbdCG snp3b4d
    gyper::VarRecord record;
    record.pos = 2;
    record.ref = {'G', 'T', 'A', 'E', 'E'};
    record.alts = {{'G'}};
    records.push_back(record);

    record.pos = 3;
    record.ref = {'T'};
    record.alts = {{'a'}, {'b'}};
    records.push_back(record);

    record.pos = 4;
    record.ref = {'A'};
    record.alts = {{'c'}, {'d'}};
    records.push_back(record);
  }

  graph = gyper::Graph();
  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());
  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  // "The graph should have the correct size"
  {
    REQUIRE(ref_nodes.size() == 2);
    REQUIRE(var_nodes.size() == 10);
    REQUIRE(graph.size() == 12);
  }

  // "The nodes should be correctly connected"
  {
    REQUIRE(ref_nodes[0].out_degree() == 10);
    REQUIRE(ref_nodes[0].get_var_index(0) == 0);
    REQUIRE(ref_nodes[0].get_var_index(1) == 1);
    REQUIRE(ref_nodes[0].get_var_index(2) == 2);
    REQUIRE(ref_nodes[0].get_var_index(3) == 3);
    REQUIRE(ref_nodes[0].get_var_index(4) == 4);
    REQUIRE(ref_nodes[0].get_var_index(5) == 5);
    REQUIRE(ref_nodes[0].get_var_index(6) == 6);
    REQUIRE(ref_nodes[0].get_var_index(7) == 7);
    REQUIRE(ref_nodes[0].get_var_index(8) == 8);
    REQUIRE(ref_nodes[0].get_var_index(9) == 9);

    REQUIRE(var_nodes[0].get_out_ref_index() == 1);
    REQUIRE(var_nodes[1].get_out_ref_index() == 1);
    REQUIRE(var_nodes[2].get_out_ref_index() == 1);
    REQUIRE(var_nodes[3].get_out_ref_index() == 1);
    REQUIRE(var_nodes[4].get_out_ref_index() == 1);
    REQUIRE(var_nodes[5].get_out_ref_index() == 1);
    REQUIRE(var_nodes[6].get_out_ref_index() == 1);
    REQUIRE(var_nodes[7].get_out_ref_index() == 1);
    REQUIRE(var_nodes[8].get_out_ref_index() == 1);
    REQUIRE(var_nodes[9].get_out_ref_index() == 1);

    REQUIRE(ref_nodes[1].out_degree() == 0);
  }

  // "The nodes should have the correct order"
  {
    REQUIRE(ref_nodes[0].get_label().order == 1);

    REQUIRE(var_nodes[0].get_label().order == 3);
    REQUIRE(var_nodes[1].get_label().order == 3);
    REQUIRE(var_nodes[2].get_label().order == 3);
    REQUIRE(var_nodes[3].get_label().order == 3);
    REQUIRE(var_nodes[4].get_label().order == 3);
    REQUIRE(var_nodes[5].get_label().order == 3);
    REQUIRE(var_nodes[6].get_label().order == 3);
    REQUIRE(var_nodes[7].get_label().order == 3);
    REQUIRE(var_nodes[8].get_label().order == 3);
    REQUIRE(var_nodes[9].get_label().order == 3);

    REQUIRE(ref_nodes[1].get_label().order == 8);
  }

  // "The nodes should have a label with the correct DNA bases"
  {
    std::vector<std::vector<char>> var_dna = get_var_dna(var_nodes);

    // for (auto const & var_node_seq : var_dna)
    //{
    //  BOOST_LOG_TRIVIAL(info) << __HERE__ << " " << std::string(var_node_seq.begin(), var_node_seq.end());
    //}

    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec("SS"));

    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GTAEE")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GTcEE")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GTdEE")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GacEE")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GadEE")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GbcEE")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GbdEE")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GaAEE")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GbAEE")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("G")) != var_dna.end());

    REQUIRE(ref_nodes[1].get_label().dna == gyper::to_vec(""));
  }

  Options::instance()->add_all_variants = false;
}

TEST_CASE("Same as above but with bases in between the variants.", "[test2]")
{
  using namespace gyper;

  gyper::print_log(gyper::log_severity::debug,
                   "[",
                   __HERE__,
                   "] TEST CASE: ",
                   "Same as above but with bases in between the variants.");

  Options::instance()->add_all_variants = true;
  std::vector<char> reference_sequence;
  char testdata[] = "GTACE";
  reference_sequence.insert(reference_sequence.end(), testdata, testdata + 5);
  std::vector<gyper::VarRecord> records;

  {
    // GTACE
    gyper::VarRecord record;
    record.pos = 0;
    record.ref = {'G', 'T', 'A', 'C'};
    record.alts = {{'G'}};
    records.push_back(record);

    record.pos = 0;
    record.ref = {'G'};
    record.alts = {{'a'}, {'b'}};
    records.push_back(record);

    record.pos = 2;
    record.ref = {'A'};
    record.alts = {{'c'}, {'d'}};
    records.push_back(record);
  }

  graph = gyper::Graph();
  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());
  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  // "The nodes should have a label with the correct DNA bases"
  {
    std::vector<std::vector<char>> var_dna = get_var_dna(var_nodes);

    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec(""));

    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GTAC")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GTcC")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GTdC")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("aTcC")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("aTdC")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("bTcC")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("bTdC")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("aTAC")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("bTAC")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("G")) != var_dna.end());

    REQUIRE(ref_nodes[1].get_label().dna == gyper::to_vec("E"));
  }

  // "The graph should have the correct size"
  {
    REQUIRE(ref_nodes.size() == 2);
    REQUIRE(var_nodes.size() == 10);
    REQUIRE(graph.size() == 12);
  }

  // "The nodes should be correctly connected"
  {
    REQUIRE(ref_nodes[0].out_degree() == 10);
    REQUIRE(ref_nodes[0].get_var_index(0) == 0);
    REQUIRE(ref_nodes[0].get_var_index(1) == 1);
    REQUIRE(ref_nodes[0].get_var_index(2) == 2);
    REQUIRE(ref_nodes[0].get_var_index(3) == 3);
    REQUIRE(ref_nodes[0].get_var_index(4) == 4);
    REQUIRE(ref_nodes[0].get_var_index(5) == 5);
    REQUIRE(ref_nodes[0].get_var_index(6) == 6);
    REQUIRE(ref_nodes[0].get_var_index(7) == 7);
    REQUIRE(ref_nodes[0].get_var_index(8) == 8);
    REQUIRE(ref_nodes[0].get_var_index(9) == 9);

    REQUIRE(var_nodes[0].get_out_ref_index() == 1);
    REQUIRE(var_nodes[1].get_out_ref_index() == 1);
    REQUIRE(var_nodes[2].get_out_ref_index() == 1);
    REQUIRE(var_nodes[3].get_out_ref_index() == 1);
    REQUIRE(var_nodes[4].get_out_ref_index() == 1);
    REQUIRE(var_nodes[5].get_out_ref_index() == 1);
    REQUIRE(var_nodes[6].get_out_ref_index() == 1);
    REQUIRE(var_nodes[7].get_out_ref_index() == 1);
    REQUIRE(var_nodes[8].get_out_ref_index() == 1);
    REQUIRE(var_nodes[9].get_out_ref_index() == 1);

    REQUIRE(ref_nodes[1].out_degree() == 0);
  }

  // "The nodes should have the correct order"
  {
    REQUIRE(ref_nodes[0].get_label().order == 1);

    REQUIRE(var_nodes[0].get_label().order == 1);
    REQUIRE(var_nodes[1].get_label().order == 1);
    REQUIRE(var_nodes[2].get_label().order == 1);
    REQUIRE(var_nodes[3].get_label().order == 1);
    REQUIRE(var_nodes[4].get_label().order == 1);
    REQUIRE(var_nodes[5].get_label().order == 1);
    REQUIRE(var_nodes[6].get_label().order == 1);
    REQUIRE(var_nodes[7].get_label().order == 1);
    REQUIRE(var_nodes[8].get_label().order == 1);
    REQUIRE(var_nodes[9].get_label().order == 1);

    REQUIRE(ref_nodes[1].get_label().order == 5);
  }

  Options::instance()->add_all_variants = false;
}

TEST_CASE("Four variants joined", "[test2]")
{
  using namespace gyper;

  gyper::print_log(gyper::log_severity::debug, "[", __HERE__, "] TEST CASE: Four variants joined.");

  Options::instance()->add_all_variants = true;
  std::vector<char> reference_sequence;
  char testdata[] = "SGTACGE";
  reference_sequence.insert(reference_sequence.end(), testdata, testdata + 7);
  std::vector<gyper::VarRecord> records;

  {
    // GTACE
    gyper::VarRecord record;
    record.pos = 1;
    record.ref = {'G', 'T', 'A', 'C', 'G'};
    record.alts = {{'G'}};
    records.push_back(record);

    record.pos = 1;
    record.ref = {'G'};
    record.alts = {{'a'}, {'b'}};
    records.push_back(record);

    record.pos = 2;
    record.ref = {'T'};
    record.alts = {{'c'}};
    records.push_back(record);

    record.pos = 4;
    record.ref = {'C'};
    record.alts = {{'d'}};
    records.push_back(record);
  }

  graph = gyper::Graph();
  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());
  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  // The nodes should have a label with the correct DNA bases
  {
    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec("S"));
    std::vector<std::vector<char>> var_dna = get_var_dna(var_nodes);
    std::sort(var_dna.begin(), var_dna.end());

    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("G")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GTACG")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GTAdG")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GcACG")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GcAdG")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("aTACG")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("aTAdG")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("acACG")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("acAdG")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("bTAdG")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("bcACG")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("bcAdG")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("bTACG")) != var_dna.end());

    REQUIRE(ref_nodes[1].get_label().dna == gyper::to_vec("E"));
  }

  // The graph should have the correct size
  {
    REQUIRE(ref_nodes.size() == 2);
    REQUIRE(var_nodes.size() == 13);
    REQUIRE(graph.size() == 15);
  }

  // The nodes should be correctly connected
  {
    REQUIRE(ref_nodes[0].out_degree() == 13);
    REQUIRE(ref_nodes[0].get_var_index(0) == 0);
    REQUIRE(ref_nodes[0].get_var_index(1) == 1);
    REQUIRE(ref_nodes[0].get_var_index(2) == 2);
    REQUIRE(ref_nodes[0].get_var_index(3) == 3);
    REQUIRE(ref_nodes[0].get_var_index(4) == 4);
    REQUIRE(ref_nodes[0].get_var_index(5) == 5);
    REQUIRE(ref_nodes[0].get_var_index(6) == 6);
    REQUIRE(ref_nodes[0].get_var_index(7) == 7);
    REQUIRE(ref_nodes[0].get_var_index(8) == 8);
    REQUIRE(ref_nodes[0].get_var_index(9) == 9);
    REQUIRE(ref_nodes[0].get_var_index(10) == 10);
    REQUIRE(ref_nodes[0].get_var_index(11) == 11);
    REQUIRE(ref_nodes[0].get_var_index(12) == 12);

    REQUIRE(var_nodes[0].get_out_ref_index() == 1);
    REQUIRE(var_nodes[1].get_out_ref_index() == 1);
    REQUIRE(var_nodes[2].get_out_ref_index() == 1);
    REQUIRE(var_nodes[3].get_out_ref_index() == 1);
    REQUIRE(var_nodes[4].get_out_ref_index() == 1);
    REQUIRE(var_nodes[5].get_out_ref_index() == 1);
    REQUIRE(var_nodes[6].get_out_ref_index() == 1);
    REQUIRE(var_nodes[7].get_out_ref_index() == 1);
    REQUIRE(var_nodes[8].get_out_ref_index() == 1);
    REQUIRE(var_nodes[9].get_out_ref_index() == 1);
    REQUIRE(var_nodes[10].get_out_ref_index() == 1);
    REQUIRE(var_nodes[11].get_out_ref_index() == 1);
    REQUIRE(var_nodes[12].get_out_ref_index() == 1);

    REQUIRE(ref_nodes[1].out_degree() == 0);
  }

  // The nodes should have the correct order
  {
    REQUIRE(ref_nodes[0].get_label().order == 1);

    for (auto const & v : var_nodes)
      REQUIRE(v.get_label().order == 2);

    REQUIRE(ref_nodes[1].get_label().order == 7);
  }

  Options::instance()->add_all_variants = false;
}

TEST_CASE("Variants of any number can be joined, here 3 are tested.", "[test]")
{
  using namespace gyper;

  gyper::print_log(gyper::log_severity::debug,
                   __HERE__,
                   " TEST CASE: Variants of any number can be joined, here 3 are tested.");

  Options::instance()->add_all_variants = true;
  std::vector<char> reference_sequence;
  char testdata[] = "SGTACGEEF";
  reference_sequence.insert(reference_sequence.end(), testdata, testdata + 9);
  std::vector<gyper::VarRecord> records;

  {
    // GTACE
    gyper::VarRecord record;
    record.pos = 1;
    record.ref = {'G', 'T', 'A', 'C', 'G'};
    record.alts = {{'G'}};
    records.push_back(record);

    record.pos = 1;
    record.ref = {'G'};
    record.alts = {{'a'}, {'b'}};
    records.push_back(record);

    record.pos = 2;
    record.ref = {'T'};
    record.alts = {{'c'}};
    records.push_back(record);

    record.pos = 4;
    record.ref = {'C'};
    record.alts = {{'d'}};
    records.push_back(record);

    record.pos = 5;
    record.ref = {'G', 'E', 'E'};
    record.alts = {{'G', 'e'}};
    records.push_back(record);
  }

  graph = gyper::Graph();
  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());
  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  // The nodes should have a label with the correct DNA bases"
  {
    std::vector<std::vector<char>> var_dna = get_var_dna(var_nodes);
    std::sort(var_dna.begin(), var_dna.end());

    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec("S"));

    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GEE")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GTACGEE")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GTACGe")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GTAdGEE")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GTAdGe")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GcACGEE")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GcACGe")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GcAdGEE")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GcAdGe")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("aTACGEE")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("aTACGe")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("aTAdGEE")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("aTAdGe")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("acACGEE")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("acACGe")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("acAdGEE")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("acAdGe")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("bTACGEE")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("bTACGe")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("bTAdGEE")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("bTAdGe")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("bcACGEE")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("bcACGe")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("bcAdGEE")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("bcAdGe")) != var_dna.end());

    REQUIRE(ref_nodes[1].get_label().dna == gyper::to_vec("F"));
  }

  // The graph should have the correct size
  {
    REQUIRE(ref_nodes.size() == 2);
    REQUIRE(var_nodes.size() == 25);
    REQUIRE(graph.size() == 27);
  }

  // The nodes should be correctly connected
  {
    REQUIRE(ref_nodes[0].out_degree() == 25);
    REQUIRE(ref_nodes[0].get_var_index(0) == 0);
    REQUIRE(ref_nodes[0].get_var_index(24) == 24);

    REQUIRE(var_nodes[0].get_out_ref_index() == 1);
    REQUIRE(var_nodes[24].get_out_ref_index() == 1);

    REQUIRE(ref_nodes[1].out_degree() == 0);
  }

  Options::instance()->add_all_variants = false;
}

// chr6  29890476  145_H GGCT  AGCG,GGCG 0 . .
// chr6  29890479  . TTGAAAGGTGAGACCTTGGGGGGCCTGATGTGTGGGGGA GTGAAAGGTGAGACCTTGGGGGGCCTGATGTGTGAGGGG . PASS  .
// chr6  29890513  . GGGGA AGGGG . PASS  .
// chr6  29890517  146_H ATGTTGGGGGGGAACAGTGG  GTGTTGGGGGGGAACAGTGG,GTGTTGGGGGGGAACAGTG  0 . .

// TEST_CASE("Check if the above creates the graph correctly")
//{
//  using namespace gyper;
//
//  BOOST_LOG_TRIVIAL(debug) << "[" << __HERE__ << "] TEST CASE: Check if the above creates the graph correctly";
//
//  Options::instance()->add_all_variants = true;
//  std::vector<char> reference_sequence;
//  char testdata[] = "CG1xGAaTGTTGGGGGGGAACAGT";
//  reference_sequence.insert(reference_sequence.end(), testdata, testdata + 24);
//
//  std::vector<gyper::VarRecord> records_step1;
//  std::vector<gyper::VarRecord> records_step2;
//
//  {
//    gyper::VarRecord record;
//    record.pos = 1;
//    record.ref = gyper::to_vec("G1");
//    record.alts = {gyper::to_vec("A2"), gyper::to_vec("G3")};
//    records_step1.push_back(record);
//    records_step2.push_back(record);
//
//    record.pos = 2;
//    record.ref = gyper::to_vec("1xGA");
//    record.alts = {gyper::to_vec("3yAG")};
//    records_step1.push_back(record);
//    records_step2.push_back(record);
//
//    record.pos = 4;
//    record.ref = gyper::to_vec("GA");
//    record.alts = {gyper::to_vec("AG")};
//    records_step2.push_back(record);
//  }
//
//  // Step 1: Simpler graph
//  {
//    graph = gyper::Graph();
//    graph.add_genomic_region(std::vector<char>(reference_sequence),
//                             std::move(records_step1),
//                             gyper::GenomicRegion());
//
//    std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
//    std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;
//
//    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec("C"));
//
//    REQUIRE(var_nodes[0].get_label().dna == gyper::to_vec("G1xGA"));
//    REQUIRE(var_nodes[1].get_label().dna == gyper::to_vec("A2xGA"));
//    REQUIRE(var_nodes[2].get_label().dna == gyper::to_vec("A3yAG"));
//    REQUIRE(var_nodes[2].get_label().dna == gyper::to_vec("G3xGA"));
//    REQUIRE(var_nodes[3].get_label().dna == gyper::to_vec("G3yAG"));
//
//    REQUIRE(var_nodes.size() == 4);
//  }
//
//  // Step 2: Simpler graph
//  {
//    graph = gyper::Graph();
//    graph.add_genomic_region(std::vector<char>(reference_sequence),
//                             std::move(records_step2),
//                             gyper::GenomicRegion());
//
//    std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
//    std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;
//
//    std::vector<std::vector<char> > var_dna;
//
//    for (auto const & v : var_nodes)
//    {
//      var_dna.push_back(v.get_label().dna);
//    }
//
//    REQUIRE(ref_nodes[0].get_label().dna ==   gyper::to_vec("C"));
//
//    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("G1xGA")) != var_dna.end());
//    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("G1xAG")) != var_dna.end());
//    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("G3yAG")) != var_dna.end());
//    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("A2xGA")) != var_dna.end());
//    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("A2xAG")) != var_dna.end());
//    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("G3xGA")) != var_dna.end());
//    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("G3xAG")) != var_dna.end());
//    //REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("A3yAG")) != var_dna.end());
//
//    REQUIRE(var_nodes.size() == 7);
//  }
//
//  Options::instance()->add_all_variants = false;
//}

TEST_CASE("Variant overlapping a N on the reference genome")
{
  using namespace gyper;
  Options::instance()->add_all_variants = true;
  std::vector<char> reference_sequence;
  char testdata[] = "GCTGCGGCGGGCGTCGCGGCCGCCCCCGGGGAGCCCGGCGGGCGCCGGCGCGNCCCCCCCCCCACCCCACGTCTCGTCGCGCGCGC";
  reference_sequence.insert(reference_sequence.end(), testdata, testdata + 86);

  {
    gyper::print_log(gyper::log_severity::info, __HERE__, " The reference allele has an N.");
    std::vector<gyper::VarRecord> records;

    {
      gyper::VarRecord record;
      record.pos = 51;
      record.ref = gyper::to_vec("GN");
      record.alts = {gyper::to_vec("GA")};
      records.push_back(record);
    }

    graph = gyper::Graph();
    graph.add_genomic_region(std::vector<char>(reference_sequence), std::move(records), gyper::GenomicRegion());

    // Here, we assume that nothing is added
    std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
    std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

    REQUIRE(ref_nodes.size() == 1);
    REQUIRE(ref_nodes[0].get_label().dna == reference_sequence);
    REQUIRE(var_nodes.size() == 0);
  }

  {
    gyper::print_log(gyper::log_severity::info, __HERE__, " An alternative allele has an N.");

    std::vector<gyper::VarRecord> records;

    {
      gyper::VarRecord record;
      record.pos = 51;
      record.ref = gyper::to_vec("G");
      record.alts = {gyper::to_vec("GN"), gyper::to_vec("GA")};
      records.push_back(record);
    }

    graph = gyper::Graph();
    graph.add_genomic_region(std::vector<char>(reference_sequence), std::move(records), gyper::GenomicRegion());

    std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
    std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

    REQUIRE(ref_nodes.size() == 2);
    REQUIRE(var_nodes.size() == 2);

    REQUIRE(var_nodes[0].get_label().dna == gyper::to_vec("G"));
    REQUIRE(var_nodes[1].get_label().dna == gyper::to_vec("GA"));
  }

  // If all alternative allele have a N, we remove the variant
  {
    std::vector<gyper::VarRecord> records;

    {
      gyper::VarRecord record;
      record.pos = 51;
      record.ref = gyper::to_vec("G");
      record.alts = {gyper::to_vec("GN"), gyper::to_vec("GNN")};
      records.push_back(record);
    }

    graph = gyper::Graph();
    graph.add_genomic_region(std::vector<char>(reference_sequence), std::move(records), gyper::GenomicRegion());

    // Here, we assume that nothing is added
    std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
    std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

    REQUIRE(ref_nodes.size() == 1);
    REQUIRE(ref_nodes[0].get_label().dna == reference_sequence);
    REQUIRE(var_nodes.size() == 0);
  }

  Options::instance()->add_all_variants = false;
}

TEST_CASE("Prior test for the next")
{
  using namespace gyper;
  Options::instance()->add_all_variants = false;
  std::vector<char> reference_sequence;
  char testdata[] = "GTTCAATG";
  reference_sequence.insert(reference_sequence.end(), testdata, testdata + 8);

  std::vector<gyper::VarRecord> records;

  {
    gyper::VarRecord record;
    record.pos = 2;
    record.ref = gyper::to_vec("TC");
    record.alts = {gyper::to_vec("T")};
    records.push_back(record);

    record.pos = 4;
    record.ref = gyper::to_vec("A");
    record.alts = {gyper::to_vec("T")};
    records.push_back(record);
  }

  graph = gyper::Graph();
  graph.add_genomic_region(std::vector<char>(reference_sequence), std::move(records), gyper::GenomicRegion());

  // Here, we assume that nothing is added
  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  // The number of nodes is correct
  {
    REQUIRE(ref_nodes.size() == 3);
    REQUIRE(var_nodes.size() == 4);
  }

  // We have the correct labels on the reference nodes
  {
    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec("GT"));
    REQUIRE(ref_nodes[1].get_label().dna == gyper::to_vec(""));
    REQUIRE(ref_nodes[2].get_label().dna == gyper::to_vec("ATG"));
  }

  // We have the correct labels on the variant nodes"
  {
    REQUIRE(var_nodes[0].get_label().dna == gyper::to_vec("TC"));
    REQUIRE(var_nodes[1].get_label().dna == gyper::to_vec("T"));
    REQUIRE(var_nodes[2].get_label().dna == gyper::to_vec("A"));
    REQUIRE(var_nodes[3].get_label().dna == gyper::to_vec("T"));
  }

  Options::instance()->add_all_variants = true;
}

TEST_CASE("Merge one path should check if we can remove the suffix of a variant before merging them")
{
  using namespace gyper;
  Options::instance()->add_all_variants = false;
  std::vector<char> reference_sequence;
  std::string testdata = "STAAAAAATF";
  reference_sequence.insert(reference_sequence.end(), testdata.begin(), testdata.end());

  std::vector<gyper::VarRecord> records;

  {
    gyper::VarRecord record;
    record.pos = 1;
    record.ref = gyper::to_vec("TAAAAAAT");   // TA
    record.alts = {gyper::to_vec("TAAAAAT")}; // T
    records.push_back(record);

    record.pos = 7;
    record.ref = gyper::to_vec("A");
    record.alts = {gyper::to_vec("T")};
    records.push_back(record);
  }

  graph = gyper::Graph();
  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());

  // Here, we assume that nothing is added
  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  REQUIRE(ref_nodes.size() == 2);
  REQUIRE(var_nodes.size() == 4);

  // We have the correct labels on the reference nodes
  {
    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec("S"));
    REQUIRE(ref_nodes[1].get_label().dna == gyper::to_vec("TF"));
  }

  // We have the correct labels on the variant nodes
  {
    std::vector<std::vector<char>> var_dna = get_var_dna(var_nodes);

    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec("S"));

    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("TAAAAAA")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("TAAAAA")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("TAAAAAT")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("TAAAAT")) != var_dna.end());
  }
}

TEST_CASE("Merge one path works with connected indel+SNP")
{
  using namespace gyper;
  Options::instance()->add_all_variants = false;
  std::vector<char> reference_sequence;
  char testdata[] = "STAAF";
  reference_sequence.insert(reference_sequence.end(), testdata, testdata + 5);

  std::vector<gyper::VarRecord> records;

  {
    gyper::VarRecord record;
    record.pos = 2;
    record.ref = gyper::to_vec("AA");
    record.alts = {gyper::to_vec("T")};
    records.push_back(record);

    record.pos = 3;
    record.ref = gyper::to_vec("A");
    record.alts = {gyper::to_vec("T")};
    records.push_back(record);
  }

  graph = gyper::Graph();
  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());

  // Here, we assume that nothing is added
  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  // The number of nodes is corrent
  {
    REQUIRE(ref_nodes.size() == 2);
    REQUIRE(var_nodes.size() == 3);
  }

  // The ref nodes have the correct DNA bases
  {
    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec("ST"));
    REQUIRE(ref_nodes[1].get_label().dna == gyper::to_vec("F"));
  }

  // The var nodes have the correct DNA bases
  {
    REQUIRE(var_nodes[0].get_label().dna == gyper::to_vec("AA"));
    REQUIRE(var_nodes[1].get_label().dna == gyper::to_vec("AT"));
    REQUIRE(var_nodes[2].get_label().dna == gyper::to_vec("T"));
  }
}

TEST_CASE("Merge path works with 3 pairs of connected SNPs")
{
  using namespace gyper;
  Options::instance()->add_all_variants = false;
  std::vector<char> reference_sequence;
  char testdata[] = "STAAAF";
  reference_sequence.insert(reference_sequence.end(), testdata, testdata + 6);

  std::vector<gyper::VarRecord> records;

  {
    gyper::VarRecord record;
    record.pos = 2;
    record.ref = gyper::to_vec("A");
    record.alts = {gyper::to_vec("T")};
    records.push_back(record);

    record.pos = 3;
    record.ref = gyper::to_vec("A");
    record.alts = {gyper::to_vec("T")};
    records.push_back(record);

    record.pos = 4;
    record.ref = gyper::to_vec("A");
    record.alts = {gyper::to_vec("T")};
    records.push_back(record);
  }

  graph = gyper::Graph();
  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());

  // Here, we assume that nothing is added
  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  // The number of nodes is corrent
  {
    REQUIRE(ref_nodes.size() == 4);
    REQUIRE(var_nodes.size() == 6);
  }

  // The ref nodes have the correct DNA bases
  {
    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec("ST"));
    REQUIRE(ref_nodes[1].get_label().dna == gyper::to_vec(""));
    REQUIRE(ref_nodes[2].get_label().dna == gyper::to_vec(""));
    REQUIRE(ref_nodes[3].get_label().dna == gyper::to_vec("F"));
  }

  // The var nodes have the correct DNA bases
  {
    REQUIRE(var_nodes[0].get_label().dna == gyper::to_vec("A"));
    REQUIRE(var_nodes[1].get_label().dna == gyper::to_vec("T"));
    REQUIRE(var_nodes[2].get_label().dna == gyper::to_vec("A"));
    REQUIRE(var_nodes[3].get_label().dna == gyper::to_vec("T"));
    REQUIRE(var_nodes[4].get_label().dna == gyper::to_vec("A"));
    REQUIRE(var_nodes[5].get_label().dna == gyper::to_vec("T"));
  }
}

TEST_CASE("Two overlapping indels")
{
  using namespace gyper;
  Options::instance()->add_all_variants = false;
  std::vector<char> reference_sequence;
  char testdata[] = "TGCAAATCTCATATATATATATATATATATATATATATATATATATTTTTTTTTTTTTTTTTTTTTTTTTA";
  reference_sequence.insert(reference_sequence.end(), testdata, testdata + 71);

  std::vector<gyper::VarRecord> records;

  {
    gyper::VarRecord record;
    record.pos = 30;
    record.ref = gyper::to_vec("ATATATATATATATATTTTTTTTTTTT");
    record.alts = {gyper::to_vec("A")};
    records.push_back(record);

    record.pos = 38;
    record.ref = gyper::to_vec("ATATATATTTTTTTTTTT");
    record.alts = {gyper::to_vec("A")};
    records.push_back(record);
  }

  graph = gyper::Graph();
  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());

  // Here, we assume that nothing is added
  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  // The number of nodes is corrent
  {
    REQUIRE(ref_nodes.size() == 2);
    REQUIRE(var_nodes.size() == 3);
  }

  // The ref nodes have the correct DNA bases
  {
    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec("TGCAAATCTCATATATATATATATATATAT"));
    REQUIRE(ref_nodes[1].get_label().dna == gyper::to_vec("TTTTTTTTTTTTTA"));
  }

  // The var nodes have the correct DNA bases
  {
    REQUIRE(var_nodes[0].get_label().dna == gyper::to_vec("ATATATATATATATATTTTTTTTTTTT"));
    REQUIRE(var_nodes[1].get_label().dna == gyper::to_vec("A"));
    REQUIRE(var_nodes[2].get_label().dna == gyper::to_vec("ATATATATAT"));
  }
}

TEST_CASE("Two deletions and one of them overlaps SNPs")
{
  using namespace gyper;

  Options::instance()->add_all_variants = false;
  std::vector<char> reference_sequence;
  std::string testdata = "SGTATATAGCTGCCGCCGTTTTTATTACCGGGGGTAGTAGTAGTAGCGCAGAGGTTTTAGAGGGCF";
  reference_sequence.insert(reference_sequence.end(), testdata.begin(), testdata.end());
  std::vector<gyper::VarRecord> records;

  {
    gyper::VarRecord record;
    record.pos = 1;
    record.ref = {'G', 'T'};
    record.alts = {{'G'}};
    records.push_back(record);

    record.pos = 1;
    record.ref = {'G', 'T', 'A', 'T', 'A', 'T', 'A', 'G', 'C', 'T', 'G', 'C', 'C', 'G', 'C', 'C', 'G', 'T', 'T', 'T'};
    record.alts = {{'G'}};
    records.push_back(record);

    record.pos = 9;
    record.ref = {'C'};
    record.alts = {{'a'}, {'b'}};
    records.push_back(record);

    record.pos = 11;
    record.ref = {'G'};
    record.alts = {{'c'}, {'d'}};
    records.push_back(record);
  }

  graph = gyper::Graph();

  {
    Contig new_contig;
    new_contig.name = "chr1";
    new_contig.length = 100000;
    graph.contigs.push_back(std::move(new_contig));

    absolute_pos.calculate_offsets(graph.contigs);
  }

  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion("chr1"));

  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  // The graph should have the correct size
  {
    REQUIRE(ref_nodes.size() == 2);
    REQUIRE(var_nodes.size() == 11);
    REQUIRE(graph.size() == 13);
  }

  // The nodes should be correctly connected
  {
    REQUIRE(ref_nodes[0].out_degree() == 11);

    for (decltype(ref_nodes[0].get_var_index(0)) i = 0; i < 11; ++i)
      REQUIRE(ref_nodes[0].get_var_index(i) == i);

    for (auto const & v : var_nodes)
      REQUIRE(v.get_out_ref_index() == 1);

    REQUIRE(ref_nodes[1].out_degree() == 0);
  }

  // The nodes should have the correct order
  {
    REQUIRE(ref_nodes[0].get_label().order == 1);

    for (auto & v : var_nodes)
      REQUIRE(v.get_label().order == 2);

    REQUIRE(ref_nodes[1].get_label().order == 22);
  }

  // The nodes should have a label with the correct DNA bases
  {
    std::vector<std::vector<char>> var_dna = get_var_dna(var_nodes);

    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec("S"));

    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("G")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GATATAGCTGCCGCCGTTT")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GATATAGCTcCCGCCGTTT")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GATATAGCTdCCGCCGTTT")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GATATAGaTGCCGCCGTTT")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GATATAGbTGCCGCCGTTT")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GTATATAGCTcCCGCCGTTT")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GTATATAGCTdCCGCCGTTT")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GTATATAGaTGCCGCCGTTT")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GTATATAGbTGCCGCCGTTT")) != var_dna.end());

    REQUIRE(ref_nodes[1].get_label().dna == gyper::to_vec("TTATTACCGGGGGTAGTAGTAGTAGCGCAGAGGTTTTAGAGGGCF"));
  }
}

TEST_CASE("Two deletions and one of them overlaps SNPs and an insertion")
{
  using namespace gyper;
  Options::instance()->add_all_variants = false;
  std::vector<char> reference_sequence;
  std::string testdata = "SGTATATAGCTGCCGCCGTTTTTATTACCGGGGGTAGTAGTAGTAGCGCAGAGGTTTTAGAGGGCF";
  reference_sequence.insert(reference_sequence.end(), testdata.begin(), testdata.end());
  std::vector<gyper::VarRecord> records;

  {
    gyper::VarRecord record;
    record.pos = 1;
    record.ref = {'G', 'T'};
    record.alts = {{'G'}};
    records.push_back(record);

    record.pos = 1;
    record.ref = {'G', 'T', 'A', 'T', 'A', 'T', 'A', 'G', 'C', 'T', 'G', 'C', 'C', 'G', 'C', 'C', 'G', 'T', 'T', 'T'};
    record.alts = {{'G'}};
    records.push_back(record);

    record.pos = 9;
    record.ref = {'C'};
    record.alts = {{'a'}, {'b'}};
    records.push_back(record);

    record.pos = 13;
    record.ref = {'C'};
    record.alts = {{'c'}};
    records.push_back(record);

    record.pos = 14;
    record.ref = {'G'};
    record.alts = {{'d', 'e'}};
    records.push_back(record);
  }

  graph = gyper::Graph();

  {
    Contig new_contig;
    new_contig.name = "chr1";
    new_contig.length = 100000;
    graph.contigs.push_back(std::move(new_contig));

    absolute_pos.calculate_offsets(graph.contigs);
  }

  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion("chr1"));
  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  // The graph should have the correct size
  {
    // graph.print();
    REQUIRE(ref_nodes.size() == 2);
    REQUIRE(var_nodes.size() == 19);
    REQUIRE(graph.size() == 21);
  }

  // The nodes should have the correct order
  {
    REQUIRE(ref_nodes[0].get_label().order == 1);

    for (auto & v : var_nodes)
      REQUIRE(v.get_label().order == 2);
  }

  // The nodes should have a label with the correct DNA bases
  {
    std::vector<std::vector<char>> var_dna = get_var_dna(var_nodes);

    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec("S"));

    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("G")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GATATAGCTGCCGCCGTTT")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GATATAGCTGCCdeCCGTTT")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GATATAGCTGCcGCCGTTT")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GATATAGaTGCCGCCGTTT")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GATATAGaTGCCdeCCGTTT")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GATATAGbTGCCGCCGTTT")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GATATAGbTGCCdeCCGTTT")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GTATATAGCTGCCdeCCGTTT")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GTATATAGCTGCcGCCGTTT")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GTATATAGaTGCCGCCGTTT")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GTATATAGaTGCCdeCCGTTT")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GTATATAGbTGCCGCCGTTT")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("GTATATAGbTGCCdeCCGTTT")) != var_dna.end());

    REQUIRE(ref_nodes[1].get_label().dna == gyper::to_vec("TTATTACCGGGGGTAGTAGTAGTAGCGCAGAGGTTTTAGAGGGCF"));
  }
}

TEST_CASE("We cant have two events that sum up to the reference")
{
  using namespace gyper;
  Options::instance()->add_all_variants = true;
  std::vector<char> reference_sequence;
  std::string testdata = "TTACTTTTTTAA";
  reference_sequence.insert(reference_sequence.end(), testdata.begin(), testdata.end());
  std::vector<gyper::VarRecord> records;

  {
    gyper::VarRecord record;
    record.pos = 3;
    record.ref = {'C'};
    record.alts = {{'C', 'T'}};
    records.push_back(record);

    record.pos = 7;
    record.ref = {'T', 'T'};
    record.alts = {{'T'}};
    records.push_back(record);
  }

  graph = gyper::Graph();

  {
    Contig new_contig;
    new_contig.name = "chr1";
    new_contig.length = 100000;
    graph.contigs.push_back(std::move(new_contig));

    absolute_pos.calculate_offsets(graph.contigs);
  }

  gyper::GenomicRegion genomic_region("chr1");
  genomic_region.check_if_var_records_match_reference_genome(records, reference_sequence);

  for (auto & var_record : records)
  {
    genomic_region.add_reference_to_record_if_they_have_a_matching_prefix(var_record, reference_sequence);
    // BOOST_LOG_TRIVIAL(debug) << __HERE__ << var_record.to_string();
  }

  graph.add_genomic_region(std::move(reference_sequence), std::move(records), std::move(genomic_region));
  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  // The nodes should have a label with the correct DNA bases
  {
    std::vector<std::vector<char>> var_dna = get_var_dna(var_nodes);

    for (auto const & var_node_seq : var_dna)
    {
      gyper::print_log(gyper::log_severity::debug,
                       __HERE__,
                       " ",
                       std::string(var_node_seq.begin(), var_node_seq.end()));
    }

    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("CT")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("C")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("CTT")) != var_dna.end());
  }

  // The graph should have the correct size
  {
    REQUIRE(ref_nodes.size() == 2);
    REQUIRE(var_nodes.size() == 3);
  }

  Options::instance()->add_all_variants = false;
}

TEST_CASE("anti events test case")
{
  using namespace gyper;
  Options::instance()->add_all_variants = true;
  std::vector<char> reference_sequence;
  std::string testdata = "TTACTTTATAAATTACTCAGTCTCGGGTATGTCC";
  reference_sequence.insert(reference_sequence.end(), testdata.begin(), testdata.end());
  std::vector<gyper::VarRecord> records;

  {
    gyper::VarRecord record;
    record.pos = 18;
    record.ref = {'A', 'G', 'T', 'C'};
    record.alts = {{'A', 'G'}};
    record.alts[0].anti_events.emplace(2);
    record.alts[0].anti_events.emplace(3);
    records.push_back(record);

    record.pos = 20;
    record.ref = {'T'};
    record.alts = {{'A'}};
    record.alts[0].events.emplace(2);
    record.alts[0].anti_events.emplace(3);
    records.push_back(record);

    record.pos = 21;
    record.ref = {'C'};
    record.alts = {{'T'}};
    record.alts[0].events.emplace(3);
    records.push_back(record);
  }

  graph = gyper::Graph();

  {
    Contig new_contig;
    new_contig.name = "chr1";
    new_contig.length = 100000;
    graph.contigs.push_back(std::move(new_contig));

    absolute_pos.calculate_offsets(graph.contigs);
  }

  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion("chr1"));
  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  // The nodes should have a label with the correct DNA bases
  {
    std::vector<std::vector<char>> var_dna = get_var_dna(var_nodes);

    for (auto const & var_node_seq : var_dna)
    {
      gyper::print_log(gyper::log_severity::info, __HERE__, " ", std::string(var_node_seq.begin(), var_node_seq.end()));
    }

    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("AG")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("AGTC")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("AGAC")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("AGTT")) != var_dna.end());
  }

  // The graph should have the correct size
  {
    REQUIRE(ref_nodes.size() == 2);
    REQUIRE(var_nodes.size() == 4);
  }

  Options::instance()->add_all_variants = false;
}

TEST_CASE("anti events test case 2 - more complex test")
{
  using namespace gyper;
  Options::instance()->add_all_variants = true;
  std::vector<char> reference_sequence;
  std::string testdata = "TCTATTTTTTTTTTTTTTTTTTTTTTGA";
  reference_sequence.insert(reference_sequence.end(), testdata.begin(), testdata.end());
  std::vector<gyper::VarRecord> records;

  {
    gyper::VarRecord record;

    record.pos = 3;
    record.ref = {'A'};
    record.alts = {{'A', 'T', 'T', 'T'}};
    record.alts[0].events.emplace(3);
    record.alts[0].anti_events.emplace(4);
    record.alts[0].anti_events.emplace(5);
    record.alts[0].anti_events.emplace(6);
    record.alts[0].anti_events.emplace(7);
    record.alts[0].anti_events.emplace(8);
    records.push_back(record);

    record.pos = 11;
    record.ref = {'T'};
    record.alts = {{'T', 'A'}};
    record.alts[0].events.emplace(4);
    record.alts[0].anti_events.emplace(5);
    record.alts[0].anti_events.emplace(6);
    record.alts[0].anti_events.emplace(7);
    record.alts[0].anti_events.emplace(8);
    records.push_back(record);

    record.pos = 15;
    record.ref = {'T'};
    record.alts = {{'C'}};
    record.alts[0].events.emplace(5);
    record.alts[0].anti_events.emplace(6);
    record.alts[0].anti_events.emplace(7);
    record.alts[0].anti_events.emplace(8);
    records.push_back(record);

    record.pos = 24;
    record.ref = {'T'};
    record.alts = {{'T', 'T', 'G'}};
    record.alts[0].events.emplace(6);
    record.alts[0].anti_events.emplace(7);
    record.alts[0].anti_events.emplace(8);
    records.push_back(record);

    record.pos = 25;
    record.ref = {'T'};
    record.alts = {{'T', 'T', 'T', 'G'}};
    record.alts[0].events.emplace(7);
    record.alts[0].anti_events.emplace(8);
    records.push_back(record);

    record.pos = 26;
    record.ref = {'G'};
    record.alts = {{'T'}};
    record.alts[0].events.emplace(8);
    records.push_back(record);
  }

  graph = gyper::Graph();

  {
    Contig new_contig;
    new_contig.name = "chr1";
    new_contig.length = 100000;
    graph.contigs.push_back(std::move(new_contig));

    absolute_pos.calculate_offsets(graph.contigs);
  }

  gyper::GenomicRegion genomic_region("chr1");

  // This is needed to prohibit a combination of alt. alleles being able to be equal to the reference
  for (auto & var_record : records)
  {
    genomic_region.add_reference_to_record_if_they_have_a_matching_prefix(var_record, reference_sequence);
  }

  graph.add_genomic_region(std::move(reference_sequence), std::move(records), std::move(genomic_region));
  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  // The nodes should have a label with the correct DNA bases
  {
    std::vector<std::vector<char>> var_dna = get_var_dna(var_nodes);

    for (auto const & var_node_seq : var_dna)
    {
      gyper::print_log(gyper::log_severity::info, __HERE__, " ", std::string(var_node_seq.begin(), var_node_seq.end()));
    }

    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("ATTTTTTTTTTTTTTTTTTTTTTG")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("ATTTTTTTTTTTTTTTTTTTTTTT")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("ATTTTTTTTTTTTTTTTTTTTTTTTGG")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("ATTTTTTTTTTTTTTTTTTTTTTGTG")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("ATTTTTTTTTTTCTTTTTTTTTTG")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("ATTTTTTTTATTTTTTTTTTTTTTG")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("ATTTTTTTTTTTTTTTTTTTTTTTTTG")) != var_dna.end());
  }

  // The graph should have the correct size
  {
    REQUIRE(ref_nodes.size() == 2);
    REQUIRE(var_nodes.size() == 7);
  }

  Options::instance()->add_all_variants = false;
}

TEST_CASE("parity events test case")
{
  using namespace gyper;
  Options::instance()->add_all_variants = true;
  std::vector<char> reference_sequence;
  std::string testdata = "TCTATTTTTTTTTTTTTTTTTTTTTTGA";
  reference_sequence.insert(reference_sequence.end(), testdata.begin(), testdata.end());
  std::vector<gyper::VarRecord> records;

  {
    gyper::VarRecord record;

    record.pos = 3;
    record.ref = {'A'};
    record.alts = {{'A', 'T', 'T'}};
    record.ref.events = {-2};
    record.ref.anti_events = {4};
    record.alts[0].events = {2};
    record.alts[0].anti_events = {3, -4};
    records.push_back(record);
    record.clear();

    record.pos = 3;
    record.ref = {'A'};
    record.alts = {{'A', 'T', 'T', 'T'}};
    record.ref.events = {-3};
    record.ref.anti_events = {};
    record.alts[0].events = {3};
    record.alts[0].anti_events = {4};
    records.push_back(record);
    record.clear();

    record.pos = 11;
    record.ref = {'T'};
    record.alts = {{'T', 'A'}};
    record.ref.events = {-4};
    record.alts[0].events = {4};
    records.push_back(record);
    record.clear();
  }

  graph = gyper::Graph();

  {
    Contig new_contig;
    new_contig.name = "chr1";
    new_contig.length = 100000;
    graph.contigs.push_back(std::move(new_contig));

    absolute_pos.calculate_offsets(graph.contigs);
  }

  gyper::GenomicRegion genomic_region("chr1");

  // This is needed to prohibit a combination of alt. alleles being able to be equal to the reference
  for (auto & var_record : records)
  {
    genomic_region.add_reference_to_record_if_they_have_a_matching_prefix(var_record, reference_sequence);
  }

  graph.add_genomic_region(std::move(reference_sequence), std::move(records), std::move(genomic_region));
  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  // The nodes should have a label with the correct DNA bases
  {
    std::vector<std::vector<char>> var_dna = get_var_dna(var_nodes);

    for (auto const & var_node_seq : var_dna)
    {
      gyper::print_log(gyper::log_severity::debug,
                       __HERE__,
                       " ",
                       std::string(var_node_seq.begin(), var_node_seq.end()));
    }

    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("ATTTTTTTT")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("ATTTTTTTT")) == var_dna.begin()); // ref

    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("ATTTTTTTTTTT")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("ATTTTTTTTA")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("ATTTTTTTTTTA")) != var_dna.end());
  }

  // The graph should have the correct size
  {
    REQUIRE(ref_nodes.size() == 2);
    REQUIRE(var_nodes.size() == 4);
  }

  Options::instance()->add_all_variants = false;
}

TEST_CASE("parity events test case 2 - snps next to each other")
{
  using namespace gyper;
  Options::instance()->add_all_variants = true;
  std::vector<char> reference_sequence;
  std::string testdata = "TCTCAGA";
  reference_sequence.insert(reference_sequence.end(), testdata.begin(), testdata.end());
  std::vector<gyper::VarRecord> records;

  {
    gyper::VarRecord record;

    record.pos = 3;
    record.ref = {'C'};
    record.alts = {{'T'}};
    record.ref.events.emplace(-1);
    record.ref.anti_events.emplace(2);
    record.ref.anti_events.emplace(3);
    record.alts[0].events.emplace(1);
    record.alts[0].anti_events.emplace(-2);
    record.alts[0].anti_events.emplace(-3);
    records.push_back(record);
    record.clear();

    record.pos = 4;
    record.ref = {'A'};
    record.alts = {{'G'}};
    record.ref.events.emplace(-2);
    record.ref.anti_events.emplace(3);
    record.alts[0].events.emplace(2);
    record.alts[0].anti_events.emplace(-3);
    records.push_back(record);
    record.clear();

    record.pos = 5;
    record.ref = {'G'};
    record.alts = {{'A'}};
    record.ref.events.emplace(-3);
    record.alts[0].events.emplace(3);
    records.push_back(record);
    record.clear();
  }

  graph = gyper::Graph();

  {
    Contig new_contig;
    new_contig.name = "chr1";
    new_contig.length = 100000;
    graph.contigs.push_back(std::move(new_contig));

    absolute_pos.calculate_offsets(graph.contigs);
  }

  gyper::GenomicRegion genomic_region("chr1");

  // This is needed to prohibit a combination of alt. alleles being able to be equal to the reference
  for (auto & var_record : records)
  {
    genomic_region.add_reference_to_record_if_they_have_a_matching_prefix(var_record, reference_sequence);
  }

  graph.add_genomic_region(std::move(reference_sequence), std::move(records), std::move(genomic_region));
  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  // graph.print();

  // The graph should have the correct size
  {
    REQUIRE(ref_nodes.size() == 2);
    REQUIRE(var_nodes.size() == 2);
  }

  {
    {
      VarNode const & ref_var_node = var_nodes[0];
      REQUIRE(ref_var_node.get_label().dna == gyper::to_vec("CAG"));
      REQUIRE(ref_var_node.events.count(-1));
      REQUIRE(ref_var_node.events.count(-2));
      REQUIRE(ref_var_node.events.count(-3));
      REQUIRE(ref_var_node.events.size() == 3);
      REQUIRE(ref_var_node.anti_events.size() == 2);
      REQUIRE(ref_var_node.anti_events.count(2));
      REQUIRE(ref_var_node.anti_events.count(3));
    }

    {
      VarNode const & var_node = var_nodes[1];
      REQUIRE(var_node.get_label().dna == gyper::to_vec("TGA"));
      REQUIRE(var_node.events.size() == 3);
      REQUIRE(var_node.events.count(1));
      REQUIRE(var_node.events.count(2));
      REQUIRE(var_node.events.count(3));
      REQUIRE(var_node.anti_events.size() == 2);
      REQUIRE(var_node.anti_events.count(-2));
      REQUIRE(var_node.anti_events.count(-3));
    }
  }

  Options::instance()->add_all_variants = false;
}
