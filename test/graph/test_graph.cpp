#include <string>
#include <vector>
#include <stdio.h>
#include <climits>
#include <cstdio>

#include <iostream>
#include <fstream>

#include <graphtyper/graph/graph.hpp>
#include <graphtyper/graph/label.hpp>
#include <graphtyper/graph/var_record.hpp>
#include <graphtyper/utilities/type_conversions.hpp>
#include <graphtyper/utilities/options.hpp>

#include <catch.hpp>


std::vector<std::vector<char> >
get_var_dna(std::vector<gyper::VarNode> const & var_nodes)
{
  std::vector<std::vector<char> > var_dna;

  for (auto const & v : var_nodes)
    var_dna.push_back(v.get_label().dna);

  return var_dna;
}


TEST_CASE("Get the reference sequence of a graph")
{
  using namespace gyper;
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

  graph = gyper::Graph(false /*use_absolute_positions*/);
  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());

  SECTION("Get the entire reference")
  {
    REQUIRE(graph.get_all_ref() == gyper::to_vec("SGTACGEEF"));
    REQUIRE(graph.get_ref(0, 999) == gyper::to_vec("SGTACGEEF"));
  }

  SECTION("Get the partial reference")
  {
    REQUIRE(graph.get_ref(0, 2) == gyper::to_vec("SG"));
    REQUIRE(graph.get_ref(1, 2) == gyper::to_vec("G"));
    REQUIRE(graph.get_ref(1, 9999) == gyper::to_vec("GTACGEEF"));
    REQUIRE(graph.get_ref(4, 7) == gyper::to_vec("CGE"));
    REQUIRE(graph.get_ref(8, 12) == gyper::to_vec("F"));
  }

  SECTION("Out of bounds reference")
  {
    REQUIRE(graph.get_ref(9, 12) == gyper::to_vec(""));
    REQUIRE(graph.get_ref(9999, 1999992) == gyper::to_vec(""));
  }

  SECTION("To is smaller")
  {
    REQUIRE(graph.get_ref(4, 1) == gyper::to_vec(""));
  }
}


TEST_CASE("Graph with a reference only.")
{
  using namespace gyper;
  std::vector<char> reference_sequence;
  char testdata[] = "ACCGGGAAAA";
  reference_sequence.insert(reference_sequence.end(), testdata, testdata + 10);

  graph = gyper::Graph(false /*use_absolute_positions*/);
  graph.add_genomic_region(std::move(reference_sequence), std::vector<gyper::VarRecord>(0), gyper::GenomicRegion());
  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  SECTION("The graph has the correct size")
  {
    REQUIRE(ref_nodes.size() == 1);
    REQUIRE(var_nodes.size() == 0);
  }

  SECTION("The only node has the correct order and dna")
  {
    REQUIRE(ref_nodes[0].out_degree() == 0);
    REQUIRE(ref_nodes[0].get_label().order == 0);
    REQUIRE(ref_nodes[0].get_label().dna.size() == 10);
    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec("ACCGGGAAAA"));
  }

  SECTION("The graph outputs the correct reference sequence")
  {
    REQUIRE(graph.get_all_ref() == gyper::to_vec("ACCGGGAAAA"));
  }
}


TEST_CASE("Graph with two variant records.")
{
  using namespace gyper;
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

  graph = gyper::Graph(false /*use_absolute_positions*/);
  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());
  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  SECTION("The graph should have the correct size")
  {
    REQUIRE(ref_nodes.size() == 3);
    REQUIRE(var_nodes.size() == 5);
    REQUIRE(graph.size() == 8);
  }

  SECTION("The nodes should be correctly connected")
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

  SECTION("The nodes should have the correct order")
  {
    REQUIRE(ref_nodes[0].get_label().order == 0);
    REQUIRE(var_nodes[0].get_label().order == 3);
    REQUIRE(var_nodes[1].get_label().order == 3);
    REQUIRE(ref_nodes[1].get_label().order == 4);
    REQUIRE(var_nodes[2].get_label().order == 6);
    REQUIRE(var_nodes[3].get_label().order == 6);
    REQUIRE(var_nodes[4].get_label().order == 6);
    REQUIRE(ref_nodes[2].get_label().order == 7);
  }

  SECTION("The nodes should have a label with the correct DNA bases")
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

  graph = gyper::Graph(false /*use_absolute_positions*/);
  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());
  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  SECTION("The graph should have the correct size")
  {
    REQUIRE(ref_nodes.size() == 3);
    REQUIRE(var_nodes.size() == 5);
    REQUIRE(graph.size() == 8);
  }

  SECTION("The nodes should be correctly connected")
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

  SECTION("The nodes should have the correct order")
  {
    REQUIRE(ref_nodes[0].get_label().order == 0);
    REQUIRE(var_nodes[0].get_label().order == 0);
    REQUIRE(var_nodes[1].get_label().order == 0);
    REQUIRE(ref_nodes[1].get_label().order == 1);
    REQUIRE(var_nodes[2].get_label().order == 6);
    REQUIRE(var_nodes[3].get_label().order == 6);
    REQUIRE(var_nodes[4].get_label().order == 6);
    REQUIRE(ref_nodes[2].get_label().order == 7);
  }

  SECTION("The nodes should have a label with the correct DNA bases")
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

  graph = gyper::Graph(false /*use_absolute_positions*/);
  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());

  REQUIRE(graph.ref_nodes.size() == 3);
  REQUIRE(graph.var_nodes.size() == 5);
  REQUIRE(graph.size() == 8);

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

  SECTION("The nodes should have the correct order")
  {
    REQUIRE(ref_nodes[0].get_label().order == 0);
    REQUIRE(var_nodes[0].get_label().order == 3);
    REQUIRE(var_nodes[1].get_label().order == 3);
    REQUIRE(ref_nodes[1].get_label().order == 4);
    REQUIRE(var_nodes[2].get_label().order == 6);
    REQUIRE(var_nodes[3].get_label().order == 6);
    REQUIRE(var_nodes[4].get_label().order == 6);
    REQUIRE(ref_nodes[2].get_label().order == 7);

  }

  SECTION("The nodes should have a label with the correct DNA bases")
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

  graph = gyper::Graph(false /*use_absolute_positions*/);
  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());

  REQUIRE(graph.ref_nodes.size() == 3);
  REQUIRE(graph.var_nodes.size() == 5);
  REQUIRE(graph.size() == 8);

  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  SECTION("The nodes should be correctly connected")
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

  SECTION("The nodes should have the correct order")
  {
    REQUIRE(ref_nodes[0].get_label().order == 0);
    REQUIRE(var_nodes[0].get_label().order == 3);
    REQUIRE(var_nodes[1].get_label().order == 3);
    REQUIRE(ref_nodes[1].get_label().order == 4);
    REQUIRE(var_nodes[2].get_label().order == 6);
    REQUIRE(var_nodes[3].get_label().order == 6);
    REQUIRE(var_nodes[4].get_label().order == 6);
    REQUIRE(ref_nodes[2].get_label().order == 7);
  }

  SECTION("The nodes should have a label with the correct DNA bases")
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

  graph = gyper::Graph(false /*use_absolute_positions*/);
  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion("chr1:2"));

  REQUIRE(graph.ref_nodes.size() == 3);
  REQUIRE(graph.var_nodes.size() == 5);
  REQUIRE(graph.size() == 8);

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

    REQUIRE(ref_nodes[1].out_degree() == 3);
    REQUIRE(ref_nodes[1].get_var_index(0) == 2);
    REQUIRE(ref_nodes[1].get_var_index(1) == 3);
    REQUIRE(ref_nodes[1].get_var_index(2) == 4);

    REQUIRE(var_nodes[2].get_out_ref_index() == 2);
    REQUIRE(var_nodes[3].get_out_ref_index() == 2);
    REQUIRE(var_nodes[4].get_out_ref_index() == 2);

    REQUIRE(ref_nodes[2].out_degree() == 0);
  }

  SECTION("The nodes should have the correct order")
  {
    REQUIRE(ref_nodes[0].get_label().order == 1);
    REQUIRE(var_nodes[0].get_label().order == 3);
    REQUIRE(var_nodes[1].get_label().order == 3);
    REQUIRE(ref_nodes[1].get_label().order == 5);
    REQUIRE(var_nodes[2].get_label().order == 6);
    REQUIRE(var_nodes[3].get_label().order == 6);
    REQUIRE(ref_nodes[2].get_label().order == 7);
  }

  SECTION("The nodes should have a label with the correct DNA bases")
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

  graph = gyper::Graph(false /*use_absolute_positions*/);
  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());
  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  SECTION("The graph should have the correct size")
  {
    REQUIRE(ref_nodes.size() == 2);
    REQUIRE(var_nodes.size() == 3);
    REQUIRE(graph.size() == 5);
  }

  SECTION("The nodes should be correctly connected")
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

  SECTION("The nodes should have the correct order")
  {
    REQUIRE(ref_nodes[0].get_label().order == 0);
    REQUIRE(var_nodes[0].get_label().order == 2);
    REQUIRE(var_nodes[1].get_label().order == 2);
    REQUIRE(var_nodes[2].get_label().order == 2);
    REQUIRE(ref_nodes[1].get_label().order == 5);
  }

  SECTION("The nodes should have a label with the correct DNA bases")
  {
    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec("AC"));
    REQUIRE(var_nodes[0].get_label().dna == gyper::to_vec("GGT"));
    REQUIRE(var_nodes[1].get_label().dna == gyper::to_vec("GATT"));
    REQUIRE(var_nodes[2].get_label().dna == gyper::to_vec("T"));
    REQUIRE(ref_nodes[1].get_label().dna == gyper::to_vec("AA"));
  }
}


TEST_CASE("Variants can overlap. Case where the second variant reaches further")
{
  using namespace gyper;
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

  graph = gyper::Graph(false /*use_absolute_positions*/);
  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());
  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  SECTION("The graph should have the correct size")
  {
    REQUIRE(ref_nodes.size() == 2);
    REQUIRE(var_nodes.size() == 3);
    REQUIRE(graph.size() == 5);
  }

  SECTION("The nodes should be correctly connected")
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

  SECTION("The nodes should have the correct order")
  {
    REQUIRE(ref_nodes[0].get_label().order == 0);
    REQUIRE(var_nodes[0].get_label().order == 2);
    REQUIRE(var_nodes[1].get_label().order == 2);
    REQUIRE(var_nodes[2].get_label().order == 2);
    REQUIRE(ref_nodes[1].get_label().order == 6);
  }

  SECTION("The nodes should have a label with the correct DNA bases")
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

  graph = gyper::Graph(false /*use_absolute_positions*/);
  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());
  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  SECTION("The graph should have the correct size")
  {
    REQUIRE(ref_nodes.size() == 3);
    REQUIRE(var_nodes.size() == 5);
    REQUIRE(graph.size() == 8);
  }

  SECTION("The nodes should be correctly connected")
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

  SECTION("The nodes should have the correct order")
  {
    REQUIRE(ref_nodes[0].get_label().order == 0);
    REQUIRE(var_nodes[0].get_label().order == 2);
    REQUIRE(var_nodes[1].get_label().order == 2);
    REQUIRE(ref_nodes[1].get_label().order == 3);
    REQUIRE(var_nodes[2].get_label().order == 3);
    REQUIRE(var_nodes[3].get_label().order == 3);
    REQUIRE(var_nodes[4].get_label().order == 3);
    REQUIRE(ref_nodes[2].get_label().order == 4);

  }

  SECTION("The nodes should have a label with the correct DNA bases")
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


TEST_CASE("When merging a deletion covering multiple short variants, all combinations of the variants need to be added.", "[test2]")
{
  using namespace gyper;
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

  graph = gyper::Graph(false /*use_absolute_positions*/);
  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());
  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  SECTION("The graph should have the correct size")
  {
    REQUIRE(ref_nodes.size() == 2);
    REQUIRE(var_nodes.size() == 10);
    REQUIRE(graph.size() == 12);
  }

  SECTION("The nodes should be correctly connected")
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

  SECTION("The nodes should have the correct order")
  {
    REQUIRE(ref_nodes[0].get_label().order == 0);

    REQUIRE(var_nodes[0].get_label().order == 2);
    REQUIRE(var_nodes[1].get_label().order == 2);
    REQUIRE(var_nodes[2].get_label().order == 2);
    REQUIRE(var_nodes[3].get_label().order == 2);
    REQUIRE(var_nodes[4].get_label().order == 2);
    REQUIRE(var_nodes[5].get_label().order == 2);
    REQUIRE(var_nodes[6].get_label().order == 2);
    REQUIRE(var_nodes[7].get_label().order == 2);
    REQUIRE(var_nodes[8].get_label().order == 2);
    REQUIRE(var_nodes[9].get_label().order == 2);

    REQUIRE(ref_nodes[1].get_label().order == 7);
  }

  SECTION("The nodes should have a label with the correct DNA bases")
  {
    std::vector<std::vector<char> > var_dna = get_var_dna(var_nodes);

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

  graph = gyper::Graph(false /*use_absolute_positions*/);
  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());
  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  SECTION("The nodes should have a label with the correct DNA bases")
  {
    std::vector<std::vector<char> > var_dna = get_var_dna(var_nodes);

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

  SECTION("The graph should have the correct size")
  {
    REQUIRE(ref_nodes.size() == 2);
    REQUIRE(var_nodes.size() == 10);
    REQUIRE(graph.size() == 12);
  }

  SECTION("The nodes should be correctly connected")
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

  SECTION("The nodes should have the correct order")
  {
    REQUIRE(ref_nodes[0].get_label().order == 0);

    REQUIRE(var_nodes[0].get_label().order == 0);
    REQUIRE(var_nodes[1].get_label().order == 0);
    REQUIRE(var_nodes[2].get_label().order == 0);
    REQUIRE(var_nodes[3].get_label().order == 0);
    REQUIRE(var_nodes[4].get_label().order == 0);
    REQUIRE(var_nodes[5].get_label().order == 0);
    REQUIRE(var_nodes[6].get_label().order == 0);
    REQUIRE(var_nodes[7].get_label().order == 0);
    REQUIRE(var_nodes[8].get_label().order == 0);
    REQUIRE(var_nodes[9].get_label().order == 0);

    REQUIRE(ref_nodes[1].get_label().order == 4);
  }

  Options::instance()->add_all_variants = false;
}


TEST_CASE("Four variants joined", "[test2]")
{
  using namespace gyper;
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

  graph = gyper::Graph(false /*use_absolute_positions*/);
  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());
  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  SECTION("The nodes should have a label with the correct DNA bases")
  {
    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec("S"));
    std::vector<std::vector<char> > var_dna = get_var_dna(var_nodes);
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

  SECTION("The graph should have the correct size")
  {
    REQUIRE(ref_nodes.size() == 2);
    REQUIRE(var_nodes.size() == 13);
    REQUIRE(graph.size() == 15);
  }

  SECTION("The nodes should be correctly connected")
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

  SECTION("The nodes should have the correct order")
  {
    REQUIRE(ref_nodes[0].get_label().order == 0);

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
    REQUIRE(var_nodes[10].get_label().order == 1);
    REQUIRE(var_nodes[11].get_label().order == 1);
    REQUIRE(var_nodes[12].get_label().order == 1);

    REQUIRE(ref_nodes[1].get_label().order == 6);
  }

  Options::instance()->add_all_variants = false;
}


TEST_CASE("Variants of any number can be joined, here 3 are tested.", "[test]")
{
  using namespace gyper;
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

  graph = gyper::Graph(false /*use_absolute_positions*/);
  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());
  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  SECTION("The nodes should have a label with the correct DNA bases")
  {
    std::vector<std::vector<char> > var_dna = get_var_dna(var_nodes);
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

    REQUIRE(ref_nodes[1].get_label().dna ==  gyper::to_vec("F"));
  }

  SECTION("The graph should have the correct size")
  {
    REQUIRE(ref_nodes.size() == 2);
    REQUIRE(var_nodes.size() == 25);
    REQUIRE(graph.size() == 27);
  }

  SECTION("The nodes should be correctly connected")
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

TEST_CASE("Check if the above creates the graph correctly")
{
  using namespace gyper;
  Options::instance()->add_all_variants = true;
  std::vector<char> reference_sequence;
  char testdata[] = "CG1xGAaTGTTGGGGGGGAACAGT";
  reference_sequence.insert(reference_sequence.end(), testdata, testdata + 24);

  std::vector<gyper::VarRecord> records_step1;
  std::vector<gyper::VarRecord> records_step2;
  std::vector<gyper::VarRecord> records_step3;

  {
    gyper::VarRecord record;
    record.pos = 1;
    record.ref = gyper::to_vec("G1");
    record.alts = {gyper::to_vec("A2"), gyper::to_vec("G3")};
    records_step1.push_back(record);
    records_step2.push_back(record);
    records_step3.push_back(record);

    record.pos = 2;
    record.ref = gyper::to_vec("1xGA");
    record.alts = {gyper::to_vec("3yAG")};
    records_step1.push_back(record);
    records_step2.push_back(record);
    records_step3.push_back(record);

    record.pos = 4;
    record.ref = gyper::to_vec("GA");
    record.alts = {gyper::to_vec("AG")};
    records_step2.push_back(record);
    records_step3.push_back(record);

    record.pos = 5;
    record.ref = gyper::to_vec("Aa");
    record.alts = {gyper::to_vec("Gb"), gyper::to_vec("Gc")};
    records_step3.push_back(record);
  }

  SECTION("Step 1: Simpler graph")
  {
    graph = gyper::Graph(false /*use_absolute_positions*/);
    graph.add_genomic_region(std::move(reference_sequence), std::move(records_step1), gyper::GenomicRegion());
    std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
    std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec("C"));

    REQUIRE(var_nodes[0].get_label().dna == gyper::to_vec("G1xGA"));
    REQUIRE(var_nodes[1].get_label().dna == gyper::to_vec("A2xGA"));
    REQUIRE(var_nodes[2].get_label().dna == gyper::to_vec("G3xGA"));
    REQUIRE(var_nodes[3].get_label().dna == gyper::to_vec("G3yAG"));

    REQUIRE(var_nodes.size() == 4);
  }

  SECTION("Step 2: Simpler graph")
  {
    graph = gyper::Graph(false /*use_absolute_positions*/);
    graph.add_genomic_region(std::move(reference_sequence), std::move(records_step2), gyper::GenomicRegion());
    std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
    std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

    std::vector<std::vector<char> > var_dna;

    for (auto v : var_nodes)
    {
      var_dna.push_back(v.get_label().dna);
    }

    REQUIRE(ref_nodes[0].get_label().dna ==   gyper::to_vec("C"));

    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("G1xGA")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("G1xAG")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("G3yAG")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("A2xGA")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("A2xAG")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("G3xGA")) != var_dna.end());
    REQUIRE(std::find(var_dna.cbegin(), var_dna.cend(), gyper::to_vec("G3xAG")) != var_dna.end());

    REQUIRE(var_nodes.size() == 7);
  }

  Options::instance()->add_all_variants = false;
}


TEST_CASE("Variant overlapping a N on the reference genome")
{
  using namespace gyper;
  Options::instance()->add_all_variants = true;
  std::vector<char> reference_sequence;
  char testdata[] = "GCTGCGGCGGGCGTCGCGGCCGCCCCCGGGGAGCCCGGCGGGCGCCGGCGCGNCCCCCCCCCCACCCCACGTCTCGTCGCGCGCGC";
  reference_sequence.insert(reference_sequence.end(), testdata, testdata + 86);

  SECTION("The reference allele has a N")
  {
    std::vector<gyper::VarRecord> records;

    {
      gyper::VarRecord record;
      record.pos = 52;
      record.ref = gyper::to_vec("GN");
      record.alts = {gyper::to_vec("GA")};
      records.push_back(record);
    }

    graph = gyper::Graph(false /*use_absolute_positions*/);
    graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());

    // Here, we assume that nothing is added
    std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
    std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

    REQUIRE(ref_nodes.size() == 1);
    REQUIRE(ref_nodes[0].get_label().dna == reference_sequence);
    REQUIRE(var_nodes.size() == 0);
  }

  SECTION("If an alternative allele has a N, we remove it")
  {
    std::vector<gyper::VarRecord> records;

    {
      gyper::VarRecord record;
      record.pos = 52;
      record.ref = gyper::to_vec("G");
      record.alts = {gyper::to_vec("GN"), gyper::to_vec("GA")};
      records.push_back(record);
    }

    graph = gyper::Graph(false /*use_absolute_positions*/);
    graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());

    std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
    std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

    REQUIRE(ref_nodes.size() == 2);
    REQUIRE(var_nodes.size() == 2);

    REQUIRE(var_nodes[0].get_label().dna == gyper::to_vec("G"));
    REQUIRE(var_nodes[1].get_label().dna == gyper::to_vec("GA"));
  }

  SECTION("If all alternative allele have a N, we remove the variant")
  {
    std::vector<gyper::VarRecord> records;

    {
      gyper::VarRecord record;
      record.pos = 52;
      record.ref = gyper::to_vec("G");
      record.alts = {gyper::to_vec("GN"), gyper::to_vec("GNN")};
      records.push_back(record);
    }

    graph = gyper::Graph(false /*use_absolute_positions*/);
    graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());

    // Here, we assume that nothing is added
    std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
    std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

    REQUIRE(ref_nodes.size() == 1);
    REQUIRE(ref_nodes[0].get_label().dna == reference_sequence);
    REQUIRE(var_nodes.size() == 0);
  }
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

  graph = gyper::Graph(false /*use_absolute_positions*/);
  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());

  // Here, we assume that nothing is added
  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  SECTION("The number of nodes is correct")
  {
    REQUIRE(ref_nodes.size() == 3);
    REQUIRE(var_nodes.size() == 4);
  }

  SECTION("We have the correct labels on the reference nodes")
  {
    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec("GT"));
    REQUIRE(ref_nodes[1].get_label().dna == gyper::to_vec(""));
    REQUIRE(ref_nodes[2].get_label().dna == gyper::to_vec("ATG"));
  }

  SECTION("We have the correct labels on the variant nodes")
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
  char testdata[] = "STAAAATF";
  reference_sequence.insert(reference_sequence.end(), testdata, testdata + 8);

  std::vector<gyper::VarRecord> records;

  {
    gyper::VarRecord record;
    record.pos = 1;
    record.ref = gyper::to_vec("TAAAAT"); // TA
    record.alts = {gyper::to_vec("TAAAT")}; // T
    records.push_back(record);

    record.pos = 4;
    record.ref = gyper::to_vec("A");
    record.alts = {gyper::to_vec("T")};
    records.push_back(record);
  }

  graph = gyper::Graph(false /*use_absolute_positions*/);
  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());

  // Here, we assume that nothing is added
  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  REQUIRE(ref_nodes.size() == 2);
  REQUIRE(var_nodes.size() == 4);

  SECTION("We have the correct labels on the reference nodes")
  {
    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec("S"));
    REQUIRE(ref_nodes[1].get_label().dna == gyper::to_vec("ATF"));
  }

  SECTION("We have the correct labels on the variant nodes")
  {
    REQUIRE(var_nodes[0].get_label().dna == gyper::to_vec("TAAA"));
    REQUIRE(var_nodes[1].get_label().dna == gyper::to_vec("TAA"));
    REQUIRE(var_nodes[2].get_label().dna == gyper::to_vec("TAAT"));
    REQUIRE(var_nodes[3].get_label().dna == gyper::to_vec("TAT"));
  }

  Options::instance()->add_all_variants = true;
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

  graph = gyper::Graph(false /*use_absolute_positions*/);
  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());

  // Here, we assume that nothing is added
  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  SECTION("The number of nodes is corrent")
  {
    REQUIRE(ref_nodes.size() == 2);
    REQUIRE(var_nodes.size() == 3);
  }

  SECTION("The ref nodes have the correct DNA bases")
  {
    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec("ST"));
    REQUIRE(ref_nodes[1].get_label().dna == gyper::to_vec("F"));
  }

  SECTION("The var nodes have the correct DNA bases")
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

  graph = gyper::Graph(false /*use_absolute_positions*/);
  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());

  // Here, we assume that nothing is added
  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  SECTION("The number of nodes is corrent")
  {
    REQUIRE(ref_nodes.size() == 4);
    REQUIRE(var_nodes.size() == 6);
  }

  SECTION("The ref nodes have the correct DNA bases")
  {
    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec("ST"));
    REQUIRE(ref_nodes[1].get_label().dna == gyper::to_vec(""));
    REQUIRE(ref_nodes[2].get_label().dna == gyper::to_vec(""));
    REQUIRE(ref_nodes[3].get_label().dna == gyper::to_vec("F"));
  }

  SECTION("The var nodes have the correct DNA bases")
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
  Options::instance()->add_all_variants = false;

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

  graph = gyper::Graph(false /*use_absolute_positions*/);
  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());

  // Here, we assume that nothing is added
  std::vector<gyper::RefNode> const & ref_nodes = graph.ref_nodes;
  std::vector<gyper::VarNode> const & var_nodes = graph.var_nodes;

  SECTION("The number of nodes is corrent")
  {
    REQUIRE(ref_nodes.size() == 2);
    REQUIRE(var_nodes.size() == 3);
  }

  SECTION("The ref nodes have the correct DNA bases")
  {
    REQUIRE(ref_nodes[0].get_label().dna == gyper::to_vec("TGCAAATCTCATATATATATATATATATAT"));
    REQUIRE(ref_nodes[1].get_label().dna == gyper::to_vec("TTTTTTTTTTTTTA"));
  }

  SECTION("The var nodes have the correct DNA bases")
  {
    REQUIRE(var_nodes[0].get_label().dna == gyper::to_vec("ATATATATATATATATTTTTTTTTTTT"));
    REQUIRE(var_nodes[1].get_label().dna == gyper::to_vec("A"));
    REQUIRE(var_nodes[2].get_label().dna == gyper::to_vec("ATATATATAT"));
  }
}
