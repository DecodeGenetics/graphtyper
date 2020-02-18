#include <catch.hpp>

#include <stdio.h>
#include <climits>
#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>

#include <graphtyper/graph/graph_serialization.hpp>
#include <graphtyper/graph/constructor.hpp>
#include <graphtyper/index/indexer.hpp>
#include <graphtyper/index/rocksdb.hpp>
#include <graphtyper/index/mem_index.hpp>
#include <graphtyper/typer/graph_swapper.hpp>


TEST_CASE("Test graph and index swapping with RocksDB index")
{
  using namespace gyper;
  // graph 1
  std::stringstream my_graph1;
  my_graph1 << gyper_SOURCE_DIRECTORY << "/test/data/graphs/index_test_chr1.grf";
  std::stringstream my_index1;
  my_index1 << gyper_SOURCE_DIRECTORY << "/test/data/graphs/index_test_chr1";

  // graph 2
  std::stringstream my_graph2;
  my_graph2 << gyper_SOURCE_DIRECTORY << "/test/data/graphs/index_test_chr2.grf";
  std::stringstream my_index2;
  my_index2 << gyper_SOURCE_DIRECTORY << "/test/data/graphs/index_test_chr2";

  gyper::load_graph(my_graph1.str());
  gyper::load_index(my_index1.str());

  REQUIRE(graph.size() > 0);
  REQUIRE(gyper::index.check());

  Graph graph2 = gyper::load_secondary_graph(my_graph2.str());
  gyper::Index<gyper::RocksDB> index2 = gyper::load_secondary_index(my_index2.str());

  REQUIRE(graph2.size() > 0);
  REQUIRE(index2.check());

  REQUIRE(graph.size() != graph2.size());
  std::size_t const N1 = graph.size();
  std::size_t const N2 = graph2.size();
  std::swap(graph, graph2);
  std::swap(gyper::index, index2);
  REQUIRE(N1 == graph2.size());
  REQUIRE(N2 == graph.size());
  REQUIRE(gyper::index.check());
  REQUIRE(index2.check());
  std::swap(graph, graph2);
  std::swap(gyper::index, index2);
  REQUIRE(N1 == graph.size());
  REQUIRE(N2 == graph2.size());
}


TEST_CASE("Test graph and index swapping with SparseHash index")
{
  using namespace gyper;
  // graph 1
  std::stringstream my_graph1;
  my_graph1 << gyper_SOURCE_DIRECTORY << "/test/data/graphs/index_test_chr1.grf";
  std::stringstream my_index1;
  my_index1 << gyper_SOURCE_DIRECTORY << "/test/data/graphs/index_test_chr1";

  // graph 2
  std::stringstream my_graph2;
  my_graph2 << gyper_SOURCE_DIRECTORY << "/test/data/graphs/index_test_chr2.grf";
  std::stringstream my_index2;
  my_index2 << gyper_SOURCE_DIRECTORY << "/test/data/graphs/index_test_chr2";

  gyper::load_graph(my_graph1.str());
  gyper::load_index(my_index1.str());
  mem_index.load(gyper::index);
  REQUIRE(mem_index.hamming0.size() > 0);

  REQUIRE(graph.size() > 0);
  REQUIRE(gyper::index.check());

  Graph graph2 = gyper::load_secondary_graph(my_graph2.str());
  MemIndex mem_index2 = load_secondary_mem_index(my_index2.str(), graph2);
  REQUIRE(graph2.size() > 0);
  REQUIRE(mem_index2.hamming0.size() > 0);

  REQUIRE(graph.size() != graph2.size());
  REQUIRE(mem_index.hamming0.size() != mem_index2.hamming0.size());

  std::size_t const N1 = graph.size();
  std::size_t const N2 = graph2.size();
  std::size_t const M1 = mem_index.hamming0.size();
  std::size_t const M2 = mem_index2.hamming0.size();

  swap_graph_and_index(graph2, mem_index2);

  REQUIRE(N1 == graph2.size());
  REQUIRE(N2 == graph.size());
  REQUIRE(M1 == mem_index2.hamming0.size());
  REQUIRE(M2 == mem_index.hamming0.size());

  swap_graph_and_index(graph2, mem_index2);

  REQUIRE(N1 == graph.size());
  REQUIRE(N2 == graph2.size());
  REQUIRE(M1 == mem_index.hamming0.size());
  REQUIRE(M2 == mem_index2.hamming0.size());
}
