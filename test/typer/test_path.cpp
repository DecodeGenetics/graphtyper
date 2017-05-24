#include <catch.hpp>

#include <stdio.h>
#include <climits>
#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>

#include <graphtyper/graph/graph.hpp>
#include <graphtyper/graph/graph_serialization.hpp>
#include <graphtyper/index/rocksdb.hpp>
#include <graphtyper/index/indexer.hpp>
#include <graphtyper/utilities/type_conversions.hpp> // to_uint64()
#include <graphtyper/typer/path.hpp>
#include <graphtyper/index/kmer_label.hpp>


TEST_CASE("Path with one label")
{
  SECTION("Label is of reference")
  {
    gyper::KmerLabel new_label(1, 32);
    gyper::Path new_path(new_label, 0, 31);

    REQUIRE(new_path.size() == 32);
    REQUIRE(new_path.start == 1);
    REQUIRE(new_path.end == 32);
    REQUIRE(new_path.var_order.size() == 0);
    REQUIRE(new_path.nums.size() == 0);
  }

  SECTION("Label is of some variant")
  {
    gyper::KmerLabel new_label(1, 32, 0);
    new_label.variant_order = 20;
    new_label.variant_num = 2;
    gyper::Path new_path(new_label, 0, 31);

    REQUIRE(new_path.size() == 32);
    REQUIRE(new_path.start == 1);
    REQUIRE(new_path.end == 32);
    REQUIRE(new_path.var_order.size() == 1);
    REQUIRE(new_path.nums.size() == 1);
    REQUIRE(new_path.var_order[0] == 20);
    REQUIRE(new_path.nums[0] == std::bitset<gyper::MAX_NUMBER_OF_HAPLOTYPES>(1 << 2));
  }
}

TEST_CASE("Sucessful merging of two paths")
{
  SECTION("Two reference paths can merge")
  {
    gyper::Path path_prev(gyper::KmerLabel(1, 32), 0, 31);
    gyper::Path path_next(gyper::KmerLabel(32, 43), 31, 62);

    gyper::Path merged_path(path_prev, path_next);

    REQUIRE(merged_path.size() == 63);
    REQUIRE(merged_path.start == 1);
    REQUIRE(merged_path.end == 43);
    REQUIRE(merged_path.var_order.size() == 0);
    REQUIRE(merged_path.nums.size() == 0);
  }

  SECTION("Two paths with the same variant id and number can successfully merge")
  {
    gyper::KmerLabel prev_label(1, 32, 0);
    gyper::KmerLabel next_label(32, 43, 0);

    prev_label.variant_num = 2;
    prev_label.variant_order = 20;
    next_label.variant_num = 2;
    next_label.variant_order = 20;

    gyper::Path path_prev(prev_label, 0, 31);
    gyper::Path path_next(next_label, 31, 62);

    gyper::Path merged_path(path_prev, path_next);

    REQUIRE(merged_path.size() == 63);
    REQUIRE(merged_path.start == 1);
    REQUIRE(merged_path.end == 43);
    REQUIRE(merged_path.var_order.size() == 1);
    REQUIRE(merged_path.nums.size() == 1);
    REQUIRE(merged_path.var_order[0] == 20);
    REQUIRE(merged_path.nums[0] == std::bitset<gyper::MAX_NUMBER_OF_HAPLOTYPES>(1 << 2));
  }

  SECTION("Two paths with the different variant id can merge")
  {
    gyper::KmerLabel label_prev(1, 32, 0);
    gyper::KmerLabel label_next(32, 43, 0);

    label_prev.variant_num = 2;
    label_prev.variant_order = 20;
    label_next.variant_num = 3;
    label_next.variant_order = 34;

    gyper::Path path_prev(label_prev, 0, 31);
    gyper::Path path_next(label_next, 31, 62);

    gyper::Path merged_path(path_prev, path_next);

    REQUIRE(merged_path.size() == 63);
    REQUIRE(merged_path.start == 1);
    REQUIRE(merged_path.end == 43);
    REQUIRE(merged_path.var_order.size() == 2);
    REQUIRE(merged_path.nums.size() == 2);
    REQUIRE(merged_path.var_order[0] == 34);
    REQUIRE(merged_path.nums[0] == std::bitset<gyper::MAX_NUMBER_OF_HAPLOTYPES>(1 << 3));
    REQUIRE(merged_path.var_order[1] == 20);
    REQUIRE(merged_path.nums[1] == std::bitset<gyper::MAX_NUMBER_OF_HAPLOTYPES>(1 << 2));
  }
}


TEST_CASE("Failed merging of two paths")
{
  SECTION("The two paths cannot merge because of different variant numbers")
  {
    gyper::KmerLabel label_prev(1, 32, 0);
    gyper::KmerLabel label_next(32, 43, 0);

    label_prev.variant_num = 2;
    label_prev.variant_order = 20;
    label_next.variant_num = 99;
    label_next.variant_order = 20;

    gyper::Path path_prev(label_prev, 0, 31);
    gyper::Path path_next(label_next, 31, 62);

    gyper::Path merged_path(path_prev, path_next);
  }
}


TEST_CASE("Test index chr3")
{
  using namespace gyper;
  // AAAACAAAATAAAACAAAATAAAAGAAAACAAAATAAAACAAAATAAAAGAAAACATTATAAAACA
  // chr3 31 rs4 A G,GA

  std::cout << "Test index chr3" << std::endl;
  std::stringstream my_graph;
  my_graph << gyper_SOURCE_DIRECTORY << "/test/data/graphs/index_test_chr3.grf";
  std::stringstream my_index;
  my_index << gyper_SOURCE_DIRECTORY << "/test/data/graphs/index_test_chr3";

  gyper::load_graph(my_graph.str().c_str());
  REQUIRE(graph.size() > 0);
  gyper::load_index(my_index.str());

  SECTION("Kmer that starts at zero is has matches positions with the next one (although this is not a real alignmment)")
  {
    // kmer 1
    std::vector<gyper::KmerLabel> labels1 = gyper::index.get(gyper::to_uint64("AAAACAAAATAAAACAAAATAAAAGAAAACGA", 0));
    REQUIRE(labels1.size() == 2);
    REQUIRE(labels1[0].start_index == 0);
    REQUIRE(labels1[1].start_index == 0);

    REQUIRE(labels1[0].end_index == 31);
    REQUIRE(graph.get_actual_pos(labels1[0].end_index) == 31);
    REQUIRE(graph.get_ref_reach_pos(labels1[0].end_index) == 31);
    REQUIRE(labels1[1].end_index == 31);

    REQUIRE(labels1[0].variant_id == 1);
    REQUIRE(labels1[1].variant_id == 2);

    // kmer 2
    std::vector<gyper::KmerLabel> labels2 = gyper::index.get(gyper::to_uint64("AAATAAAACAAAATAAAAGAAAACATTATAAA", 0));
    REQUIRE(labels2.size() == 3);
    REQUIRE(labels2[0].start_index == 31);
    REQUIRE(labels2[0].end_index == 62);
    REQUIRE(labels2[0].variant_id == 0);
    REQUIRE(labels2[1].start_index == 31);
    REQUIRE(labels2[1].end_index == 62);
    REQUIRE(labels2[1].variant_id == 1);
    REQUIRE(labels2[2].start_index == 32);
    REQUIRE(labels2[2].end_index == 62);
    REQUIRE(labels2[2].variant_id == 2);
  }
}
