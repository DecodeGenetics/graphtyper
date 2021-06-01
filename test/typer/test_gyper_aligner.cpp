#include <climits>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>

#include <graphtyper/graph/graph.hpp>
#include <graphtyper/index/rocksdb.hpp>
#include <graphtyper/typer/aligner.hpp>
#include <graphtyper/utilities/type_conversions.hpp> // to_uint64()

#include <catch.hpp>

typedef gyper::Aligner<gyper::SimpleGraph, gyper::Index<gyper::RocksDB>> TAligner;

TEST_CASE("Align read from start", "[align]")
{
  std::stringstream my_graph;
  my_graph << gyper_SOURCE_DIRECTORY << "/test/data/graphs/index_test_chr1.grf";
  std::stringstream my_index;
  my_index << gyper_SOURCE_DIRECTORY << "/test/data/graphs/index_test_chr1";
  TAligner aligner(my_graph.str().c_str(), my_index.str().c_str());
  seqan::String<seqan::Dna> query = "TTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTTTGCCCAGGTTTCCCCAGGTTTCCCCTTTGGA";
  aligner.align_read_from_start(query);
}

TEST_CASE("GENERAL PURPOSE", "[align]")
{
  std::stringstream my_graph;
  my_graph << gyper_SOURCE_DIRECTORY << "/test/data/graphs/index_test_chr1.grf";
  std::stringstream my_index;
  my_index << gyper_SOURCE_DIRECTORY << "/test/data/graphs/index_test_chr1";

  TAligner aligner(my_graph.str().c_str(), my_index.str().c_str());
  REQUIRE(aligner.index.size() > 0);
  REQUIRE(aligner.writer.graph.size() > 0);

  SECTION("Common k-mer on the reference")
  {
    seqan::String<seqan::Dna> query = "TTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTT";
    uint64_t tst = to_uint64(query);
    REQUIRE(aligner.index.get(to_uint64(query)).size() == 3);
    std::vector<gyper::KmerLabel> matching_kmers = aligner.index.get(to_uint64(query));

    for (auto match_it = matching_kmers.begin(); match_it != matching_kmers.end(); ++match_it)
    {
      if (match_it->start_index == 23 || match_it->start_index == 13 || match_it->start_index == 3)
        REQUIRE(true);
      else
        REQUIRE(false);

      if (match_it->end_index == 54 || match_it->end_index == 44 || match_it->end_index == 34)
        REQUIRE(true);
      else
        REQUIRE(false);
    }
  }

  SECTION("Unique k-mer on the reference")
  {
    seqan::String<seqan::Dna> query = "TTCCCCAGGTTTCCCCAGGTTTCCCCTTTGGA";
    REQUIRE(aligner.index.get(to_uint64(query)).size() == 1);
    std::vector<gyper::KmerLabel> matching_kmers = aligner.index.get(to_uint64(query));

    auto match_it = matching_kmers.begin();
    REQUIRE(match_it->start_index == 34);
    REQUIRE(match_it->end_index == 65);
    // REQUIRE(match_it->variant_num == 1);
    ++match_it;
    REQUIRE(match_it == matching_kmers.end());
  }

  SECTION("Unique k-mer on the variant")
  {
    seqan::String<seqan::Dna> query = "TTGCCCAGGTTTCCCCAGGTTTCCCCTTTGGA";
    REQUIRE(aligner.index.get(to_uint64(query)).size() == 1);
    std::vector<gyper::KmerLabel> matching_kmers = aligner.index.get(to_uint64(query));

    auto match_it = matching_kmers.begin();
    REQUIRE(match_it->start_index == 34);
    REQUIRE(match_it->end_index == 65);
    // REQUIRE(match_it->variant_list.front().bits.to_ulong() == 2);
    ++match_it;
    REQUIRE(match_it == matching_kmers.end());
  }

  SECTION("Non-existing k-mer")
  {
    seqan::String<seqan::Dna> query = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    REQUIRE(aligner.index.get(to_uint64(query)).size() == 0);
    std::vector<gyper::KmerLabel> matching_kmers = aligner.index.get(to_uint64(query));

    auto match_it = matching_kmers.begin();
    REQUIRE(match_it == matching_kmers.end());
  }
}
