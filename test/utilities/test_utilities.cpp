#include <climits>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/graph_serialization.hpp>
#include <graphtyper/utilities/kmer_help_functions.hpp>
#include <graphtyper/utilities/type_conversions.hpp>

#include <catch2/catch.hpp>

TEST_CASE("Converting reads", "[utils]")
{
  using namespace gyper;

  {
    std::string read = "TTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTTTGCCCAGGTTTCCCCAGGTTTCCCCTTTGGA";
    std::string kmer1 = "TTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTT";
    std::string kmer2 = "TTGCCCAGGTTTCCCCAGGTTTCCCCTTTGGA";
    REQUIRE(to_uint64(read, 0) == to_uint64(kmer1, 0));
    REQUIRE(to_dna_str(to_uint64(read, 0)) == kmer1);
    REQUIRE(to_uint64(read, 31) == to_uint64(kmer2, 0));
    REQUIRE(to_dna_str(to_uint64(read, 31)) == kmer2);
  }
}

TEST_CASE("Mismatches of the last base")
{
  using namespace gyper;

  SECTION("Last base is 'A'")
  {
    std::string kmer = "ATTCCCCAGGTTTCCCCAGGTTTCCCCAGGTA";
    std::array<uint64_t, 3> mismatches = get_mismatches_of_last_base(to_uint64(kmer, 0));

    REQUIRE(to_dna_str(mismatches[0]) == "ATTCCCCAGGTTTCCCCAGGTTTCCCCAGGTC");
    REQUIRE(to_dna_str(mismatches[1]) == "ATTCCCCAGGTTTCCCCAGGTTTCCCCAGGTG");
    REQUIRE(to_dna_str(mismatches[2]) == "ATTCCCCAGGTTTCCCCAGGTTTCCCCAGGTT");
  }

  SECTION("Last base is 'C'")
  {
    std::string kmer = "TTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTC";
    std::array<uint64_t, 3> mismatches = get_mismatches_of_last_base(to_uint64(kmer, 0));

    REQUIRE(to_dna_str(mismatches[0]) == "TTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTA");
    REQUIRE(to_dna_str(mismatches[1]) == "TTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTG");
    REQUIRE(to_dna_str(mismatches[2]) == "TTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTT");
  }

  SECTION("Last base is 'G'")
  {
    std::string kmer = "CTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTG";
    std::array<uint64_t, 3> mismatches = get_mismatches_of_last_base(to_uint64(kmer, 0));

    REQUIRE(to_dna_str(mismatches[0]) == "CTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTA");
    REQUIRE(to_dna_str(mismatches[1]) == "CTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTC");
    REQUIRE(to_dna_str(mismatches[2]) == "CTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTT");
  }

  SECTION("Last base is 'T'")
  {
    std::string kmer = "GATCCCCAGGTTTCCCCAGGTTTCCCCAGGTT";
    std::array<uint64_t, 3> mismatches = get_mismatches_of_last_base(to_uint64(kmer, 0));

    REQUIRE(to_dna_str(mismatches[0]) == "GATCCCCAGGTTTCCCCAGGTTTCCCCAGGTA");
    REQUIRE(to_dna_str(mismatches[1]) == "GATCCCCAGGTTTCCCCAGGTTTCCCCAGGTC");
    REQUIRE(to_dna_str(mismatches[2]) == "GATCCCCAGGTTTCCCCAGGTTTCCCCAGGTG");
  }
}

TEST_CASE("Mismatches of the first base")
{
  using namespace gyper;

  // SECTION("First base is 'A'")
  {
    std::string kmer = "ATTCCCCAGGTTTCCCCAGGTTTCCCCAGGTA";
    std::array<uint64_t, 3> mismatches = get_mismatches_of_first_base(to_uint64(kmer, 0));

    REQUIRE(to_dna_str(mismatches[0]) == "CTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTA");
    REQUIRE(to_dna_str(mismatches[1]) == "GTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTA");
    REQUIRE(to_dna_str(mismatches[2]) == "TTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTA");
  }

  // SECTION("First base is 'C'")
  {
    std::string kmer = "CTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTC";
    std::array<uint64_t, 3> mismatches = get_mismatches_of_first_base(to_uint64(kmer, 0));

    REQUIRE(to_dna_str(mismatches[0]) == "ATTCCCCAGGTTTCCCCAGGTTTCCCCAGGTC");
    REQUIRE(to_dna_str(mismatches[1]) == "GTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTC");
    REQUIRE(to_dna_str(mismatches[2]) == "TTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTC");
  }

  // SECTION("First base is 'G'")
  {
    std::string kmer = "GTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTG";
    std::array<uint64_t, 3> mismatches = get_mismatches_of_first_base(to_uint64(kmer, 0));

    REQUIRE(to_dna_str(mismatches[0]) == "ATTCCCCAGGTTTCCCCAGGTTTCCCCAGGTG");
    REQUIRE(to_dna_str(mismatches[1]) == "CTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTG");
    REQUIRE(to_dna_str(mismatches[2]) == "TTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTG");
  }

  // SECTION("First base is 'T'")
  {
    std::string kmer = "TTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTA";
    std::array<uint64_t, 3> mismatches = get_mismatches_of_first_base(to_uint64(kmer, 0));

    REQUIRE(to_dna_str(mismatches[0]) == "ATTCCCCAGGTTTCCCCAGGTTTCCCCAGGTA");
    REQUIRE(to_dna_str(mismatches[1]) == "CTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTA");
    REQUIRE(to_dna_str(mismatches[2]) == "GTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTA");
  }
}

TEST_CASE("Hamming distance 1")
{
  using namespace gyper;

  SECTION("A homopolymer")
  {
    std::string kmer = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    auto hamming1_array = to_uint64_vec_hamming_distance_1(to_uint64(kmer, 0));

    // We should find 32-mers in hamming distance 1
    REQUIRE(std::find(hamming1_array.begin(), hamming1_array.end(), to_uint64("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 0)) ==
            hamming1_array.end());
    REQUIRE(std::find(hamming1_array.begin(), hamming1_array.end(), to_uint64("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAC", 0)) !=
            hamming1_array.end());
    REQUIRE(std::find(hamming1_array.begin(), hamming1_array.end(), to_uint64("AAAAAAAAAAAAAAAAAAAAAAAAAAAAGAAA", 0)) !=
            hamming1_array.end());
    REQUIRE(std::find(hamming1_array.begin(), hamming1_array.end(), to_uint64("AAAAAAAAAAAAAAAAAAAAAAAAATAAAAAA", 0)) !=
            hamming1_array.end());
    REQUIRE(std::find(hamming1_array.begin(), hamming1_array.end(), to_uint64("AAAAAAAAAAAAAAAAATAAAAAAAAAAAAAA", 0)) !=
            hamming1_array.end());
    REQUIRE(std::find(hamming1_array.begin(), hamming1_array.end(), to_uint64("AAAAAAAAAAAAAAAACAAAAAAAAAAAAAAA", 0)) !=
            hamming1_array.end());
    REQUIRE(std::find(hamming1_array.begin(), hamming1_array.end(), to_uint64("AAAAAAAAAAAAAAAAGAAAAAAAAAAAAAAA", 0)) !=
            hamming1_array.end());
    REQUIRE(std::find(hamming1_array.begin(), hamming1_array.end(), to_uint64("AAAAAAAAAAAAAAAATAAAAAAAAAAAAAAA", 0)) !=
            hamming1_array.end());
    REQUIRE(std::find(hamming1_array.begin(), hamming1_array.end(), to_uint64("AAAAAAAAAAAAAAATAAAAAAAAAAAAAAAA", 0)) !=
            hamming1_array.end());
    REQUIRE(std::find(hamming1_array.begin(), hamming1_array.end(), to_uint64("TAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 0)) !=
            hamming1_array.end());
    REQUIRE(std::find(hamming1_array.begin(), hamming1_array.end(), to_uint64("AAAAGAAAAAAAAAAAAAAAAAAAAAAAAAAA", 0)) !=
            hamming1_array.end());
    REQUIRE(std::find(hamming1_array.begin(), hamming1_array.end(), to_uint64("AAAAAAAAAAACAAAAAAAAAAAAAAAAAAAA", 0)) !=
            hamming1_array.end());

    // The array should never contain the same 32-mer twice
    std::sort(hamming1_array.begin(), hamming1_array.end());
    REQUIRE(std::adjacent_find(hamming1_array.begin(), hamming1_array.end()) == hamming1_array.end());
  }
}
