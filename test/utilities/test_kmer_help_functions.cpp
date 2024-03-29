#include <climits>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/graph_serialization.hpp>
#include <graphtyper/index/kmer_label.hpp>
#include <graphtyper/utilities/kmer_help_functions.hpp>
#include <graphtyper/utilities/type_conversions.hpp>

#include <catch2/catch.hpp>

TEST_CASE("Get the number of kmers in a dna string")
{
  using namespace gyper;

  SECTION("K == the length of dna string")
  {
    REQUIRE(get_num_kmers(seqan::Dna5String("AAAACAAAAGAAAACAAAAGAAAACAAAAGAT")) == 1);
  }

  SECTION("Exactly 2 kmers (32 + 31 = 63)")
  {
    REQUIRE(get_num_kmers(seqan::Dna5String("AAAACAAAAGAAAACAAAAGAAAACAAAAGATAAAACAAAAGAAAACAAAAGAAAACAAAAG")) == 1);
    REQUIRE(get_num_kmers(seqan::Dna5String("AAAACAAAAGAAAACAAAAGAAAACAAAAGATAAAACAAAAGAAAACAAAAGAAAACAAAAGA")) == 2);
    REQUIRE(get_num_kmers(seqan::Dna5String("AAAACAAAAGAAAACAAAAGAAAACAAAAGATAAAACAAAAGAAAACAAAAGAAAACAAAAGAT")) == 2);
  }

  SECTION("3 kmers (32 + 31 + 31 = 94)")
  {
    REQUIRE(get_num_kmers(seqan::Dna5String(
              "AAAACAAAAGAAAACAAAAGAAAACAAAAGATAAAACAAAAGAAAACAAAAGAAAACAAAAGATAAAACAAAAGAAAACAAAAGAAAACAAAA")) == 2);
    REQUIRE(get_num_kmers(seqan::Dna5String(
              "AAAACAAAAGAAAACAAAAGAAAACAAAAGATAAAACAAAAGAAAACAAAAGAAAACAAAAGATAAAACAAAAGAAAACAAAAGAAAACAAAAG")) == 3);
    REQUIRE(get_num_kmers(seqan::Dna5String(
              "AAAACAAAAGAAAACAAAAGAAAACAAAAGATAAAACAAAAGAAAACAAAAGAAAACAAAAGATAAAACAAAAGAAAACAAAAGAAAACAAAAGA")) == 3);
  }
}

TEST_CASE("Get the ith kmers in a dna string")
{
  using namespace gyper;

  SECTION("K == the length of dna string")
  {
    seqan::Dna5String kmer1 = "AAAACAAAAGAAAACAAAAGAAAACAAAAGAT";
    REQUIRE(get_ith_kmer(kmer1, 0u) == kmer1);
    REQUIRE(get_ith_kmer(seqan::Dna5String("AAAACAAAAGAAAACAAAAGAAAACAAAAGATT"), 0u) == kmer1);
    REQUIRE(get_ith_kmer(seqan::Dna5String("AAAACAAAAGAAAACAAAAGAAAACAAAAGATTT"), 0u) ==
            seqan::Dna5String("AAACAAAAGAAAACAAAAGAAAACAAAAGATT"));
  }

  SECTION("2 kmers (32 + 31 = 63)")
  {
    REQUIRE(get_ith_kmer(seqan::Dna5String("AAAACAAAAGAAACCAAAAGAAAACAAAAGATAAAACAAAAGAAAACAAAAGAAAACAAAAG"), 0u) ==
            seqan::Dna5String("AAAAGAAAACAAAAGATAAAACAAAAGAAAAC"));
    REQUIRE(get_ith_kmer(seqan::Dna5String("AAAACAAAAGAAACCAAAAGAAAACAAAAGATAAAACAAAAGAAAACAAAAGAAAACAAAAGA"), 0u) ==
            seqan::Dna5String("AAAACAAAAGAAACCAAAAGAAAACAAAAGAT"));
    REQUIRE(get_ith_kmer(seqan::Dna5String("AAAACAAAAGAAACCAAAAGAAAACAAAAGATAAAACAAAAGAAAACAAAAGAAAACAAAAGA"), 1u) ==
            seqan::Dna5String("TAAAACAAAAGAAAACAAAAGAAAACAAAAGA"));
  }
}

TEST_CASE("Get multiple uint64_t with Iupac reads")
{
  SECTION("Non ACGT")
  {
    seqan::IupacString read1 = "ACCGGGGTTAAAATTGAAAACCCCTAAAATTGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    std::vector<uint64_t> keys = gyper::to_uint64_vec(read1, 0);

    REQUIRE(keys.size() == 1);
    REQUIRE(gyper::to_dna_str(keys[0]) == "ACCGGGGTTAAAATTGAAAACCCCTAAAATTG");

    keys = gyper::to_uint64_vec(read1, 10);
    REQUIRE(keys.size() == 1);
    REQUIRE(gyper::to_dna_str(keys[0]) == "AAATTGAAAACCCCTAAAATTGAAAAAAAAAA");
  }

  SECTION("One and two non ACGT")
  {
    seqan::IupacString read1 =
      "ACCGGGGTTAAAATTGAAAACCCCTAAAATTNAAAAAAAAAAAAAAAAAAAAAAAAAWAAAAAAAAAATTTTTTTBTTTTTTTTTTTTTTTTTTT";
    std::vector<uint64_t> keys = gyper::to_uint64_vec(read1, 0);

    REQUIRE(keys.size() == 4);
    REQUIRE(gyper::to_dna_str(keys[0]) == "ACCGGGGTTAAAATTGAAAACCCCTAAAATTT");
    REQUIRE(gyper::to_dna_str(keys[1]) == "ACCGGGGTTAAAATTGAAAACCCCTAAAATTA");
    REQUIRE(gyper::to_dna_str(keys[2]) == "ACCGGGGTTAAAATTGAAAACCCCTAAAATTC");
    REQUIRE(gyper::to_dna_str(keys[3]) == "ACCGGGGTTAAAATTGAAAACCCCTAAAATTG");

    keys = gyper::to_uint64_vec(read1, 32);
    REQUIRE(keys.size() == 2);
    REQUIRE(gyper::to_dna_str(keys[0]) == "AAAAAAAAAAAAAAAAAAAAAAAAATAAAAAA");
    REQUIRE(gyper::to_dna_str(keys[1]) == "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");

    keys = gyper::to_uint64_vec(read1, 63);
    REQUIRE(keys.size() == 3);
    REQUIRE(gyper::to_dna_str(keys[0]) == "AAAAATTTTTTTTTTTTTTTTTTTTTTTTTTT");
    REQUIRE(gyper::to_dna_str(keys[1]) == "AAAAATTTTTTTCTTTTTTTTTTTTTTTTTTT");
    REQUIRE(gyper::to_dna_str(keys[2]) == "AAAAATTTTTTTGTTTTTTTTTTTTTTTTTTT");
  }

  SECTION("High amounts of Ns result in no keys")
  {
    seqan::IupacString read1 = "NNNNNNNNNNNNAAAAAAAAAAAAAAAAAAAAAA";
    std::vector<uint64_t> keys = gyper::to_uint64_vec(read1, 0);
    REQUIRE(keys.size() == 0);
  }
}
