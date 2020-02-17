#include <catch.hpp>

#include <stdio.h>
#include <climits>
#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>

#include <graphtyper/graph/graph_serialization.hpp> // load_graph()
#include <graphtyper/index/indexer.hpp> // load_index()
#include <graphtyper/typer/vcf.hpp>


TEST_CASE("Create a VCF and add samples")
{
  using namespace gyper;

  Vcf vcf(WRITE_BGZF_MODE, "-");

  SECTION("Initially there are no samples")
  {
    REQUIRE(vcf.sample_names.size() == 0);
  }

  vcf.sample_names.push_back("TEST_SAMP1");

  SECTION("Here we have added a single sample")
  {
    REQUIRE(vcf.sample_names.size() == 1);
    REQUIRE(vcf.sample_names[0] == std::string("TEST_SAMP1"));
  }

  SECTION("..however, adding another one will")
  {
    vcf.sample_names.push_back("TEST_SAMP2");

    REQUIRE(vcf.sample_names.size() == 2);
    REQUIRE(vcf.sample_names[1] == std::string("TEST_SAMP2"));
  }
}


TEST_CASE("Create a VCF and add variants")
{
  using namespace gyper;

  // Set up
  {
    std::stringstream my_graph;
    my_graph << gyper_SOURCE_DIRECTORY << "/test/data/graphs/index_test_chr2.grf";
    std::stringstream my_index;
    my_index << gyper_SOURCE_DIRECTORY << "/test/data/graphs/index_test_chr2";
    // CCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTTTGGACCC
    // chr2    2       rs2     C       A       0       .       .
    // chr2    3       rs3     C       A       0       .       .

    gyper::load_graph(my_graph.str());
    gyper::load_index(my_index.str());
  }


  std::vector<gyper::Haplotype> haps = graph.get_all_haplotypes(32 /*variant distance*/);
  assert(haps.size() == 1);
  Vcf vcf(WRITE_BGZF_MODE, "-");

  SECTION("Initially there are no variants")
  {
    REQUIRE(vcf.variants.size() == 0);
  }

  vcf.add_haplotype(haps.at(0), false);

  SECTION("Now both variant should have been added")
  {
    REQUIRE(vcf.variants.size() == 2);

    REQUIRE(vcf.variants[0].seqs.size() == 2);
    REQUIRE(vcf.variants[0].seqs[0] == gyper::to_vec("CC"));
    REQUIRE(vcf.variants[0].seqs[1] == gyper::to_vec("CA"));

    REQUIRE(vcf.variants[1].seqs.size() == 2);
    REQUIRE(vcf.variants[1].seqs[0] == gyper::to_vec("CC"));
    REQUIRE(vcf.variants[1].seqs[1] == gyper::to_vec("CA"));
  }
}
