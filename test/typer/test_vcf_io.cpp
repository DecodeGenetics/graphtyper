#include <catch.hpp>

#include <stdio.h>
#include <climits>
#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>

#include <graphtyper/constants.hpp>
#include <graphtyper/typer/vcf.hpp>


TEST_CASE("Read the index test VCF file")
{
  using gyper::Vcf;
  using gyper::READ_UNCOMPRESSED_MODE;

  std::string vcf_filename;

  {
    std::stringstream vcf_ss;
    vcf_ss << gyper_SOURCE_DIRECTORY << "/test/data/reference/index_test.vcf";
    vcf_filename = vcf_ss.str();
  }

  Vcf vcf(READ_UNCOMPRESSED_MODE, vcf_filename);
  vcf.read();

  SECTION("This VCF does now has all variants")
  {
    REQUIRE(vcf.variants.size() == 5);
    // chr1	37	rs1	C	G	0	.	.
    // chr2	2	rs2	C	A	0	.	.
    // chr2	3	rs3	C	A	0	.	.
    // chr3	31	rs4	A	G,GA	0	.	.
    // chr4	2	.	A	T	0	.	.

    auto const & vars = vcf.variants;
    REQUIRE(vars[0].abs_pos == 37);

    REQUIRE(vars[0].seqs.size() == 2);
    REQUIRE(vars[1].seqs.size() == 2);
    REQUIRE(vars[2].seqs.size() == 2);
    REQUIRE(vars[3].seqs.size() == 3);
    REQUIRE(vars[4].seqs.size() == 2);
  }

  SECTION("This VCF does not have any samples")
  {
    REQUIRE(vcf.sample_names.size() == 0);
  }
}
