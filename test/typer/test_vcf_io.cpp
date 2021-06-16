#include <climits>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>

#include <graphtyper/constants.hpp>
#include <graphtyper/typer/vcf.hpp>

#include <catch.hpp>

TEST_CASE("Read the index test VCF file")
{
  using gyper::Vcf;

  std::string vcf_filename;

  {
    std::stringstream vcf_ss;
    vcf_ss << gyper_SOURCE_DIRECTORY << "/test/data/reference/index_test.vcf.gz";
    vcf_filename = vcf_ss.str();
  }

  Vcf vcf(gyper::READ_BGZF_MODE, vcf_filename);
  vcf.read();

  SECTION("This VCF does now has all variants")
  {
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
