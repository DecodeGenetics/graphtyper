#include <string> //std::string

#include <graphtyper/graph/absolute_position.hpp>
#include <graphtyper/graph/genomic_region.hpp>
#include <graphtyper/graph/graph.hpp>
#include <graphtyper/graph/label.hpp>
#include <graphtyper/utilities/type_conversions.hpp>

#include <catch.hpp>

TEST_CASE("Getting absolute and contig positions from genomic region")
{
  gyper::GenomicRegion genomic_region;

  SECTION("Get absolute position")
  {
    REQUIRE(genomic_region.get_absolute_position("chr1", 1) == 1);
    REQUIRE(genomic_region.get_absolute_position("chr1", 100) == 100);
    REQUIRE(genomic_region.get_absolute_position("chr2", 100) == 100 + gyper::CHR01_LENGTH);
    REQUIRE(genomic_region.get_absolute_position("chr4", 1) == 1 + gyper::CHR01_LENGTH + gyper::CHR02_LENGTH + gyper::CHR03_LENGTH);
  }

  SECTION("Get contig position")
  {
    REQUIRE(genomic_region.get_contig_position(1).first == "chr1");
    REQUIRE(genomic_region.get_contig_position(1).second == 1ul);

    REQUIRE(genomic_region.get_contig_position(100).first == "chr1");
    REQUIRE(genomic_region.get_contig_position(100).second == 100ul);

    REQUIRE(genomic_region.get_contig_position(1 + gyper::CHR01_LENGTH).first == "chr2");
    REQUIRE(genomic_region.get_contig_position(1 + gyper::CHR01_LENGTH).second == 1ul);

    REQUIRE(genomic_region.get_contig_position(gyper::CHR01_LENGTH).first == "chr1");
    REQUIRE(genomic_region.get_contig_position(gyper::CHR01_LENGTH).second == gyper::CHR01_LENGTH);

    uint32_t const abs_pos_chrY = genomic_region.get_absolute_position("chrY", gyper::CHR0Y_LENGTH);
    REQUIRE(genomic_region.get_contig_position(abs_pos_chrY).first == "chrY");
    REQUIRE(genomic_region.get_contig_position(abs_pos_chrY).second == gyper::CHR0Y_LENGTH);
  }
}
