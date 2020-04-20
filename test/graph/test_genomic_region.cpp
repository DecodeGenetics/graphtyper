#include <string> //std::string

#include <boost/log/trivial.hpp>

#include <graphtyper/graph/absolute_position.hpp>
#include <graphtyper/graph/genomic_region.hpp>
#include <graphtyper/graph/graph.hpp>
#include <graphtyper/graph/label.hpp>
#include <graphtyper/utilities/type_conversions.hpp>

#include <catch.hpp>

/*
TEST_CASE("Getting absolute and contig positions from genomic region")
{
  BOOST_LOG_TRIVIAL(debug) << "[" << __HERE__
                           << "] TEST CASE: Getting absolute and contig positions from genomic region.";

  gyper::GenomicRegion genomic_region;

  uint32_t const CHR01_LENGTH = 66ul;
  uint32_t const CHR02_LENGTH = 66ul;
  uint32_t const CHR03_LENGTH = 66ul;

  // Get absolute position
  {
    REQUIRE(genomic_region.get_absolute_position("chr1", 1) == 1);
    REQUIRE(genomic_region.get_absolute_position("chr1", 100) == 100);
    REQUIRE(genomic_region.get_absolute_position("chr2", 100) == 100 + CHR01_LENGTH);
    REQUIRE(genomic_region.get_absolute_position("chr4", 1) == 1 + CHR01_LENGTH + CHR02_LENGTH + CHR03_LENGTH);
  }

  // Get contig position
  {
    REQUIRE(genomic_region.get_contig_position(1, gyper::graph).first == "chr1");
    REQUIRE(genomic_region.get_contig_position(1, gyper::graph).second == 1ul);

    REQUIRE(genomic_region.get_contig_position(3, gyper::graph).first == "chr1");
    REQUIRE(genomic_region.get_contig_position(3, gyper::graph).second == 3ul);

    REQUIRE(genomic_region.get_contig_position(1 + CHR01_LENGTH, gyper::graph).first == "chr2");
    REQUIRE(genomic_region.get_contig_position(1 + CHR01_LENGTH, gyper::graph).second == 1ul);

    REQUIRE(genomic_region.get_contig_position(CHR01_LENGTH, gyper::graph).first == "chr1");
    REQUIRE(genomic_region.get_contig_position(CHR01_LENGTH, gyper::graph).second == CHR01_LENGTH);
  }
}
*/
