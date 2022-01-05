#include <string>
#include <vector>

#include <graphtyper/graph/graph.hpp>
#include <graphtyper/graph/haplotype.hpp>
#include <graphtyper/graph/label.hpp>
#include <graphtyper/graph/var_record.hpp>
#include <graphtyper/utilities/type_conversions.hpp>

#include <catch2/catch.hpp>

TEST_CASE("Haplotype with one genotype")
{
  using namespace gyper;
  std::vector<char> reference_sequence;
  char testdata[] = "SGTACGEEF";
  reference_sequence.insert(reference_sequence.end(), testdata, testdata + 9);
  std::vector<gyper::VarRecord> records;

  {
    gyper::VarRecord record;
    record.pos = 1;
    record.ref = {'G', 'T', 'A', 'C', 'G'};
    record.alts = {{'G'}};
    records.push_back(record);

    record.pos = 1;
    record.ref = {'G'};
    record.alts = {{'K'}};
    records.push_back(record);
  }

  graph = gyper::Graph();
  graph.add_genomic_region(std::move(reference_sequence), std::move(records), gyper::GenomicRegion());
  graph.create_special_positions();

  std::vector<gyper::Haplotype> haps = graph.get_all_haplotypes();
  REQUIRE(haps.size() == 1);
  REQUIRE(haps[0].get_genotype_num() == 3);
}
