#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>

#include <boost/serialization/access.hpp>

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/var_record.hpp>

namespace gyper
{

static const uint32_t AS_LONG_AS_POSSIBLE = 0xFFFFFFFFULL;


class GenomicRegion
{
  friend class boost::serialization::access;

public:
  uint16_t rID;
  std::string chr;
  uint32_t begin;
  uint32_t end;
  uint32_t region_to_refnode;

  // Not serialized
  std::array<uint32_t, 24ul> offsets;

  /****************
   * CONSTRUCTORS *
   ****************/
  GenomicRegion(uint32_t _region_to_refnode = 0);
  GenomicRegion(std::string region, uint32_t _region_to_refnode = 0);
  GenomicRegion(uint16_t && r, std::string && c, uint32_t && b, uint32_t && e, uint32_t _region_to_refnode = 0);

  uint32_t get_absolute_begin_position() const;
  uint32_t get_absolute_end_position() const;

  uint32_t get_absolute_position(std::string const & chromosome, uint32_t contig_position) const;
  uint32_t get_absolute_position(uint32_t contig_position) const;
  std::pair<std::string, uint32_t> get_contig_position(uint32_t absolute_position) const;

  std::string to_string() const;

  void check_if_var_records_match_reference_genome(std::vector<VarRecord> const & var_records, std::vector<char> const & reference);
  void add_reference_to_record_if_they_have_a_matching_prefix(VarRecord & var_record, std::vector<char> const & reference);

private:
  template <typename Archive>
  void serialize(Archive & ar, const unsigned int);
  void gather_contig_positions();

  std::unordered_map<std::string, uint32_t> chromosome_to_offset;
};

// Accumulate
static const uint32_t CHR01_LENGTH = 248956422ul; // 248956422
static const uint32_t CHR02_LENGTH = 242193529ul; // 491149951
static const uint32_t CHR03_LENGTH = 198295559ul; // 689445510
static const uint32_t CHR04_LENGTH = 190214555ul; // 879660065
static const uint32_t CHR05_LENGTH = 181538259ul; // 1061198324
static const uint32_t CHR06_LENGTH = 170805979ul; // 1232004303
static const uint32_t CHR07_LENGTH = 159345973ul;
static const uint32_t CHR08_LENGTH = 145138636ul;
static const uint32_t CHR09_LENGTH = 138394717ul;
static const uint32_t CHR10_LENGTH = 133797422ul;
static const uint32_t CHR11_LENGTH = 135086622ul;
static const uint32_t CHR12_LENGTH = 133275309ul;
static const uint32_t CHR13_LENGTH = 114364328ul;
static const uint32_t CHR14_LENGTH = 107043718ul;
static const uint32_t CHR15_LENGTH = 101991189ul;
static const uint32_t CHR16_LENGTH = 90338345ul;
static const uint32_t CHR17_LENGTH = 83257441ul;
static const uint32_t CHR18_LENGTH = 80373285ul;
static const uint32_t CHR19_LENGTH = 58617616ul; // 2'713'028'904 2713028904ul
static const uint32_t CHR20_LENGTH = 64444167ul; // 2'777'473'071 2777473071ul
static const uint32_t CHR21_LENGTH = 46709983ul;
static const uint32_t CHR22_LENGTH = 50818468ul;
static const uint32_t CHR0X_LENGTH = 156040895ul;
static const uint32_t CHR0Y_LENGTH = 57227415ul; // 3'088'269'832
static const uint32_t CHR0M_LENGTH = 16569;
static const uint32_t CHRUN_LENGTH = 400000000ul;
// Special position start is at 3489660928

static const std::array<uint32_t, 26ul> chromosome_lengths =
{{
  CHR01_LENGTH,
  CHR02_LENGTH,
  CHR03_LENGTH,
  CHR04_LENGTH,
  CHR05_LENGTH,
  CHR06_LENGTH,
  CHR07_LENGTH,
  CHR08_LENGTH,
  CHR09_LENGTH,
  CHR10_LENGTH,
  CHR11_LENGTH,
  CHR12_LENGTH,
  CHR13_LENGTH,
  CHR14_LENGTH,
  CHR15_LENGTH,
  CHR16_LENGTH,
  CHR17_LENGTH,
  CHR18_LENGTH,
  CHR19_LENGTH,
  CHR20_LENGTH,
  CHR21_LENGTH,
  CHR22_LENGTH,
  CHR0X_LENGTH,
  CHR0Y_LENGTH,
  CHR0M_LENGTH,
  CHRUN_LENGTH
}};


static const std::array<std::string, 26ul> chromosome_names =
{{
  "chr1",
  "chr2",
  "chr3",
  "chr4",
  "chr5",
  "chr6",
  "chr7",
  "chr8",
  "chr9",
  "chr10",
  "chr11",
  "chr12",
  "chr13",
  "chr14",
  "chr15",
  "chr16",
  "chr17",
  "chr18",
  "chr19",
  "chr20",
  "chr21",
  "chr22",
  "chrX",
  "chrY",
  "chrM",
  "chrUn"
}};

} // namespace gyper
