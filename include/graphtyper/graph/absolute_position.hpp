#pragma once

#include <array>
#include <string>
#include <unordered_map>
#include <vector>


namespace gyper
{

// Maximum allowed length for each chromosomes (current lengths allow both GRChb37 and GRCh38)
uint32_t static const CHR01_LENGTH = 249250621ul;
uint32_t static const CHR02_LENGTH = 243199373ul;
uint32_t static const CHR03_LENGTH = 198295559ul;
uint32_t static const CHR04_LENGTH = 191154276ul;
uint32_t static const CHR05_LENGTH = 181538259ul;
uint32_t static const CHR06_LENGTH = 171115067ul;
uint32_t static const CHR07_LENGTH = 159345973ul;
uint32_t static const CHR08_LENGTH = 146364022ul;
uint32_t static const CHR09_LENGTH = 141213431ul;
uint32_t static const CHR10_LENGTH = 135534747ul;
uint32_t static const CHR11_LENGTH = 135086622ul;
uint32_t static const CHR12_LENGTH = 133851895ul;
uint32_t static const CHR13_LENGTH = 115169878ul;
uint32_t static const CHR14_LENGTH = 107349540ul;
uint32_t static const CHR15_LENGTH = 102531392ul;
uint32_t static const CHR16_LENGTH = 90354753ul;
uint32_t static const CHR17_LENGTH = 83257441ul;
uint32_t static const CHR18_LENGTH = 80373285ul;
uint32_t static const CHR19_LENGTH = 59128983ul;
uint32_t static const CHR20_LENGTH = 64444167ul;
uint32_t static const CHR21_LENGTH = 48129895ul;
uint32_t static const CHR22_LENGTH = 51304566ul;
uint32_t static const CHR0X_LENGTH = 156040895ul;
uint32_t static const CHR0Y_LENGTH = 59373566ul;
uint32_t static const CHR0M_LENGTH = 16569;
uint32_t static const CHRUN_LENGTH = 400000000ul;
// Special position start after all the other chromosomes

std::array<uint32_t, 26u> const chromosome_lengths =
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

std::array<std::string, 26u> const chromosome_names =
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


class AbsolutePosition
{
public:
  AbsolutePosition();

  uint32_t get_absolute_position(std::string const & chromosome, uint32_t const contig_position) const;
  std::pair<std::string, uint32_t> get_contig_position(uint32_t const absolute_position) const;

  std::array<uint32_t, 26u> offsets;
  std::unordered_map<std::string, uint32_t> chromosome_to_offset;

};

extern AbsolutePosition absolute_pos;

} // namespace gyper
