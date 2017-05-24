#pragma once

#include <array>
#include <string>
#include <unordered_map>
#include <vector>


namespace gyper
{

class AbsolutePosition
{
public:
  AbsolutePosition();

  uint32_t get_absolute_position(std::string const & chromosome, uint32_t const contig_position) const;
  std::pair<std::string, uint32_t> get_contig_position(uint32_t const absolute_position) const;

  std::array<uint32_t, 24ul> offsets;
  std::unordered_map<std::string, uint32_t> chromosome_to_offset;

  // Accumulate
  uint32_t const CHR01_LENGTH = 248956422ul; // 248 956 422
  uint32_t const CHR02_LENGTH = 242193529ul; // 491 149 951
  uint32_t const CHR03_LENGTH = 198295559ul; // 689 445 510
  uint32_t const CHR04_LENGTH = 190214555ul;
  uint32_t const CHR05_LENGTH = 181538259ul;
  uint32_t const CHR06_LENGTH = 170805979ul;
  uint32_t const CHR07_LENGTH = 159345973ul;
  uint32_t const CHR08_LENGTH = 145138636ul;
  uint32_t const CHR09_LENGTH = 138394717ul;
  uint32_t const CHR10_LENGTH = 133797422ul;
  uint32_t const CHR11_LENGTH = 135086622ul;
  uint32_t const CHR12_LENGTH = 133275309ul;
  uint32_t const CHR13_LENGTH = 114364328ul;
  uint32_t const CHR14_LENGTH = 107043718ul;
  uint32_t const CHR15_LENGTH = 101991189ul;
  uint32_t const CHR16_LENGTH = 90338345ul;
  uint32_t const CHR17_LENGTH = 83257441ul;
  uint32_t const CHR18_LENGTH = 80373285ul;
  uint32_t const CHR19_LENGTH = 58617616ul;  // 2'713'028'904 2713028904ull
  uint32_t const CHR20_LENGTH = 64444167ul;  // 2'777'473'071 2777473071ull
  uint32_t const CHR21_LENGTH = 46709983ul;
  uint32_t const CHR22_LENGTH = 50818468ul;
  uint32_t const CHR0X_LENGTH = 156040895ul;
  uint32_t const CHR0Y_LENGTH = 57227415ul;  //
  // 3'088'269'832

  std::array<uint32_t, 24ul> const chromosome_lengths =
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
    CHR0Y_LENGTH
}};

  std::array<std::string, 24ul> const chromosome_names =
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
    "chrY"
}};

};

extern AbsolutePosition absolute_pos;

} // namespace gyper
