#pragma once

#include <cstdint> // int32_t
#include <string> // std::string

#include <boost/serialization/access.hpp>

#include <graphtyper/graph/alu_sequences.hpp>


namespace gyper
{

class Variant;


enum SVTYPE
{
  NOT_SV = 0, DEL, DEL_ALU, DUP, INS, INS_ALU, INV, BND, OTHER
};


enum INVTYPE
{
  NOT_INV = 0, INV3, INV5, BOTH_BREAKPOINTS
};


class SV
{
  friend class boost::serialization::access; // boost is my friend

public:
  SVTYPE type = NOT_SV;
  std::string chrom = "";
  int32_t begin = 0;
  int32_t length = 0;  // Length of alt. allele minus ref. allele
  int32_t size = 0;  // Total size, i.e. length of larger allele excluding padding base
  int32_t end = 0;
  int16_t n_clusters = 0;
  int16_t num_merged_svs = -1;
  int32_t or_start = -1; // Start coordinate of the sequence origin
  int32_t or_end = -1; // End coordinate of the sequence origin
  int32_t related_sv = -1;
  std::string model = "AGGREGATED";
  INVTYPE inv_type = NOT_INV;
  std::vector<char> seq{};
  std::vector<char> hom_seq{};
  std::vector<char> ins_seq{};
  std::vector<char> ins_seq_left{};
  std::vector<char> ins_seq_right{};
  std::vector<char> original_alt{};

  SV() = default;
  std::string get_type() const;
  std::string get_allele() const;

private:
  template<typename Archive>
  void serialize(Archive &ar, unsigned int);
};


void reformat_sv_vcf_records(std::vector<Variant> & variant);

} // namespace gyper