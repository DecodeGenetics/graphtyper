#pragma once

#include <cstdint> // int32_t
#include <string> // std::string

#include <boost/serialization/access.hpp>

#include <graphtyper/graph/alu_sequences.hpp>


enum SVTYPE {NOT_SV = 0, DEL, DEL_ALU, DUP, INS, INS_ALU, INV, OTHER};


class SV
{
  friend class boost::serialization::access; // boost is my friend

public:
  SVTYPE type = NOT_SV;
  int32_t length = 0;  // Length of alt. allele minus ref. allele
  int32_t size = 0;  // Total size, i.e. length of larger allele excluding padding base
  int32_t end = 0;
  int16_t n_clusters = 0;
  int32_t or_start = -1; // Start coordinate of the sequence origin
  int32_t or_end = -1; // End coordinate of the sequence origin
  std::vector<char> seq{};

  SV() = default;

  std::string get_type() const;
  std::string get_allele() const;

private:
  template <typename Archive>
  void serialize(Archive & ar, unsigned int);
};
