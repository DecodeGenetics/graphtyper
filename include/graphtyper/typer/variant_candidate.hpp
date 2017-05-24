#pragma once

#include <string> // std::string
#include <vector> // std::vector

#include <graphtyper/typer/variant.hpp>

namespace gyper
{

/**
Variant candidate is almost the same as Variant, but I was not sure how to use polymorphism effectively here,
so I ended up repeating some code in both classes.
*/

class VariantCandidate
{
public:
  uint32_t abs_pos = 0;
  std::vector<std::vector<char> > seqs;
  bool is_low_qual = false;

  /******************
   * CLASS MODIFERS *
   ******************/
  bool add_base_in_front(bool const add_N = false);
  bool add_base_in_back(bool const add_N = false);
  void normalize(); /** \brief Defined here: http://genome.sph.umich.edu/wiki/Variant_Normalization */

  /*********************
   * CLASS INFORMATION *
   ********************/
  bool is_normalized() const;
  bool is_snp_or_snps() const;
  std::string print() const;

  /*********************
   * OPERATOR OVERLOAD *
   *********************/
  bool operator==(VariantCandidate const & b) const;
  bool operator!=(VariantCandidate const & b) const;
  bool operator<(VariantCandidate const & b) const;
};

// Hash function for VariantCandidate
struct VariantCandidateHash
{
  std::size_t operator()(VariantCandidate const & v) const;
};

} // namespace gyper
