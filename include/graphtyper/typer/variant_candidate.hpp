#pragma once

#include <cstdint>
#include <string> // std::string
#include <vector> // std::vector

#include <cereal/access.hpp>

#include <graphtyper/typer/variant.hpp>

namespace gyper
{
/**
Variant candidate is almost the same as Variant, but I was not sure how to use polymorphism effectively here,
so I ended up repeating some code in both classes.
*/

class VariantCandidate
{
  friend class cereal::access;

public:
  uint32_t abs_pos = 0;
  uint32_t original_pos = 0u;
  std::vector<std::vector<char>> seqs;
  uint16_t flags = 0;

  /******************
   * CLASS MODIFERS *
   ******************/
  bool add_base_in_front(bool add_N = false);
  bool add_base_in_back(bool add_N = false);
  void expanded_normalized();
  void normalize(); /** \brief Defined here: http://genome.sph.umich.edu/wiki/Variant_Normalization */

  /*********************
   * CLASS INFORMATION *
   ********************/
  bool is_normalized() const;
  bool is_snp_or_snps() const;

  // returns 0 for false, 1 for transition, 2 for transversion
  int is_transition_or_transversion() const;
  std::string print() const;

  /*********************
   * OPERATOR OVERLOAD *
   *********************/
  bool operator==(VariantCandidate const & b) const;
  bool operator!=(VariantCandidate const & b) const;
  bool operator<(VariantCandidate const & b) const;

private:
  template <class Archive>
  void serialize(Archive & ar, const unsigned int)
  {
    ar & abs_pos;
    ar & original_pos;
    ar & flags;
    ar & seqs;
  }
};

// Hash function for VariantCandidate
struct VariantCandidateHash
{
  std::size_t operator()(VariantCandidate const & v) const;
};

} // namespace gyper
