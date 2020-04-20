#pragma once

#include <cstdint>
#include <vector>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <graphtyper/constants.hpp>


namespace gyper
{

class Graph;

class KmerLabel
{
private:
  friend class boost::serialization::access;

public:
  uint32_t start_index{0}; /** \brief The index where the variant starts on the reference genome. */
  uint32_t end_index{0};   /** \brief The index where the variant ends on the reference genome. */
  uint32_t variant_id{INVALID_ID};  /** \brief The variant id, i.e. the index of the variant on the graph. */


  KmerLabel() noexcept = default;

  KmerLabel(uint32_t const s, uint32_t const e) noexcept
    : start_index(s)
    , end_index(e)
    , variant_id(INVALID_ID)
  {}

  // Used when add a KmerLabel to the index.
  KmerLabel(uint32_t const s, uint32_t const e, uint32_t const i) noexcept
    : start_index(s)
    , end_index(e)
    , variant_id(i)
  {}

  // Comparison operators
  bool operator==(KmerLabel const & c2) const;
  bool operator!=(KmerLabel const & c2) const;


private:
  template <class Archive>
  void
  serialize(Archive & ar, const unsigned int)
  {
    ar & start_index;
    ar & end_index;
    ar & variant_id;
  }


};

// Type aliases
using TKmerLabels = std::vector<std::vector<KmerLabel> >;

} // namespace gyper
