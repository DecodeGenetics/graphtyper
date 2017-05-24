#pragma once

#include <unordered_map>
#include <vector>

#include <boost/optional.hpp>
#include <boost/optional/optional_io.hpp>
#include <boost/unordered/unordered_map.hpp>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/forward_list.hpp>

#include <graphtyper/constants.hpp>


namespace gyper
{

class KmerLabel
{
  friend class boost::serialization::access;

public:
  uint32_t start_index; /** \brief The index where the variant starts on the reference genome. */
  uint32_t end_index;   /** \brief The index where the variant ends on the reference genome. */
  uint32_t variant_id = INVALID_ID;  /** \brief The variant id, i.e. the index of the variant on the graph. */
  uint16_t variant_num = INVALID_NUM; /** \brief The number of the variant (filled when fetching the K-mer from the index) */
  uint32_t variant_order = INVALID_ID; /** \brief The order of the variant (filled when fetching the K-mer from the index) */

  KmerLabel() noexcept
    : start_index(INVALID_ID)
    , end_index(INVALID_ID)
    , variant_id(INVALID_ID)
  {}

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

  // Used when creating a new KmerLabel from the graph.
  KmerLabel(uint32_t const s, uint32_t const e, uint32_t const i, uint16_t const n, uint32_t const o) noexcept
    : start_index(s)
    , end_index(e)
    , variant_id(i)
    , variant_num(n)
    , variant_order(o)
  {}

  bool
  operator==(KmerLabel const & c2) const
  {
    return start_index == c2.start_index &&
           end_index == c2.end_index &&
           variant_id == c2.variant_id &&
           variant_num == c2.variant_num &&
           variant_order == c2.variant_order
    ;
  }


  bool
  operator!=(KmerLabel const & c2) const
  {
    return !(*this == c2);
  }


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
