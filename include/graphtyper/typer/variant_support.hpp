#pragma once

#include <cstdint> // uint32_t
#include <set> // std::set<T>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>


namespace gyper
{

class VariantSupport
{
  friend class boost::serialization::access;

public:
  uint16_t hq_support = 0u;
  uint16_t lq_support = 0u;
  uint16_t proper_pairs = 0u;
  uint16_t depth = 0u;
  uint16_t mapq0 = 0u;
  uint16_t unaligned = 0u;
  uint16_t clipped = 0u;
  uint16_t first_in_pairs = 0u;
  uint16_t sequence_reversed = 0u;
  std::set<uint32_t> unique_positions;
  uint32_t pn_index = 0u;

  VariantSupport() = default;

  void set_depth(uint16_t _depth);

  uint32_t get_support() const;
  double get_ratio() const;
  bool is_above_cutoff() const;
  bool is_ratio_above_cutoff() const;
  bool is_support_above_cutoff() const;

private:
  template <class Archive>
  void
  serialize(Archive & ar, const unsigned int)
  {
    ar & hq_support;
    ar & lq_support;
    ar & proper_pairs;
    ar & depth;
    ar & mapq0;
    ar & unaligned;
    ar & clipped;
    ar & first_in_pairs;
    ar & sequence_reversed;
    ar & unique_positions;
    ar & pn_index;
  }
};

} // namespace gyper
