#pragma once

#include <cstdint> // uint32_t
#include <set> // std::set<T>

#include <cereal/access.hpp>


namespace gyper
{

class VariantSupport
{
  friend class cereal::access;

public:
  uint16_t hq_support{0u};
  uint16_t lq_support{0u};
  uint16_t proper_pairs{0u};
  uint16_t depth{0u};
  uint16_t first_in_pairs{0u};
  uint16_t sequence_reversed{0u};
  uint16_t clipped{0u};
  uint16_t var_size{0u};
  uint16_t growth{0u};
  std::set<uint32_t> unique_positions;
  bool is_indel{false};
  bool is_any_mapq_good{false};

#ifndef NDEBUG
  uint32_t pn_index{0u};
#endif

  VariantSupport() = default;

  void set_depth(uint16_t _depth);
  long get_score() const;
  double get_corrected_support() const;
  double get_ratio() const;
  bool is_above_cutoff(long const min_support, double const min_ratio) const;
  bool is_ratio_above_cutoff(double const min_ratio) const;
  bool is_support_above_cutoff(long const min_support) const;

private:
  template <class Archive>
  void
  serialize(Archive & ar, const unsigned int)
  {
    ar & hq_support;
    ar & lq_support;
    ar & proper_pairs;
    ar & depth;
    ar & first_in_pairs;
    ar & sequence_reversed;
    ar & clipped;
    ar & is_indel;
    ar & var_size;
    ar & growth;

// We only need to serialize things that are used in "filter_varmap_for_all"
//    ar & unique_positions;
//    ar & is_any_mapq_good;

#ifndef NDEBUG
    ar & pn_index;
#endif
  }


};

} // namespace gyper
