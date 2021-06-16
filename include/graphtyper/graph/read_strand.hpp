#pragma once

#include <cstdint>

#include <cereal/access.hpp>

namespace gyper
{
class ReadStrand
{
  friend class cereal::access;

public:
  uint32_t r1_forward{0u};
  uint32_t r1_reverse{0u};
  uint32_t r2_forward{0u};
  uint32_t r2_reverse{0u};

  long get_weight() const;
  long get_r1_count() const;
  long get_reverse_count() const;
  long get_max_bias() const;

  void merge_with(ReadStrand const & other_rs);

  std::string str() const;

private:
  template <class Archive>
  void serialize(Archive & ar, unsigned int version);
};

} // namespace gyper
