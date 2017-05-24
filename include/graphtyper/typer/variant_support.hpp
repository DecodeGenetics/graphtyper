#pragma once

#include <cstdint> // uint32_t

namespace gyper
{

struct VariantSupport
{
  uint16_t support = 0;
  uint16_t lq_support = 0;
  uint16_t depth = 0;

  VariantSupport() = default;
  VariantSupport(uint16_t const _support, uint16_t const _depth)
    : support(_support)
    {
      depth = std::max(_support, _depth); // Make sure depth is never less than the support
    }

  double get_ratio() const;
  bool is_above_cutoff() const;
  bool is_ratio_above_cutoff() const;
  bool is_support_above_cutoff() const;
  bool is_highly_certainly_real() const;
};

} // namespace gyper
