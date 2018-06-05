#include <cstdint>

#include <graphtyper/utilities/options.hpp>
#include <graphtyper/typer/variant_support.hpp>


namespace gyper
{


void
VariantSupport::set_depth(uint16_t const _depth)
{
  // Make sure that the depth is never reduced
  depth = std::max(depth, _depth);
}

double
VariantSupport::get_ratio() const
{
  double const effective_depth = static_cast<int32_t>(depth) - static_cast<int32_t>(lq_support / 2);

  if (effective_depth < 0.9) // Should never happend, but just in case there is some extreme edge case
    return 1.0;
  else
    return static_cast<double>(get_support()) / effective_depth;
}


uint32_t
VariantSupport::get_support() const
{
  return hq_support + lq_support / 2u;
}


bool
VariantSupport::is_ratio_above_cutoff() const
{
  return get_ratio() > Options::instance()->minimum_variant_support_ratio;
}


bool
VariantSupport::is_support_above_cutoff() const
{
  return unique_positions.size() > 1 &&
         proper_pairs > 2 &&
         hq_support > 1 &&
         first_in_pairs > 0 &&
         first_in_pairs < depth &&
         sequence_reversed > 0 &&
         sequence_reversed < depth &&
         get_support() >= Options::instance()->minimum_variant_support;
}


bool
VariantSupport::is_above_cutoff() const
{
  return is_support_above_cutoff() && is_ratio_above_cutoff();
}

} // namespace gyper
