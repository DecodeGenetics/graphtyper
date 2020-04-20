#include <cstdint>
#include <iostream>

#include <graphtyper/typer/variant_support.hpp>
#include <graphtyper/utilities/options.hpp>


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
  double const effective_depth = static_cast<long>(depth) -
                                 static_cast<long>(lq_support / 2);

  if (effective_depth < 0.99) // Should never happend, but just in case there is some extreme edge case
  {
    return 1.0;
  }
  else
  {
    double const growth_correction = static_cast<double>(growth / 3 + 10.0) / 10.0;
    return get_corrected_support() / effective_depth * growth_correction;
  }
}


long
VariantSupport::get_score() const
{
  long score = (get_corrected_support() * get_ratio() * 10.0 + 0.5);

  if (hq_support >= 4 && proper_pairs >= 4 && (hq_support + lq_support - clipped >= 3))
    score += 20;

  if (hq_support >= 8 && proper_pairs >= 8 && (hq_support + lq_support - clipped >= 6))
    score += 30;

  return score;
}


double
VariantSupport::get_corrected_support() const
{
  double const correction = static_cast<double>(var_size / 3 + 10.0) / 10.0;
  return correction * (static_cast<double>(hq_support) + static_cast<double>(lq_support) / 2.0);
}


bool
VariantSupport::is_ratio_above_cutoff(double const min_ratio) const
{
  return get_ratio() > min_ratio;
}


bool
VariantSupport::is_support_above_cutoff(long const min_support) const
{
  int const _depth = hq_support + lq_support;
  bool const is_promising = unique_positions.size() >= 3 &&
                            hq_support >= 4 &&
                            proper_pairs >= 3 &&
                            (_depth - clipped >= 3);
  Options const & copts = *(Options::const_instance());

  return (copts.no_filter_on_begin_pos || unique_positions.size() > 1)
         &&
         (!copts.filter_on_mapq || is_any_mapq_good)
         &&
         (!copts.filter_on_proper_pairs || (proper_pairs >= 2 || (proper_pairs >= 1 && is_indel)))
         &&
         (hq_support >= 3 || (hq_support >= 2 && is_indel))
         &&
         (!copts.filter_on_read_bias ||
          is_indel ||
          is_promising ||
          (first_in_pairs > 0 && first_in_pairs < _depth))
         &&
         (!copts.filter_on_strand_bias ||
          is_indel ||
          (is_promising && sequence_reversed > 0 && sequence_reversed < _depth) ||
          (sequence_reversed > 1 && sequence_reversed < (_depth - 1)))
         &&
         (clipped <= (_depth - 3) ||
          (is_indel && clipped <= (_depth - 1)))
         &&
         static_cast<long>(get_corrected_support() + 0.5) >= min_support;
}


bool
VariantSupport::is_above_cutoff(long const min_support, double const min_ratio) const
{
  return is_support_above_cutoff(min_support) && is_ratio_above_cutoff(min_ratio);
}


} // namespace gyper
