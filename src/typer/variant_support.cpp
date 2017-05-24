#include <graphtyper/utilities/options.hpp>
#include <graphtyper/typer/variant_support.hpp>


namespace gyper
{

double
VariantSupport::get_ratio() const
{
  if (depth == 0) // Should neve happend, but just in case there is some extreme edge case
    return 1.0;
  else
    return static_cast<double>(support) / static_cast<double>(depth);
}


bool
VariantSupport::is_ratio_above_cutoff() const
{
  return get_ratio() > Options::instance()->minimum_variant_support_ratio;
}


bool
VariantSupport::is_support_above_cutoff() const
{
  return support >= Options::instance()->minimum_variant_support;
}


bool
VariantSupport::is_above_cutoff() const
{
  return is_support_above_cutoff() && is_ratio_above_cutoff();
}


bool
VariantSupport::is_highly_certainly_real() const
{
  return support >= Options::instance()->certain_variant_support && get_ratio() >= Options::instance()->certain_variant_support_ratio;
}

} // namespace gyper
