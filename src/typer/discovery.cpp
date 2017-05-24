#include <vector>

#include <graphtyper/typer/discovery.hpp>
#include <graphtyper/typer/genotype_paths.hpp>
#include <graphtyper/typer/variant_candidate.hpp>


namespace gyper
{

std::vector<VariantCandidate>
discover_variants(std::vector<GenotypePaths> const & genos)
{
  std::vector<VariantCandidate> variants;

  for (auto const & geno : genos)
  {
    // Require all matches to be unique
    if (geno.all_paths_unique())
    {
      std::vector<VariantCandidate> new_vars = geno.find_new_variants();
      std::move(new_vars.begin(), new_vars.end(), std::back_inserter(variants));
    }
  }

  return variants;
}


std::vector<VariantCandidate>
discover_variants(std::vector<std::pair<GenotypePaths, GenotypePaths> > const & genos)
{
  std::vector<VariantCandidate> variants;

  for (auto const & geno : genos)
  {
    // Require all matches to be unique
    if (geno.first.all_paths_unique())
    {
      std::vector<VariantCandidate> new_vars = geno.first.find_new_variants();
      std::move(new_vars.begin(), new_vars.end(), std::back_inserter(variants));
    }

    if (geno.second.all_paths_unique())
    {
      std::vector<VariantCandidate> new_vars = geno.second.find_new_variants();
      std::move(new_vars.begin(), new_vars.end(), std::back_inserter(variants));
    }
  }

  return variants;
}


} // namespace gyper
