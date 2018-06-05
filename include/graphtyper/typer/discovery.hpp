#pragma once

#include <vector> // std::vector
#include <utility> // std::pair

#include <graphtyper/typer/genotype_paths.hpp> // gyper::GenotypePaths
#include <graphtyper/typer/variant_candidate.hpp> // gyper::VariantCandidate


namespace gyper
{

std::vector<VariantCandidate>
discover_variants(std::vector<GenotypePaths> const & genos);

std::vector<VariantCandidate>
discover_variants(std::vector<std::pair<GenotypePaths, GenotypePaths> > const & genos);

} // namespace gyper
