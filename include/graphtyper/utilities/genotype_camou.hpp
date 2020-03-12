#pragma once

#include <string>
#include <vector>

#include <graphtyper/graph/genomic_region.hpp>


namespace gyper
{

void
genotype_camou(std::string const & interval_fn,
               std::string const & ref_path,
               std::vector<std::string> const & sams,
               std::string const & output_path,
               std::vector<double> const & avg_cov_by_readlen);

} // namespace gyper
