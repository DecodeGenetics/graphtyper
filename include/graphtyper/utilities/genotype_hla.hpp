#pragma once

#include <string>
#include <vector>

#include <graphtyper/graph/genomic_region.hpp>


namespace gyper
{

void
genotype_hla_regions(std::string ref_path,
                     std::string const & hla_vcf,
                     std::string const & interval_fn,
                     std::vector<std::string> const & sams,
                     std::vector<std::string> const & sam_index_paths,
                     std::vector<double> const & avg_cov_by_readlen,
                     std::vector<GenomicRegion> const & regions,
                     std::string const & output_path,
                     bool const is_copy_reference);

} // namespace gyper
