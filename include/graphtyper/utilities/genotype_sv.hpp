#pragma once

#include <string>
#include <vector>

#include <graphtyper/graph/genomic_region.hpp>


namespace gyper
{

void
genotype_sv(std::string const & interval_fn,
            std::string const & ref_path,
            std::vector<std::string> const & sams,
            gyper::GenomicRegion const & genomic_region,
            std::string const & output_path,
            std::vector<double> const & avg_cov_by_readlen);

void
genotype_sv_regions(std::string ref_path,
                    std::string const & sv_vcf,
                    std::vector<std::string> const & sams,
                    std::vector<GenomicRegion> const & regions,
                    std::string const & output_path,
                    bool const is_copy_reference);

} // namespace gyper
