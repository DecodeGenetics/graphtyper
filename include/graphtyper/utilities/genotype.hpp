#pragma once

#include <string>
#include <vector>

#include <graphtyper/graph/genomic_region.hpp>

namespace bamshrink
{

struct Options;

}


namespace gyper
{

std::vector<std::string>
run_bamshrink(std::vector<std::string> const & sams,
              std::vector<std::string> const & sams_index,
              std::string const & ref_fn,
              GenomicRegion const & region,
              std::vector<double> const & avg_cov_by_readlen,
              std::string const & tmp);

// bamshrink variant for multiple regions
std::vector<std::string>
run_bamshrink_multi(std::vector<std::string> const & sams,
                    std::string const & ref_fn,
                    std::string const & interval_fn,
                    std::vector<double> const & avg_cov_by_readlen,
                    std::string const & tmp);

void
run_samtools_merge(std::vector<std::string> & shrinked_sams, std::string const & tmp);

void
genotype(std::string ref_path,
         std::vector<std::string> const & sams,
         std::vector<std::string> const & sams_index,
         GenomicRegion const & region,
         std::string const & output_path,
         std::vector<double> const & avg_cov_by_readlen,
         bool const is_copy_reference);


void
genotype_regions(std::string const & ref_path,
                 std::vector<std::string> const & sams,
                 std::vector<std::string> const & sams_index,
                 std::vector<GenomicRegion> const & regions,
                 std::string const & output_path,
                 std::vector<double> const & avg_cov_by_readlen,
                 bool const is_copy_reference);

} // namespace gyper
