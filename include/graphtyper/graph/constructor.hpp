#pragma once

#include <string>
#include <vector>

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/graph.hpp>

namespace gyper
{

class GenomicRegion;
class Variant;

void construct_graph(std::string const & reference_filename,
                     std::string const & vcf_filename,
                     std::string const & region,
                     bool is_sv_graph = false,
                     bool use_absolute_positions = true,
                     bool check_index = true);

void
open_and_read_reference_genome(std::vector<char> & reference_sequence,
                               std::string const & reference_fn,
                               GenomicRegion const & genomic_region);

std::vector<Variant>
get_variants_using_tabix(std::string const & vcf, GenomicRegion const & genomic_region);

} // namespace gyper
