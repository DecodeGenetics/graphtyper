#pragma once

#include <cstdint> // uint32_t
#include <string>  // std::string
#include <vector>  // std::vector

#include <graphtyper/graph/genotype.hpp> // Genotype
#include <graphtyper/typer/variant.hpp>  // Variant

namespace gyper
{
class Vcf;

void extract_to_vcf(std::string const & graph_path,
                    std::string const & haps_path,
                    std::string const & output_vcf,
                    std::string const & region,
                    bool const is_splitting_vars);

void extract_to_vcf(Vcf & haps_vcf,
                    std::vector<std::string> const & haps_paths,
                    std::string const & output_vcf,
                    bool const is_splitting_vars);

Variant make_variant_of_gapped_strings(std::string & gapped_ref,
                                       std::string & gapped_alt,
                                       long pos,
                                       long & ref_to_seq_offset);

std::vector<VariantCandidate> find_variants_in_alignment(uint32_t pos,
                                                         std::vector<char> const & ref,
                                                         std::vector<char> const & seq,
                                                         std::vector<char> const & qual);

} // namespace gyper
