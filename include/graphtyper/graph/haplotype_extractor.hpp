#pragma once

#include <string>
#include <vector>

#include <graphtyper/graph/genotype.hpp> // Genotype
#include <graphtyper/typer/variant.hpp> // Variant

#include <seqan/sequence.h> // Dna5String


namespace gyper
{

void extract_to_vcf(std::string graph_path, std::string haps_path, std::string output_folder, std::string const & region);

std::vector<VariantCandidate>
find_variants_in_alignment(uint32_t pos,
                           std::vector<char> const & ref,
                           seqan::Dna5String const & read
                           );

}
