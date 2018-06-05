#pragma once

#include <string>
#include <vector>

#include <graphtyper/graph/genotype.hpp> // Genotype
#include <graphtyper/typer/variant.hpp> // Variant

#include <seqan/sequence.h> // Dna5String


namespace gyper
{

void extract_to_vcf(std::string const & graph_path,
                    std::string const & haps_path,
                    std::string const & output_folder,
                    std::string const & region
  );

Variant
make_variant_of_gapped_strings(std::string & gapped_ref,
                               std::string & gapped_alt,
                               long pos,
                               long & ref_to_seq_offset
  );

std::vector<VariantCandidate>
find_variants_in_alignment(uint32_t pos,
                           std::vector<char> const & ref,
                           seqan::Dna5String const & read,
                           std::vector<char> const & qual
                           );

}
