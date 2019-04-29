#pragma once
#include <vector>

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/graph.hpp>
#include <graphtyper/index/mem_index.hpp>
#include <graphtyper/typer/genotype_paths.hpp>
#include <graphtyper/utilities/sam_reader.hpp>

#include <seqan/sequence.h>
#include <seqan/bam_io.h>


namespace gyper
{

void
align_unpaired_read_pairs(TReads & reads, std::vector<GenotypePaths> & genos);


GenotypePaths find_genotype_paths_of_a_single_sequence(seqan::IupacString const & read, seqan::CharString const & qual, int const mismatches = -1, gyper::Graph const & graph = gyper::graph);

std::vector<std::pair<GenotypePaths, GenotypePaths> >
align_paired_reads(std::vector<TReadPair> const & records);

std::vector<GenotypePaths>
find_haplotype_paths(std::vector<seqan::Dna5String> const & sequences);

} // namespace gyper
