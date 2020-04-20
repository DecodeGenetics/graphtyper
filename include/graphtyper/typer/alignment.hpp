#pragma once
#include <vector>

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/graph.hpp>
#include <graphtyper/index/ph_index.hpp>
#include <graphtyper/typer/genotype_paths.hpp>
#include <graphtyper/utilities/sam_reader.hpp>

#include <seqan/sequence.h>
#include <seqan/bam_io.h>


namespace gyper
{

std::pair<GenotypePaths, GenotypePaths>
align_read(bam1_t * rec,
           seqan::IupacString const & seq,
           seqan::IupacString const & rseq,
           gyper::PHIndex const & ph_index);

GenotypePaths *
update_unpaired_read_paths(std::pair<GenotypePaths, GenotypePaths> & geno_paths, bam1_t * rec);

void
further_update_unpaired_read_paths_for_discovery(GenotypePaths & geno,
                                                 seqan::IupacString const & seq,
                                                 seqan::IupacString const & rseq,
                                                 bam1_t * rec
                                                 );

#ifndef NDEBUG
void
update_paths(std::pair<GenotypePaths, GenotypePaths> & geno_paths,
             seqan::IupacString const & seq,
             seqan::IupacString const & rseq,
             bam1_t * rec);
#endif // NDEBUG

#ifdef NDEBUG
void
update_paths(std::pair<GenotypePaths, GenotypePaths> & geno_paths, bam1_t * rec);
#endif // DEBUG


void
further_update_paths_for_discovery(std::pair<GenotypePaths, GenotypePaths> & geno_paths,
                                   seqan::IupacString const & seq,
                                   seqan::IupacString const & rseq,
                                   bam1_t * rec
                                   );

std::pair<GenotypePaths *, GenotypePaths *>
get_better_paths(std::pair<GenotypePaths, GenotypePaths> & geno_paths1,
                 std::pair<GenotypePaths, GenotypePaths> & geno_paths2
                 );

} // namespace gyper
