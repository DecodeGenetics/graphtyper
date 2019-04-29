#pragma once

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <graphtyper/graph/graph.hpp>
#include <graphtyper/index/kmer_label.hpp>
#include <graphtyper/index/rocksdb.hpp>
#include <graphtyper/index/mem_index.hpp>
#include <graphtyper/typer/genotype_paths.hpp>
#include <graphtyper/utilities/type_conversions.hpp>

namespace gyper
{

template <typename TSequence>
std::size_t get_num_kmers(TSequence const & dna);
template <typename TSequence>
TSequence get_ith_kmer(TSequence const & dna, std::size_t i);

template <typename TSequence>
uint32_t read_offset(TSequence const & dna);

template <typename TSeq>
std::vector<KmerLabel>
query_index_for_first_kmer(TSeq const & read, MemIndex const & _mem_index = gyper::mem_index);

template <typename TSeq>
std::vector<KmerLabel>
query_index_for_last_kmer(TSeq const & read, MemIndex const & _mem_index = gyper::mem_index);

template <typename TSeq>
std::vector<std::vector<KmerLabel> >
query_index(TSeq const & read, gyper::MemIndex const & mem_index = gyper::mem_index);

template <typename TSeq>
std::vector<std::vector<KmerLabel> >
query_index_hamming_distance1(TSeq const & read, gyper::MemIndex const & mem_index = gyper::mem_index);

template <typename TSeq>
std::vector<std::vector<KmerLabel> >
query_index_hamming_distance1_without_index(TSeq const & read, gyper::MemIndex const & mem_index = gyper::mem_index);

} // namespace gyper
