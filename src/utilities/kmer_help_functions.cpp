#include <iostream>
#include <vector>

#include <graphtyper/index/ph_index.hpp>
#include <graphtyper/utilities/kmer_help_functions.hpp>
#include <graphtyper/utilities/type_conversions.hpp>


namespace gyper
{

template <typename TSequence>
std::size_t
get_num_kmers(TSequence const & dna)
{
  if (seqan::length(dna) < K)
    return 0;
  else
    return 1 + (seqan::length(dna) - K) / (K - 1);
}


template <typename TSequence>
TSequence
get_ith_kmer(TSequence const & dna, std::size_t i)
{
  assert(seqan::length(dna) >= K);
  // assert (std::count_if(seqan::begin(dna), seqan::end(dna), seqan::Dna5('N')) == 0);
  TSequence new_dna(dna);
  seqan::erase(new_dna, 0, (seqan::length(dna) - K) % (K - 1) / 2 + (K - 1) * i);
  assert(seqan::length(new_dna) >= K);
  seqan::resize(new_dna, K);
  return new_dna;
}


// Explicit intantations
template std::size_t get_num_kmers(seqan::Dna5String const &);
template std::size_t get_num_kmers(seqan::IupacString const &);
template seqan::Dna5String get_ith_kmer(seqan::Dna5String const &, std::size_t);
template seqan::IupacString get_ith_kmer(seqan::IupacString const &, std::size_t);


template <typename TSeq>
std::vector<KmerLabel>
query_index_for_first_kmer(TSeq const & read, PHIndex const & ph_index)
{
  std::vector<uint64_t> keys = to_uint64_vec(read, 0);
  return ph_index.get(keys);
}


template <typename TSeq>
std::vector<KmerLabel>
query_index_for_last_kmer(TSeq const & read, PHIndex const & ph_index)
{
  std::vector<uint64_t> keys = to_uint64_vec(read, seqan::length(read) - K);
  return ph_index.get(keys);
}


template <typename TSeq>
std::vector<std::vector<KmerLabel> >
query_index(TSeq const & read, PHIndex const & ph_index)
{
  std::vector<std::vector<uint64_t> > multi_keys;
  long const num_keys = get_num_kmers(read);

  for (long i = 0; i < num_keys; ++i)
    multi_keys.push_back(to_uint64_vec(read, (K - 1) * i));

  return ph_index.multi_get(multi_keys);
}


// Explicit instantation
template std::vector<KmerLabel> query_index_for_first_kmer(seqan::IupacString const & read,
                                                           PHIndex const & ph_index);
template std::vector<KmerLabel> query_index_for_last_kmer(seqan::IupacString const & read, PHIndex const & ph_index);
template std::vector<std::vector<KmerLabel> > query_index<seqan::Dna5String>(seqan::Dna5String const &,
                                                                             PHIndex const & ph_index);
template std::vector<std::vector<KmerLabel> > query_index<seqan::IupacString>(seqan::IupacString const &,
                                                                              PHIndex const & ph_index);


/*
template <typename TSeq>
std::vector<std::vector<KmerLabel> >
query_index_hamming_distance1(TSeq const & read, gyper::PHIndex const & ph_index)
{
  std::vector<std::vector<uint64_t> > multi_keys;
  std::size_t const num_keys = get_num_kmers(read);

  for (unsigned i = 0; i < num_keys; ++i)
    multi_keys.push_back(to_uint64_vec(read, (K - 1) * i));

  return ph_index.multi_get_hamming1(multi_keys);
}


// Explicit instantation
template std::vector<std::vector<KmerLabel> >
query_index_hamming_distance1<seqan::Dna5String>(seqan::Dna5String const &, gyper::PHIndex const &);
template std::vector<std::vector<KmerLabel> >
query_index_hamming_distance1<seqan::IupacString>(seqan::IupacString const &, gyper::PHIndex const &);
*/

template <typename TSeq>
std::vector<std::vector<KmerLabel> >
query_index_hamming_distance1_without_index(TSeq const & read, gyper::PHIndex const & ph_index)
{
  std::vector<std::vector<uint64_t> > multi_keys;
  std::size_t const num_keys = get_num_kmers(read);

  // Find exact keys
  for (std::size_t i = 0; i < num_keys; ++i)
    multi_keys.push_back(to_uint64_vec(read, (K - 1) * i));

  // Find keys in hamming distance 1 to the exact keys
  for (std::size_t i = 0; i < num_keys; ++i)
  {
    // If the key is not unique, do not add any keys with hamming distance 1
    if (multi_keys[i].size() != 1)
      continue;

    multi_keys[i].reserve(96);
    std::array<uint64_t, 96> ham1_keys = to_uint64_vec_hamming_distance_1(multi_keys[i][0]);
    multi_keys[i][0] = ham1_keys[0]; // Remove the exact key from these results
    std::move(ham1_keys.begin() + 1, ham1_keys.end(), std::back_inserter(multi_keys[i]));
    assert(multi_keys[i].size() == 96);
  }

  return ph_index.multi_get(multi_keys);
}


// Explicit instantation
template std::vector<std::vector<KmerLabel> >
query_index_hamming_distance1_without_index<seqan::Dna5String>(seqan::Dna5String const &, gyper::PHIndex const &);
template std::vector<std::vector<KmerLabel> >
query_index_hamming_distance1_without_index<seqan::IupacString>(seqan::IupacString const &, gyper::PHIndex const &);


} // namespace gyper
