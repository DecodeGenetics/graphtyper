#include <algorithm>
#include <array>
#include <chrono>
#include <ctime>
#include <iterator>

#include <boost/log/trivial.hpp>

#include <htslib/sam.h>

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/graph_serialization.hpp>
#include <graphtyper/typer/alignment.hpp>
#include <graphtyper/utilities/kmer_help_functions.hpp>
#include <graphtyper/utilities/io.hpp>
#include <graphtyper/utilities/options.hpp>
#include <graphtyper/utilities/type_conversions.hpp>


namespace
{

void
find_genotype_paths_of_one_of_the_sequences(seqan::IupacString const & read,
                                            gyper::GenotypePaths & geno,
                                            gyper::Graph const & graph = gyper::graph,
                                            gyper::MemIndex const & mem_index = gyper::mem_index
                                            )
{
  using namespace gyper;

  TKmerLabels r_hamming0 = query_index(read, mem_index);
  TKmerLabels r_hamming1 = query_index_hamming_distance1_without_index(read, mem_index);

  // Stop if all kmer are extremely common
  for (auto it = r_hamming0.cbegin();;)
  {
    if (it->size() < MAX_UNIQUE_KMER_POSITIONS)
    {
      break;
    }
    else
    {
      ++it;

      // We found no k-mers with less than MAX_UNIQUE_KMER_POSITIONS locations!
      if (it == r_hamming0.cend())
        return;
    }
  }

  {
    uint32_t read_start_index = 0;

    {
      assert(r_hamming0.size() == r_hamming1.size());

      for (long i = 0; i < static_cast<long>(r_hamming0.size()); ++i)
      {
        geno.add_next_kmer_labels(r_hamming0[i],
                                  read_start_index,
                                  read_start_index + (K - 1),
                                  0 /*mismatches*/
                                  );

        geno.add_next_kmer_labels(r_hamming1[i],
                                  read_start_index,
                                  read_start_index + (K - 1),
                                  1 /*mismatches*/
                                  );

        read_start_index += (K - 1);
      }
    }
  }

  geno.remove_short_paths();

  // Extend the paths
  geno.walk_read_starts(read, -1, graph);
  geno.walk_read_ends(read, -1, graph);
  geno.update_longest_path_size();
  geno.remove_short_paths();

  // Filter the bad stuff
  geno.remove_paths_with_too_many_mismatches();
  geno.remove_fully_special_paths();
  geno.remove_non_ref_paths_when_read_matches_ref(); // Should be the last check

  geno.update_longest_path_size();
  geno.remove_short_paths();

  if (graph.is_sv_graph)
    geno.remove_support_from_read_ends();
  //else
  //  geno.force_full_overlap();

  // Can fail in new SV indel alignment :'(
#ifndef NDEBUG
  /*
  if (!geno.check_no_variant_is_missing())
  {
    std::cerr << "ERROR: Variant missing in read:\n";
    std::cerr << std::string(geno.read2.begin(), geno.read2.end()) << std::endl;
    assert(false);
  } //*/
#endif // NDEBUG
}


bool
is_clipped(bam1_t const & b)
{
  if (b.core.n_cigar == 0)
    return false;

  auto it = b.data + b.core.l_qname;

  // Check first
  uint32_t opAndCnt;
  memcpy(&opAndCnt, it, sizeof(uint32_t));

  if ((opAndCnt & 15) == 4)
  {
    //std::cerr << "MIDNSHP=X*******"[opAndCnt & 15] << " first ";
    return true;
  }

  // Check last
  memcpy(&opAndCnt, it + sizeof(uint32_t) * (b.core.n_cigar - 1), sizeof(uint32_t));

  if ((opAndCnt & 15) == 4)
  {
    //std::cerr << "MIDNSHP=X*******"[opAndCnt & 15] << " last ";
    return true;
  }

  return false;
}


} // anon namespace


namespace gyper
{


std::pair<GenotypePaths, GenotypePaths>
align_read(bam1_t * rec, seqan::IupacString const & seq, seqan::IupacString const & rseq)
{
  auto const & core = rec->core;

  std::pair<GenotypePaths, GenotypePaths> geno_paths = std::make_pair<GenotypePaths, GenotypePaths>(
    GenotypePaths(core.flag, core.l_qseq),
    GenotypePaths(core.flag, core.l_qseq));

  find_genotype_paths_of_one_of_the_sequences(seq, geno_paths.first);
  find_genotype_paths_of_one_of_the_sequences(rseq, geno_paths.second);
  return geno_paths;
}


GenotypePaths *
update_unpaired_read_paths(std::pair<GenotypePaths, GenotypePaths> & geno_paths, bam1_t * rec)
{
  assert(rec);

  auto const & core = rec->core;

  switch (compare_pair_of_genotype_paths(geno_paths.first, geno_paths.second))
  {
  case 1:
  {
    GenotypePaths * geno = &geno_paths.first;
    geno->flags = core.flag & ~IS_PROPER_PAIR; // clear proper pair bit
    geno->mapq = core.qual;

    if (!(core.flag & IS_UNMAPPED)) // True if read was mapped
      geno->original_pos = core.pos;

    if (core.qual < 25) // True if MQ is below 25
      set_bit(geno->flags, IS_MAPQ_BAD);

    if (is_clipped(*rec))
      set_bit(geno->flags, IS_CLIPPED);

#ifndef NDEBUG
    if (Options::instance()->stats.size() > 0)
      geno->details->query_name = std::string(reinterpret_cast<char *>(rec->data));
#endif // NDEBUG

    return geno;
  }

  case 2:
  {
    GenotypePaths * geno = &geno_paths.second;
    geno->flags = (core.flag ^ IS_SEQ_REVERSED) & ~IS_PROPER_PAIR; // clear proper pair bit
    geno->mapq = core.qual;

    if (!(core.flag & IS_UNMAPPED)) // True if read was mapped
      geno->original_pos = core.pos;

    if (core.qual < 25) // True if MQ is below 25
      set_bit(geno->flags, IS_MAPQ_BAD);

    if (is_clipped(*rec))
      set_bit(geno->flags, IS_CLIPPED);

#ifndef NDEBUG
    if (Options::instance()->stats.size() > 0)
      geno->details->query_name = std::string(reinterpret_cast<char *>(rec->data));
#endif // NDEBUG

    return geno;
  }

  default: break;
  }

  return nullptr;
}


void
further_update_unpaired_read_paths_for_discovery(GenotypePaths & geno,
                                                 seqan::IupacString const & seq,
                                                 seqan::IupacString const & rseq,
                                                 bam1_t * rec
                                                 )
{
  assert(rec);
  auto const & core = rec->core;

  if ((geno.flags & IS_SEQ_REVERSED) == (core.flag & IS_SEQ_REVERSED))
  {
    // seq reversed bit has not been flipped
    geno.read2 = std::vector<char>(seqan::begin(seq), seqan::end(seq));
    geno.qual2.resize(core.l_qseq);
    auto it = bam_get_qual(rec);

    for (int i = 0; i < core.l_qseq; ++it, ++i)
    {
      assert(i < static_cast<long>(geno.qual2.size()));
      geno.qual2[i] = static_cast<char>(*it + 33);
    }
  }
  else
  {
    geno.read2 = std::vector<char>(seqan::begin(rseq), seqan::end(rseq));
    geno.qual2.resize(core.l_qseq);
    auto it = bam_get_qual(rec);

    for (int i = 0; i < core.l_qseq; ++it, ++i)
    {
      assert(core.l_qseq - 1 - i >= 0l);
      assert(core.l_qseq - 1 - i < static_cast<long>(geno.qual2.size()));
      geno.qual2[core.l_qseq - 1 - i] = static_cast<char>(*it + 33);
    }
  }
}


#ifndef NDEBUG
void
update_paths(std::pair<GenotypePaths, GenotypePaths> & geno_paths,
             seqan::IupacString const & seq,
             seqan::IupacString const & rseq,
             bam1_t * rec)
#else
void
update_paths(std::pair<GenotypePaths, GenotypePaths> & geno_paths, bam1_t * rec)
#endif
{
  assert(rec);
  auto const & core = rec->core;

  GenotypePaths & geno1 = geno_paths.first;
  GenotypePaths & geno2 = geno_paths.second;

  geno1.flags = core.flag & ~IS_PROPER_PAIR; // clear proper pair bit
  geno1.mapq = core.qual;
  geno1.ml_insert_size = std::abs(core.isize);

  if (!(core.flag & IS_UNMAPPED))
  {
    geno1.original_pos = core.pos;
    geno2.original_pos = geno1.original_pos;
  }

  if (core.qual < 25)
    set_bit(geno1.flags, IS_MAPQ_BAD);

  if (is_clipped(*rec))
  {
    set_bit(geno1.flags, IS_CLIPPED);
    set_bit(geno2.flags, IS_CLIPPED);
  }

  geno2.flags = (core.flag ^ IS_SEQ_REVERSED) & ~IS_PROPER_PAIR;
  geno2.mapq = geno1.mapq;
  geno2.ml_insert_size = geno1.ml_insert_size;

#ifndef NDEBUG
  if (Options::instance()->stats.size() > 0)
  {
    geno1.details->query_name = std::string(reinterpret_cast<char *>(rec->data));
    geno2.details->query_name = geno1.details->query_name;

    geno1.read2 = std::vector<char>(seqan::begin(seq), seqan::end(seq));
    geno1.qual2.resize(core.l_qseq);
    auto it = bam_get_qual(rec);

    for (int i = 0; i < core.l_qseq; ++it, ++i)
      geno1.qual2[i] = static_cast<char>(*it + 33);

    geno2.read2 = std::vector<char>(seqan::begin(rseq), seqan::end(rseq));
    geno2.qual2 = std::vector<char>(geno1.qual2.rbegin(), geno1.qual2.rend());
  }
#endif // NDEBUG
}


void
further_update_paths_for_discovery(std::pair<GenotypePaths, GenotypePaths> & geno_paths,
                                   seqan::IupacString const & seq,
                                   seqan::IupacString const & rseq,
                                   bam1_t * rec
                                   )
{
  auto const & core = rec->core;
  GenotypePaths & geno1 = geno_paths.first;
  GenotypePaths & geno2 = geno_paths.second;

  assert(seqan::length(seq) == geno1.read_length);
  assert(seqan::length(seq) == geno2.read_length);

  geno1.read2 = std::vector<char>(seqan::begin(seq), seqan::end(seq));
  geno1.qual2.resize(core.l_qseq);
  auto it = bam_get_qual(rec);

  for (int i = 0; i < core.l_qseq; ++it, ++i)
    geno1.qual2[i] = static_cast<char>(*it + 33);

  geno2.read2 = std::vector<char>(seqan::begin(rseq), seqan::end(rseq));
  geno2.qual2 = std::vector<char>(geno1.qual2.rbegin(), geno1.qual2.rend());
}


std::pair<GenotypePaths *, GenotypePaths *>
get_better_paths(std::pair<GenotypePaths, GenotypePaths> & geno_paths1,
                 std::pair<GenotypePaths, GenotypePaths> & geno_paths2
                 )
{
  // contains first read w. forward strand and second read reverse strand
  std::pair<GenotypePaths *, GenotypePaths *> genos1;
  std::pair<GenotypePaths *, GenotypePaths *> genos2;

  {
    std::array<GenotypePaths *, 4> arr = {nullptr, nullptr, nullptr, nullptr};
    // index in array:
    //  0 is second in pair, reverse strand
    //  1 is first in pair, reverse strand
    //  2 is second in pair, forward strand
    //  3 is first in pair, forward strand
    auto get_index = [](uint16_t const flags) -> int
                     {
                       return ((flags & IS_FIRST_IN_PAIR) != 0) + 2 * ((flags & IS_SEQ_REVERSED) == 0);
                     };

    arr[get_index(geno_paths1.first.flags)] = &geno_paths1.first;
    arr[get_index(geno_paths1.second.flags)] = &geno_paths1.second;
    arr[get_index(geno_paths2.first.flags)] = &geno_paths2.first;
    arr[get_index(geno_paths2.second.flags)] = &geno_paths2.second;

    assert(arr[0]);
    assert(arr[1]);
    assert(arr[2]);
    assert(arr[3]);

    // Make sure we got a ptr in every field in release mode as well
    if (!arr[0] || !arr[1] || !arr[2] || !arr[3])
    {
      BOOST_LOG_TRIVIAL(warning) << "Unexpected read orientation, ptr array: "
                                 << arr[0] << " " << arr[1] << " " << arr[2] << " " << arr[3];
      genos1.first = nullptr;
      return genos1;
    }

    genos1 = {arr[3], arr[0]};
    genos2 = {arr[1], arr[2]};
  }

  switch (compare_pair_of_genotype_paths(genos1, genos2))
  {
  case 1:
    set_bit(genos1.first->flags, IS_PROPER_PAIR);
    set_bit(genos1.second->flags, IS_PROPER_PAIR);
    return genos1;

  case 2:
    set_bit(genos2.first->flags, IS_PROPER_PAIR);
    set_bit(genos2.second->flags, IS_PROPER_PAIR);
    return genos2;

  default:
    genos1.first = nullptr;
    return genos1;
  }

  return genos1;
}


} // namespace gyper
