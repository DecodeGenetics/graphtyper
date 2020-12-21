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
                                            gyper::PHIndex const & ph_index,
                                            gyper::Graph const & graph)
{
  using namespace gyper;

  TKmerLabels r_hamming0 = query_index(read, ph_index);
  TKmerLabels r_hamming1 = query_index_hamming_distance1_without_index(read, ph_index);
  assert(r_hamming0.size() > 0);

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
    uint32_t read_start_index{0};

    {
      assert(r_hamming0.size() == r_hamming1.size());

      for (long i = 0; i < static_cast<long>(r_hamming0.size()); ++i)
      {
        geno.add_next_kmer_labels(r_hamming0[i],
                                  read_start_index,
                                  read_start_index + (K - 1),
                                  0 /*mismatches*/);

        geno.add_next_kmer_labels(r_hamming1[i],
                                  read_start_index,
                                  read_start_index + (K - 1),
                                  1 /*mismatches*/);

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

  if (graph.is_sv_graph)
    geno.remove_fully_special_paths();

  geno.remove_non_ref_paths_when_read_matches_ref(); // Should be the last check
  geno.update_longest_path_size();
  geno.remove_short_paths();

  if (graph.is_sv_graph)
    geno.remove_support_from_read_ends();

  // store read sequence
  geno.read2 = std::vector<char>(seqan::begin(read), seqan::end(read));

  // store score diff


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


long
clipped_count(bam1_t const & b)
{
  long count{0};

  if (b.core.n_cigar > 0)
  {
    auto it = b.data + b.core.l_qname;

    // Check first
    uint32_t opAndCnt;
    memcpy(&opAndCnt, it, sizeof(uint32_t));

    if ((opAndCnt & 15) == 4)
    {
      uint32_t const cigar_count = opAndCnt >> 4;
      count += cigar_count;
      //std::cerr << "MIDNSHP=X*******"[opAndCnt & 15] << " first ";
      return true;
    }

    // Check last
    memcpy(&opAndCnt, it + sizeof(uint32_t) * (b.core.n_cigar - 1), sizeof(uint32_t));

    if ((opAndCnt & 15) == 4)
    {
      //std::cerr << "MIDNSHP=X*******"[opAndCnt & 15] << " last ";
      uint32_t const cigar_count = opAndCnt >> 4;
      count += cigar_count;
      return true;
    }
  }

  return count;
}


uint8_t
get_score_diff(bam1_t * rec)
{
  uint8_t * it = bam_get_aux(rec);
  auto const l_aux = bam_get_l_aux(rec);
  unsigned i{0};
  int64_t as{-1};
  int64_t xs{-1};

  while (i < l_aux)
  {
    i += 3;
    char type = *(it + i - 1);
    bool is_as = false;
    bool is_xs = false;

    // Check for AS and XS
    if (*(it + i - 2) == 'S')
    {
      if (*(it + i - 3) == 'A')
        is_as = true;
      else if (*(it + i - 3) == 'X')
        is_xs = true;
    }

    switch (type)
    {
    case 'A':
    {
      // A printable character
      ++i;
      break;
    }

    case 'Z':
    {
      // A string! Let's loop it until qNULL
      // Check for RG
      while (*(it + i) != '\0' && *(it + i) != '\n')
        ++i;

      ++i;
      break;
    }

    case 'c':
    {
      if (is_as)
      {
        int8_t num;
        memcpy(&num, it + i, sizeof(int8_t));
        as = num;
      }
      else if (is_xs)
      {
        int8_t num;
        memcpy(&num, it + i, sizeof(int8_t));
        xs = num;
      }

      i += sizeof(int8_t);
      break;
    }

    case 'C':
    {
      if (is_as)
      {
        uint8_t num;
        memcpy(&num, it + i, sizeof(uint8_t));
        as = num;
      }
      else if (is_xs)
      {
        uint8_t num;
        memcpy(&num, it + i, sizeof(uint8_t));
        xs = num;
      }

      i += sizeof(uint8_t);
      break;
    }

    case 's':
    {
      if (is_as)
      {
        int16_t num;
        memcpy(&num, it + i, sizeof(int16_t));
        as = num;
      }
      else if (is_xs)
      {
        int16_t num;
        memcpy(&num, it + i, sizeof(int16_t));
        xs = num;
      }

      i += sizeof(int16_t);
      break;
    }

    case 'S':
    {
      if (is_as)
      {
        uint16_t num;
        memcpy(&num, it + i, sizeof(uint16_t));
        as = num;
      }
      else if (is_xs)
      {
        uint16_t num;
        memcpy(&num, it + i, sizeof(uint16_t));
        xs = num;
      }

      i += sizeof(uint16_t);
      break;
    }

    case 'i':
    {
      if (is_as)
      {
        int32_t num;
        memcpy(&num, it + i, sizeof(int32_t));
        as = num;
      }
      else if (is_xs)
      {
        int32_t num;
        memcpy(&num, it + i, sizeof(int32_t));
        xs = num;
      }

      i += sizeof(int32_t);
      break;
    }

    case 'I':
    {
      if (is_as)
      {
        uint32_t num;
        memcpy(&num, it + i, sizeof(uint32_t));
        as = num;
      }
      else if (is_xs)
      {
        uint32_t num;
        memcpy(&num, it + i, sizeof(uint32_t));
        xs = num;
      }

      i += sizeof(uint32_t);
      break;
    }

    case 'f':
    {
      i += sizeof(float);
      break;
    }

    default:
    {
      i = l_aux;       // Unkown tag, stop
      break;
    }
    }
  }


  if (as == -1 || as < xs)
  {
    return 0;
  }
  else
  {
    if (xs == -1)
      xs = 0;

    // check for overflow
    long const diff = as - xs;
    return diff < std::numeric_limits<uint8_t>::max() ? diff : std::numeric_limits<uint8_t>::max();
  }
}


} // anon namespace


namespace gyper
{


std::pair<GenotypePaths, GenotypePaths>
align_read(bam1_t * rec,
           seqan::IupacString const & seq,
           seqan::IupacString const & rseq,
           gyper::PHIndex const & ph_index)
{
  auto const & core = rec->core;

  std::pair<GenotypePaths, GenotypePaths> geno_paths =
    std::make_pair<GenotypePaths, GenotypePaths>(
      GenotypePaths(core.flag, core.l_qseq),
      GenotypePaths(core.flag, core.l_qseq));

  // Hard restriction on read length is 63 bp (2*32 - 1)
  if (seqan::length(seq) < (2 * K - 1))
    return geno_paths;

  if ((core.flag & IS_PAIRED) == 0u ||
      (core.tid == core.mtid &&
       core.isize > -1200 &&
       core.isize < 1200 &&
       (((core.flag & IS_SEQ_REVERSED) != 0u) != ((core.flag & IS_MATE_SEQ_REVERSED) != 0u))))
  {
    find_genotype_paths_of_one_of_the_sequences(seq, geno_paths.first, ph_index, graph);

    // skip reverse orientation graph alignment if the read fulfills the above.. unless
    if (Options::const_instance()->force_align_both_orientations)
      find_genotype_paths_of_one_of_the_sequences(rseq, geno_paths.second, ph_index, graph);
  }
  else
  {
    find_genotype_paths_of_one_of_the_sequences(seq, geno_paths.first, ph_index, graph);
    find_genotype_paths_of_one_of_the_sequences(rseq, geno_paths.second, ph_index, graph);
  }

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
    GenotypePaths & geno = geno_paths.first;
    geno.flags = core.flag & ~IS_PROPER_PAIR; // clear proper pair bit
    geno.mapq = core.qual;

    if (!(core.flag & IS_UNMAPPED)) // True if read was mapped
      geno.original_pos = core.pos;

    if (core.qual < 25) // True if MQ is below 25
      set_bit(geno.flags, IS_MAPQ_BAD);

    {
      long const clipped_count_bp = clipped_count(*rec);

      if (clipped_count_bp > 3)
        set_bit(geno.flags, IS_CLIPPED);
    }

    geno.score_diff = get_score_diff(rec);

#ifndef NDEBUG
    if (Options::instance()->stats.size() > 0)
    {
      geno.details->query_name = std::string(reinterpret_cast<char *>(rec->data));
      geno.qual2.resize(core.l_qseq);
      auto it = bam_get_qual(rec);

      for (int i = 0; i < core.l_qseq; ++it, ++i)
        geno.qual2[i] = static_cast<char>(*it + 33);
    }
#endif // NDEBUG

    return &geno;
  }

  case 2:
  {
    GenotypePaths & geno = geno_paths.second;
    geno.flags = (core.flag ^ IS_SEQ_REVERSED) & ~IS_PROPER_PAIR; // clear proper pair bit
    geno.mapq = core.qual;

    if (!(core.flag & IS_UNMAPPED)) // True if read was mapped
      geno.original_pos = core.pos;

    if (core.qual < 25) // True if MQ is below 25
      set_bit(geno.flags, IS_MAPQ_BAD);

    {
      long const clipped_count_bp = clipped_count(*rec);

      if (clipped_count_bp > 3)
        set_bit(geno.flags, IS_CLIPPED);
    }

    geno.score_diff = get_score_diff(rec);

#ifndef NDEBUG
    if (Options::instance()->stats.size() > 0)
    {
      geno.details->query_name = std::string(reinterpret_cast<char *>(rec->data));
      geno.qual2.resize(core.l_qseq);
      auto it = bam_get_qual(rec);

      for (int i = (core.l_qseq - 1); i >= 0; ++it, --i)
        geno.qual2[i] = static_cast<char>(*it + 33);
    }
#endif // NDEBUG

    return &geno;
  }

  default: break;
  }

  return nullptr;
}


void
further_update_unpaired_read_paths_for_discovery(GenotypePaths & geno, bam1_t * rec)
{
  assert(rec);
  auto const & core = rec->core;

  if ((geno.flags & IS_SEQ_REVERSED) == (core.flag & IS_SEQ_REVERSED))
  {
    // seq reversed bit has not been flipped
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


void
update_paths(std::pair<GenotypePaths, GenotypePaths> & geno_paths, bam1_t * rec)
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

  {
    long const clipped_count_bp = clipped_count(*rec);

    if (clipped_count_bp > 3)
    {
      set_bit(geno1.flags, IS_CLIPPED);
      set_bit(geno2.flags, IS_CLIPPED);
    }
  }

  {
    auto const score_diff = get_score_diff(rec);
    geno1.score_diff = score_diff;
    geno2.score_diff = score_diff;
  }

  geno2.flags = (core.flag ^ IS_SEQ_REVERSED) & ~IS_PROPER_PAIR;
  geno2.mapq = geno1.mapq;
  geno2.ml_insert_size = geno1.ml_insert_size;

#ifndef NDEBUG
  if (Options::instance()->stats.size() > 0)
  {
    geno1.details->query_name = std::string(reinterpret_cast<char *>(rec->data));
    geno2.details->query_name = geno1.details->query_name;

    geno1.qual2.resize(core.l_qseq);
    auto it = bam_get_qual(rec);

    for (int i = 0; i < core.l_qseq; ++it, ++i)
      geno1.qual2[i] = static_cast<char>(*it + 33);

    geno2.qual2 = std::vector<char>(geno1.qual2.rbegin(), geno1.qual2.rend());
  }
#endif // NDEBUG
}


void
further_update_paths_for_discovery(std::pair<GenotypePaths, GenotypePaths> & geno_paths, bam1_t * rec)
{
  assert(rec);

  auto const & core = rec->core;
  GenotypePaths & geno1 = geno_paths.first;
  GenotypePaths & geno2 = geno_paths.second;

  geno1.qual2.resize(core.l_qseq);
  auto it = bam_get_qual(rec);

  for (int i = 0; i < core.l_qseq; ++it, ++i)
    geno1.qual2[i] = static_cast<char>(*it + 33);

  geno2.qual2 = std::vector<char>(geno1.qual2.rbegin(), geno1.qual2.rend());
}


std::pair<GenotypePaths *, GenotypePaths *>
get_better_paths(std::pair<GenotypePaths, GenotypePaths> & geno_paths1,
                 std::pair<GenotypePaths, GenotypePaths> & geno_paths2)
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
