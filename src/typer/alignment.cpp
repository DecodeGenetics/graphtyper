#include <array>
#include <chrono>
#include <ctime>
#include <mutex>

#include <boost/log/trivial.hpp>

#include <seqan/basic.h>
#include <seqan/modifier.h>
#include <seqan/stream.h>

#include <graphtyper/graph/graph_serialization.hpp>
#include <graphtyper/typer/alignment.hpp>
#include <graphtyper/utilities/kmer_help_functions.hpp>
#include <graphtyper/utilities/io.hpp>
#include <graphtyper/utilities/options.hpp>
#include <graphtyper/utilities/type_conversions.hpp>


namespace
{

bool
is_clipped(seqan::String<seqan::CigarElement<> > const & cigar)
{
  return seqan::length(cigar) == 0 ||
    seqan::front(cigar).operation == 'S' ||
    seqan::front(cigar).operation == 'H' ||
    seqan::back(cigar).operation == 'S' ||
    seqan::back(cigar).operation == 'H';
}

void
merge_index_queries(seqan::IupacString const & read,
                    gyper::GenotypePaths & geno,
                    gyper::TKmerLabels const & r_hamming0,
                    gyper::TKmerLabels const & r_hamming1,
                    gyper::Graph const & graph = gyper::graph
                    )
{
  using namespace gyper;

  {
    auto min_it = std::min_element(r_hamming0.begin(),
                                   r_hamming0.end(),
                                   [](std::vector<KmerLabel> const & l1, std::vector<KmerLabel> const & l2)
                                   {
                                     return l1.size() < l2.size();
                                   });

    assert (r_hamming0.size() > 0);
    assert(min_it != r_hamming0.end());

    if (min_it->size() > Options::instance()->MAX_UNIQUE_KMER_POSITIONS)
      return;

    uint32_t read_start_index = 0;

    if (r_hamming1.size() == 0)
    {
      for (unsigned i = 0; i < r_hamming0.size(); ++i)
      {
        geno.add_next_kmer_labels(r_hamming0[i], read_start_index, read_start_index + (K - 1), 0 /*mismatches*/);
        read_start_index += (K - 1);
      }
    }
    else
    {
      assert (r_hamming0.size() == r_hamming1.size());

      for (std::size_t i = 0; i < r_hamming0.size(); ++i)
      {
        geno.add_next_kmer_labels(r_hamming0[i], read_start_index, read_start_index + (K - 1), 0 /*mismatches*/);
        geno.add_next_kmer_labels(r_hamming1[i], read_start_index, read_start_index + (K - 1), 1 /*mismatches*/);
        read_start_index += (K - 1);
      }
    }
  }

  geno.remove_short_paths();
  geno.walk_read_starts(read, -1, graph);
  geno.walk_read_ends(read, -1, graph);
  geno.remove_short_paths();
  geno.remove_paths_with_too_many_mismatches();
  geno.remove_short_paths();

#ifndef NDEBUG
  if (!geno.check_no_variant_is_missing())
  {
    std::cout << geno.read << std::endl;
    assert(false);
  }
#endif // NDEBUG
}


void
find_genotype_paths_of_one_of_the_sequences(seqan::IupacString const & read, gyper::GenotypePaths & geno, bool const hamming_distance1_index_available, gyper::Graph const & graph = gyper::graph, gyper::MemIndex const & mem_index = gyper::mem_index) // By default index should be available
{
  using namespace gyper;

  TKmerLabels r_hamming0 = query_index(read, mem_index);
  TKmerLabels r_hamming1;

  if (Options::instance()->always_query_hamming_distance_one)
  {
    if (hamming_distance1_index_available)
      r_hamming1 = query_index_hamming_distance1(read);
    else
      r_hamming1 = query_index_hamming_distance1_without_index(read, mem_index);

    merge_index_queries(read, geno, r_hamming0, r_hamming1, graph);
  }
  else
  {
    merge_index_queries(read, geno, r_hamming0, r_hamming1, graph);

    // Check if a sufficiently good path was found, otherwise allow hamming distance one when quering index
    if (geno.longest_path_size() < ((3 * K) - 2))
    {
      if (hamming_distance1_index_available)
        r_hamming1 = query_index_hamming_distance1(read);
      else
        r_hamming1 = query_index_hamming_distance1_without_index(read, mem_index);

      assert (r_hamming0.size() == r_hamming1.size());
      merge_index_queries(read, geno, r_hamming0, r_hamming1, graph);
    }
  }
}


std::string
get_read_group(seqan::CharString tag_string)
{
  seqan::BamTagsDict tags_dict(tag_string);
  unsigned tagIdx = 0;
  std::string read_group("NA");

  if (seqan::findTagKey(tagIdx, tags_dict, "RG"))
  {
    seqan::CharString read_group_raw;

    if (seqan::extractTagValue(read_group_raw, tags_dict, tagIdx))
      read_group = std::string(seqan::toCString(read_group_raw));
  }

  return read_group;

  /*
  std::string tag_string(seqan::toCString(tag_seqan_string));
  std::size_t const n = tag_string.size();
  std::string const key = "RG:Z:";
  std::size_t i = 0;
  std::string read_group("NA");

  while (i != std::string::npos && i + key.size() < n)
  {
    std::cerr << "reading tag: " << std::string(tag_string.begin() + i, tag_string.begin() + i + n) << "\n";

    if (std::equal(key.begin(), key.end(), tag_string.begin() + i))
    {
      auto find_it = std::find(tag_string.begin() + i + key.size(), tag_string.end(), '\t');
      return std::string(tag_string.begin() + i, find_it);
    }

    i = tag_string.find('\t', i + 1);
  }

  return read_group;
  */
}


} // anon namespace


namespace gyper
{

GenotypePaths
align_a_single_sequence_without_hamming_distance1_index(seqan::BamAlignmentRecord const & record, Graph const & graph, MemIndex const & mem_index)
{
  GenotypePaths geno(record.seq, record.qual, record.mapQ);
  find_genotype_paths_of_one_of_the_sequences(record.seq, geno, false /*no hamming distance 1 index available*/, graph, mem_index);
  return geno;
}


void
align_unpaired_read_pairs(TReads & reads, std::vector<GenotypePaths> & genos)
{
  for (auto read_it = reads.begin(); read_it != reads.end(); ++read_it)
  {
    GenotypePaths geno1(read_it->first.seq, read_it->first.qual, read_it->first.mapQ);
    find_genotype_paths_of_one_of_the_sequences(read_it->first.seq, geno1, false /*No hamming1 distance index*/);
    //seqan::BamAlignmentRecord reverse(read_it->first);
    seqan::reverseComplement(read_it->first.seq);
    seqan::reverse(read_it->first.qual);
    GenotypePaths geno2(read_it->first.seq, read_it->first.qual, read_it->first.mapQ);
    find_genotype_paths_of_one_of_the_sequences(read_it->first.seq, geno2, false /*No hamming1 distance index*/);

    switch (compare_pair_of_genotype_paths(geno1, geno2))
    {
    case 1:
      // geno1.forward_strand is true by default
      geno1.is_first_in_pair = seqan::hasFlagFirst(read_it->first);
      geno1.is_originally_unaligned = seqan::hasFlagUnmapped(read_it->first);
      geno1.original_pos = read_it->first.beginPos + 1; // Change to 1-based system
      geno1.is_originally_clipped = seqan::hasFlagUnmapped(read_it->first) || is_clipped(read_it->first.cigar);

      if (Options::instance()->stats.size() > 0)
      {
        geno1.query_name = seqan::toCString(read_it->first.qName);
        geno1.details->query_name = seqan::toCString(read_it->first.qName);
        geno1.details->read_group = ::get_read_group(read_it->first.tags);
      }

      genos.push_back(std::move(geno1));
      break;

    case 2:
      geno2.forward_strand = false;
      geno2.is_first_in_pair = seqan::hasFlagFirst(read_it->first);
      geno2.is_originally_unaligned = seqan::hasFlagUnmapped(read_it->first);
      geno2.original_pos = read_it->first.beginPos + 1; // Change to 1-based system
      geno2.is_originally_clipped = seqan::hasFlagUnmapped(read_it->first) || is_clipped(read_it->first.cigar);

      if (Options::instance()->stats.size() > 0)
      {
        geno2.query_name = seqan::toCString(read_it->first.qName);
        geno2.details->query_name = seqan::toCString(read_it->first.qName);
        geno2.details->read_group = ::get_read_group(read_it->first.tags);
      }

      genos.push_back(std::move(geno2));
      break;
    }
  }
}


int64_t
get_insert_size(std::vector<Path>::const_iterator it1,
                std::vector<Path>::const_iterator it2,
                uint32_t const OPTIMAL,
                bool const REVERSE_COMPLEMENT
  )
{
  std::unordered_set<int64_t> distances;

  if (REVERSE_COMPLEMENT)
  {
    std::vector<Location> ll1 = graph.get_locations_of_a_position(it2->start, *it2);
    std::vector<Location> ll2 = graph.get_locations_of_a_position(it1->end, *it1);
    distances = graph.reference_distance_between_locations(ll1, ll2);
  }
  else
  {
    std::vector<Location> ll1 = graph.get_locations_of_a_position(it1->start, *it1);
    std::vector<Location> ll2 = graph.get_locations_of_a_position(it2->end, *it2);
    distances = graph.reference_distance_between_locations(ll1, ll2);
  }

  // Find the distance which is closest to the OPTIMAL
  int64_t best_distance = 0x00000000FFFFFFFFll;

  for (auto it = distances.begin(); it != distances.end(); ++it)
  {
    if (std::llabs(*it - OPTIMAL) < std::llabs(best_distance - OPTIMAL))
      best_distance = *it;
  }

  return best_distance;
}


int64_t
find_shortest_distance(GenotypePaths const & geno1,
                       GenotypePaths const & geno2,
                       uint32_t const OPTIMAL,
                       bool const REVERSE_COMPLEMENT
  )
{
  int64_t shortest_distance_diff = 0x00000000FFFFFFFFll;
  int64_t shortest_distance = 0x00000000FFFFFFFFll;

  for (auto it1 = geno1.paths.cbegin(); it1 != geno1.paths.cend(); ++it1)
  {
    for (auto it2 = geno2.paths.cbegin(); it2 != geno2.paths.cend(); ++it2)
    {
      int64_t distance = get_insert_size(it1, it2, OPTIMAL, REVERSE_COMPLEMENT);
      int64_t distance_diff = std::llabs(distance - static_cast<int64_t>(OPTIMAL));

      if (distance_diff < shortest_distance_diff)
      {
        shortest_distance_diff = distance_diff;
        shortest_distance = distance;
      }
    }
  }

  return shortest_distance;
}


void
remove_distant_paths(GenotypePaths & geno1,
                     GenotypePaths & geno2,
                     int64_t const SHORTEST_DISTANCE,
                     uint32_t const OPTIMAL,
                     bool const REVERSE_COMPLEMENT
  )
{
  {
    auto it1 = geno1.paths.begin();

    while (it1 != geno1.paths.end())
    {
      bool found_close_match = false;

      for (auto it2 = geno2.paths.cbegin(); it2 != geno2.paths.cend(); ++it2)
      {
        int64_t distance = get_insert_size(it1, it2, OPTIMAL, REVERSE_COMPLEMENT);
        distance = std::abs(distance - static_cast<int64_t>(OPTIMAL));

        if (distance <= SHORTEST_DISTANCE)
        {
          found_close_match = true;
          break;
        }
      }

      if (not found_close_match)
        it1 = geno1.paths.erase(it1); // Didn't find anyone close
      else
        ++it1;
    }
  }

  if (geno1.paths.size() == 0)
  {
    geno1.clear_paths(); // call clear_paths to reset the maximum path length
    geno2.clear_paths();
  }
  else
  {
    auto it2 = geno2.paths.begin();

    while (it2 != geno2.paths.end())
    {
      bool found_close_match = false;

      for (auto it1 = geno1.paths.cbegin(); it1 != geno1.paths.cend(); ++it1)
      {
        int64_t distance = get_insert_size(it1, it2, OPTIMAL, REVERSE_COMPLEMENT);
        distance = std::abs(distance - static_cast<int64_t>(OPTIMAL));

        if (distance <= SHORTEST_DISTANCE)
        {
          found_close_match = true;
          break;
        }
      }

      if (not found_close_match)
        it2 = geno2.paths.erase(it2); // Didn't find anyone close
      else
        ++it2;
    }
  }
}


GenotypePaths
find_genotype_paths_of_a_single_sequence(seqan::IupacString const & read, seqan::CharString const & qual, int const mismatches, gyper::Graph const & graph)
{
  uint32_t read_start_index = 0;
  TKmerLabels r1 = query_index(read);
  GenotypePaths geno(read, qual);

  for (unsigned i = 0; i < r1.size(); ++i)
  {
    geno.add_next_kmer_labels(r1[i], read_start_index, read_start_index + (K - 1), 0 /*mismatches*/);
    read_start_index += (K - 1);
  }

  // Compare read ends to the graph
  geno.walk_read_starts(read, mismatches /*max mismatches*/, graph);
  geno.walk_read_ends(read, mismatches /*max mismatches*/, graph);
  geno.walk_read_starts(read, mismatches /*max mismatches*/, graph);
  geno.remove_short_paths();
  return geno;
}


std::vector<GenotypePaths>
find_haplotype_paths(std::vector<seqan::Dna5String> const & sequences)
{
  std::vector<GenotypePaths> hap_paths;
  uint32_t count_too_short_sequences = 0;

  for (unsigned i = 0; i < sequences.size(); ++i)
  {
    if (seqan::length(sequences[i]) < 50)
    {
      // std::cout << "Skipped a sequence" << std::endl;
      GenotypePaths new_geno;
      new_geno.longest_path_length = 0;
      hap_paths.push_back(std::move(new_geno));
      continue;
    }

    GenotypePaths new_geno = find_genotype_paths_of_a_single_sequence(sequences[i], "" /*qual*/, 0); // We allow no mismatches

    // Merge the two genotype paths
    {
      if (new_geno.longest_path_length != seqan::length(sequences[i])) // Everything must align
      {
        std::cout << i << " " << new_geno.longest_path_length << " vs. " << seqan::length(sequences[i]) << std::endl;
        std::cout << sequences[i] << std::endl;

        // if (new_geno.paths.size() > 0)
        // {
        //   for (auto const & path : new_geno.paths)
        //     std::cout << "new_geno @" << (path.start_pos() - 1061198324ul) << "-" << (path.end_pos() - 1061198324ul)
        //                << " starting from " << path.read_start_index << std::endl;
        // }

        new_geno.longest_path_length = 0;
        ++count_too_short_sequences;
      }

      hap_paths.push_back(std::move(new_geno));
    }
  }

  if (count_too_short_sequences > 0)
    BOOST_LOG_TRIVIAL(info) << "[graphtyper::alignment] Could not align " << count_too_short_sequences << " sequences.";

  return hap_paths;
}


bool
support_same_path(std::pair<GenotypePaths, GenotypePaths> const & genos)
{
  using ReadExplain = std::unordered_map<uint32_t, std::bitset<MAX_NUMBER_OF_HAPLOTYPES> >; // Variant order to explain
  ReadExplain read_explain1;
  ReadExplain read_explain2;

  for (auto const & path : genos.first.paths)
  {
    for (unsigned i = 0; i < path.var_order.size(); ++i)
      read_explain1[path.var_order[i]] |= path.nums[i];
  }

  for (auto const & path : genos.second.paths)
  {
    for (unsigned i = 0; i < path.var_order.size(); ++i)
      read_explain2[path.var_order[i]] |= path.nums[i];
  }

  for (auto it1 = read_explain1.begin(); it1 != read_explain1.end(); ++it1)
  {
    auto it2 = read_explain2.find(it1->first);

    if (it2 != read_explain2.end())
    {
      if ((it1->second & it2->second).none())
        return false;
    }
  }

  return true;
}


std::pair<GenotypePaths, GenotypePaths>
find_genotype_paths_of_a_sequence_pair(seqan::BamAlignmentRecord const & record1,
                                       seqan::BamAlignmentRecord const & record2,
                                       bool const REVERSE_COMPLEMENT
                                       )
{
  // Create two empty paths, one for each read
  std::pair<GenotypePaths, GenotypePaths> genos =
    std::make_pair(GenotypePaths(record1.seq, record1.qual, record1.mapQ),
                   GenotypePaths(record2.seq, record2.qual, record2.mapQ)
                   );

  find_genotype_paths_of_one_of_the_sequences(record1.seq, genos.first, false /*true when hamming distance 1 index is available*/);
  find_genotype_paths_of_one_of_the_sequences(record2.seq, genos.second, false /*true when hamming distance 1 index is available*/);

  // Remove distant paths (from optimal insert size)
  if (genos.first.paths.size() > 0 && genos.second.paths.size() > 0) // Both reads aligned
  {
    uint32_t const OPTIMAL = Options::instance()->optimal_insert_size;
    uint32_t const THRESHOLD = Options::instance()->max_insert_size_threshold;
    int64_t const INSERT_SIZE = find_shortest_distance(genos.first, genos.second, OPTIMAL, REVERSE_COMPLEMENT);
    int64_t const SHORTEST_DISTANCE = std::min(static_cast<int64_t>(Options::instance()->max_insert_size),
                                               static_cast<int64_t>(std::llabs(INSERT_SIZE - static_cast<int64_t>(OPTIMAL)))
                                               );

    if (REVERSE_COMPLEMENT)
    {
      genos.first.ml_insert_size = -INSERT_SIZE;
      genos.second.ml_insert_size = INSERT_SIZE;
    }
    else
    {
      genos.first.ml_insert_size = INSERT_SIZE;
      genos.second.ml_insert_size = -INSERT_SIZE;
    }

    // Store total mismatches in each read pair
    {
      uint8_t const READ_PAIR_MISMATCHES = genos.first.paths[0].mismatches + genos.second.paths[0].mismatches;
      genos.first.read_pair_mismatches = READ_PAIR_MISMATCHES;
      genos.second.read_pair_mismatches = READ_PAIR_MISMATCHES;
    }

    if (Options::instance()->hq_reads)
    {
      if (SHORTEST_DISTANCE > static_cast<int64_t>(THRESHOLD) || !support_same_path(genos))
      {
        genos.first.clear_paths();
        genos.second.clear_paths();
      }
      else if (genos.first.paths.size() > 1 || genos.second.paths.size() > 1)
      {
        int64_t HQ_THRESHOLD = std::min(static_cast<int64_t>(THRESHOLD), static_cast<int64_t>(SHORTEST_DISTANCE + 5ll));
        remove_distant_paths(genos.first, genos.second, HQ_THRESHOLD, OPTIMAL, REVERSE_COMPLEMENT);
        assert(genos.first.paths.size() > 0 && genos.second.paths.size() > 1);
      }
    }
    else
    {
      remove_distant_paths(genos.first, genos.second, SHORTEST_DISTANCE + THRESHOLD, OPTIMAL, REVERSE_COMPLEMENT);
    }
  }

  return genos;
}


std::vector<std::pair<GenotypePaths, GenotypePaths> >
get_best_genotype_paths(std::vector<TReadPair> const & records)
{
  std::vector<std::pair<GenotypePaths, GenotypePaths> > genos;

  for (auto record_it = records.cbegin(); record_it != records.cend(); ++record_it)
  {
    assert(std::distance(seqan::begin(record_it->first.seq), seqan::end(record_it->first.seq)) >= 2 * K - 1);
    assert(std::distance(seqan::begin(record_it->second.seq), seqan::end(record_it->second.seq)) >= 2 * K - 1);
    std::pair<GenotypePaths, GenotypePaths> genos1 =
      find_genotype_paths_of_a_sequence_pair(record_it->first, record_it->second, false /*REVERSE_COMPLEMENT*/);

    seqan::BamAlignmentRecord rec_first(record_it->first);
    seqan::BamAlignmentRecord rec_second(record_it->second);
    seqan::reverseComplement(rec_first.seq);
    seqan::reverse(rec_first.qual);
    seqan::reverseComplement(rec_second.seq);
    seqan::reverse(rec_second.qual);
    std::pair<GenotypePaths, GenotypePaths> genos2 =
      find_genotype_paths_of_a_sequence_pair(rec_first, rec_second, true /*REVERSE_COMPLEMENT*/);

    switch (compare_pair_of_genotype_paths(genos1, genos2))
    {
    case 1:
      genos1.second.forward_strand = false; // The second read in the pair has been reverse complemented
      genos1.first.is_originally_unaligned = seqan::hasFlagUnmapped(record_it->first);
      genos1.second.is_originally_unaligned = seqan::hasFlagUnmapped(record_it->second);
      genos1.first.is_originally_clipped = seqan::hasFlagUnmapped(record_it->first) || is_clipped(record_it->first.cigar);
      genos1.second.is_originally_clipped = seqan::hasFlagUnmapped(record_it->second) || is_clipped(record_it->second.cigar);
      genos1.second.is_first_in_pair = false;
      genos1.first.original_pos = record_it->first.beginPos + 1;
      genos1.second.original_pos = record_it->second.beginPos + 1;

      if (Options::instance()->stats.size() > 0)
      {
        // First read
        genos1.first.query_name = seqan::toCString(record_it->first.qName);
        genos1.first.details->query_name = seqan::toCString(record_it->first.qName);
        genos1.first.details->read_group = ::get_read_group(record_it->first.tags);

        // Second read
        genos1.second.query_name = seqan::toCString(record_it->second.qName);
        genos1.second.details->query_name = seqan::toCString(record_it->second.qName);
        genos1.second.details->read_group = ::get_read_group(record_it->second.tags);
      }

      genos.push_back(std::move(genos1));
      break;

    case 2:
      genos2.first.forward_strand = false; // The first read in the pair has been reverse complemented
      genos2.first.is_originally_unaligned = seqan::hasFlagUnmapped(record_it->first);
      genos2.second.is_originally_unaligned = seqan::hasFlagUnmapped(record_it->second);
      genos2.first.is_originally_clipped = seqan::hasFlagUnmapped(record_it->first) || is_clipped(record_it->first.cigar);
      genos2.second.is_originally_clipped = seqan::hasFlagUnmapped(record_it->second) || is_clipped(record_it->second.cigar);
      genos2.second.is_first_in_pair = false;
      genos2.first.original_pos = record_it->first.beginPos + 1;
      genos2.second.original_pos = record_it->second.beginPos + 1;

      if (Options::instance()->stats.size() > 0)
      {
        genos2.first.query_name = seqan::toCString(record_it->first.qName);
        genos2.first.details->query_name = seqan::toCString(record_it->first.qName);
        genos2.first.details->read_group = ::get_read_group(record_it->first.tags);

        genos2.second.query_name = seqan::toCString(record_it->second.qName);
        genos2.second.details->query_name = seqan::toCString(record_it->second.qName);
        genos2.second.details->read_group = ::get_read_group(record_it->second.tags);
      }

      genos.push_back(std::move(genos2));
      break;
    }
  }

  return genos;
}


} // namespace gyper
