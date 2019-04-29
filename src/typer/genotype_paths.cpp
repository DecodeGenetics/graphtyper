#include <cassert>
#include <iostream>
#include <string>
#include <vector>

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/graph.hpp>
#include <graphtyper/graph/haplotype_extractor.hpp>
#include <graphtyper/index/kmer_label.hpp>
#include <graphtyper/typer/genotype_paths.hpp>
#include <graphtyper/typer/variant_candidate.hpp>
#include <graphtyper/utilities/kmer_help_functions.hpp>
#include <graphtyper/utilities/options.hpp>

#include <boost/log/trivial.hpp>

#include <seqan/basic.h>
#include <seqan/hts_io.h> // BamAlignmentRecord


namespace
{

bool
merge_current_compatibility(gyper::KmerLabel const & l, gyper::Path const & p)
{
  // It is not compatible if the variant order is not already in the path
  return l.start_index == p.start && l.end_index == p.end;
}


bool
merge_forward_compatibility(gyper::Path const & prev, gyper::Path const & next)
{
  return prev.end == next.start && prev.read_end_index == next.read_start_index;
}


std::vector<gyper::Path>
find_all_nonduplicated_paths(std::vector<gyper::KmerLabel> const & ll,
                             uint32_t const read_start_index,
                             uint32_t const read_end_index,
                             uint16_t const mismatches)
{
  if (ll.size() == 0)
    return std::vector<gyper::Path>(0);

  assert(read_end_index > read_start_index);
  std::vector<gyper::Path> paths(1, gyper::Path(ll.front(), read_start_index, read_end_index, mismatches));

  for (unsigned i = 1; i < ll.size(); ++i)
  {
    bool nothing_found = true;

    for (unsigned d = 0; d < paths.size(); ++d)
    {
      if (merge_current_compatibility(ll[i], paths[d]))
      {
        paths[d].merge_with_current(ll[i]);
        nothing_found = false;
        break;
      }
    }

    if (nothing_found)
    {
      // Nothing was found. This means the label isn't a duplicate, we can add it.
      paths.push_back(gyper::Path(ll[i], read_start_index, read_end_index, mismatches));
    }
  }

  return paths;
}

} // anon namespace


namespace gyper
{

GenotypePaths::GenotypePaths()
  : read(0), paths(0), longest_path_length(0)
{
  if (Options::instance()->stats.size() > 0)
  {
    details = std::unique_ptr<GenotypePathsDetails>(new GenotypePathsDetails);
  }
}


GenotypePaths::GenotypePaths(seqan::IupacString const & _read,
                             seqan::CharString const & _qual,
                             uint8_t const _mapq
                             )
  : read(seqan::begin(_read), seqan::end(_read))
  , qual(seqan::begin(_qual), seqan::end(_qual))
  , longest_path_length(0)
  , mapq(_mapq)
{
  if (Options::instance()->stats.size() > 0)
  {
    details = std::unique_ptr<GenotypePathsDetails>(new GenotypePathsDetails);
  }
}


bool
GenotypePaths::all_paths_unique() const
{
  for (std::size_t i = 1; i < paths.size(); ++i)
  {
    if (paths[0].start_ref_reach_pos() != paths[i].start_ref_reach_pos() &&
        paths[0].end_ref_reach_pos() != paths[i].end_ref_reach_pos()
        )
    {
      return false;
    }
  }

  return true;
}


void
GenotypePaths::add_prev_kmer_labels(std::vector<KmerLabel> const & ll,
                                    uint32_t const read_start_index,
                                    uint32_t const read_end_index,
                                    int const mismatches
                                    )
{
  assert(read_end_index > read_start_index);
  assert(mismatches >= 0);
  std::vector<Path> const pp = find_all_nonduplicated_paths(ll,
                                                            read_start_index,
                                                            read_end_index,
                                                            (uint16_t)mismatches);
  std::size_t const original_size = paths.size();

  // Keep track of which new paths did matches with any previous paths
  std::vector<uint8_t> new_path_matched(pp.size(), 0u); // 1 == true, 0 == false

  for (unsigned i = 0; i < original_size; ++i)
  {
    if (paths[i].read_start_index != read_end_index)
      continue;

    bool matched_at_least_once = false;
    Path original_path(paths[i]);

    for (unsigned j = 0; j < pp.size(); ++j)
    {
      // Check if the two paths are compatible
      if (merge_forward_compatibility(pp[j], original_path))
      {
        Path np(pp[j], original_path);

        // This checks if the path could be merged
        if (np.read_start_index != pp[j].read_start_index)
          continue;

        // Merge success
        new_path_matched[j] = 1u;

        if (matched_at_least_once)
        {
          assert(np.size() <= longest_path_length);
          paths.push_back(std::move(np));
        }
        else
        {
          longest_path_length = std::max(np.size(), longest_path_length);
          paths[i] = std::move(np);
          matched_at_least_once = true;
        }
      }
    }
  }

  for (unsigned j = 0; j < new_path_matched.size(); ++j)
  {
    if (new_path_matched[j] == 0u)
    {
      longest_path_length = std::max(pp[j].size(), longest_path_length);
      paths.push_back(pp[j]);
    }
  }
}


void
GenotypePaths::add_next_kmer_labels(std::vector<KmerLabel> const & ll,
                                    uint32_t const read_start_index,
                                    uint32_t const read_end_index,
                                    int const mismatches
                                    )
{
  assert(read_end_index > read_start_index);
  std::vector<Path> const pp = find_all_nonduplicated_paths(ll,
                                                            read_start_index,
                                                            read_end_index,
                                                            (uint16_t)mismatches
                                                            );

  std::size_t const original_size = paths.size();

  // Keep track of which new paths did matches with any previous paths
  std::vector<uint8_t> new_path_matched(pp.size(), 0u); // 1 == true, 0 == false

  for (unsigned i = 0; i < original_size; ++i)
  {
    if (paths[i].read_end_index != read_start_index)
      continue;

    bool matched_at_least_once = false;
    Path original_path(paths[i]);

    for (unsigned j = 0; j < pp.size(); ++j)
    {
      // Check if the two paths are compatible
      if (merge_forward_compatibility(original_path, pp[j]))
      {
        Path np(original_path, pp[j]);

        // This checks if the path could be merged
        if (np.start != original_path.start || np.read_start_index != original_path.read_start_index)
          continue;

        new_path_matched[j] = 1u; // Merging was a success

        if (matched_at_least_once)
        {
          assert(np.size() <= longest_path_length);
          paths.push_back(std::move(np));
        }
        else
        {
          longest_path_length = std::max(np.size(), longest_path_length);
          paths[i] = std::move(np);
          matched_at_least_once = true;
        }
      }
    }
  }

  for (unsigned j = 0; j < new_path_matched.size(); ++j)
  {
    if (new_path_matched[j] == 0u)
    {
      longest_path_length = std::max(pp[j].size(), longest_path_length);
      paths.push_back(pp[j]);
    }
  }
}


void
GenotypePaths::clear_paths()
{
  paths.clear();
  longest_path_length = 0;
}


void
GenotypePaths::remove_paths_with_too_many_mismatches()
{
  assert(paths.size() > 0 || longest_path_length == 0);

  if (paths.size() == 0)
    return;

  // Maximum mismatches allowed
  uint16_t min_mismatches = Options::instance()->is_segment_calling ? (short)2 : (short)10;

  // Find the minimum number of mismatches in the aligned paths
  for (auto const & path : paths)
    min_mismatches = std::min(path.mismatches, min_mismatches);

  // Erase paths with more mismatches than the minimum number of mismatches
  paths.erase(std::remove_if(paths.begin(), paths.end(), [min_mismatches](Path const & path){
      return path.mismatches > min_mismatches;
    }), paths.end());

  // Finally, update the longest path size
  update_longest_path_size();
}


void
GenotypePaths::remove_support_from_read_ends()
{
  long constexpr MIN_OFFSET = 4;

  for (Path & path : paths)
  {
    if (path.var_order.size() == 0)
      continue;

    if (!graph.is_special_pos(path.start) && !graph.is_special_pos(path.end))
      continue;

    auto min_max_elements = std::minmax_element(path.var_order.begin(), path.var_order.end());
    assert(min_max_elements.first != path.var_order.end());
    assert(min_max_elements.second != path.var_order.end());

    // Check end position
    if (graph.is_special_pos(path.end) && path.end_correct_pos() <= (*min_max_elements.second) + MIN_OFFSET)
    {
      long const index = std::distance(path.var_order.begin(), min_max_elements.second);
      assert(index < static_cast<long>(path.nums.size()));
      path.nums[index].set();
    }

    // Check start position
    if (graph.is_special_pos(path.start))
    {
      bool is_ambigous = false;

      if (graph.is_special_pos(path.start + static_cast<uint32_t>(MIN_OFFSET)))
      {
        long const start_ref_reach_pos = path.start_ref_reach_pos();
        long const start_offset_ref_reach_pos =
          graph.get_ref_reach_pos(path.start + static_cast<uint32_t>(MIN_OFFSET));
        is_ambigous = start_ref_reach_pos != start_offset_ref_reach_pos;
      }
      else
      {
        is_ambigous = true;
      }

      if (is_ambigous)
      {
        long const index = std::distance(path.var_order.begin(), min_max_elements.first);
        assert(index < static_cast<long>(path.nums.size()));
        path.nums[index].set();
      }
    }
  }
}


void
GenotypePaths::extend_paths_at_sv_breakpoints(seqan::IupacString const & seqan_read)
{
  if (paths.size() > 5)
    return;

  for (Path & path : paths)
  {
    if (path.size() < 60 || path.mismatches > 2)
    {
      // Require not too small matches or long matches with many mismatches
    }
    else if (path.read_start_index > 0 && path.read_end_index < (read.size() - 1))
    {
      // Ignore reads that are clipped on both sides
    }
    else if (path.read_start_index > 0)
    {
      std::vector<Location> s_locs = graph.get_locations_of_a_position(path.start, path);
      std::vector<Location> e_locs = graph.get_locations_of_a_position(path.end, path);
      bool is_match_found = false;

      for (auto const & s : s_locs)
      {
        if (s.node_type != 'V')
          continue;

        for (auto const & e : e_locs)
        {
          if (e.node_type == 'V' &&
            s.node_order == e.node_order &&
            path.read_start_index > s.offset
             )
          {
            // Require a kmer match before the variant node with first kmer in the read
            std::vector<KmerLabel> labels = query_index_for_first_kmer(seqan_read);

            for (KmerLabel const & label : labels)
            {
              if (label.start_index < s.node_order && label.start_index + 100 >= s.node_order)
              {
                // We found a kmer before the variant node, now we trust that the read belongs here!
                is_match_found = true;
                break; // We only need to find a single match
              }
            }

            // Assume that the breakpoint is incorrect and align to it (but not further)
            if (is_match_found)
            {
              assert (s.offset < e.offset);
              path.read_start_index -= s.offset;
              assert(path.read_start_index > 0);
              path.start = s.node_order;
              path.mismatches += 2;
              break;
            } // is_match_found
          }
        }

        if (is_match_found)
          break;
      }
    }
    else if (path.read_end_index < (read.size() - 1))
    {
      std::vector<Location> e_locs = graph.get_locations_of_a_position(path.end, path);

      for (auto const & e : e_locs)
      {
        if (e.node_type != 'R')
          continue;

        RefNode const & ref_node = graph.ref_nodes[e.node_index];

        if (ref_node.out_degree() <= 1)
          continue; // Final reference node

        // Check if the read can possibly reach the breakpoint
        long const remaining_read = read.size() - 1 - path.read_end_index;
        long const remaining_ref = ref_node.get_label().dna.size() - e.offset;

        if (remaining_read < remaining_ref)
          continue;

        auto const next_var_index = ref_node.get_var_index(1);
        VarNode const & next_var_node = graph.var_nodes[next_var_index];

        if (next_var_node.get_label().dna.size() < 150) // Not an SV
          continue;

        // Require a kmer match on the variant node with the very last kmer in the read
        std::vector<KmerLabel> labels = query_index_for_last_kmer(seqan_read);
        bool is_match_found = false;

        for (KmerLabel const & label : labels)
        {
          if (label.variant_id == next_var_index)
          {
            // We found a kmer in the next variant node, trust that the read belongs here!
            is_match_found = true;
            break; // We only need to find a single match
          }
        }

        if (is_match_found)
        {
          assert(read.size() > 0);
          path.read_end_index += remaining_ref;
          assert(path.read_end_index <= read.size() - 1);
          path.end = next_var_node.get_label().order;
          path.var_order.push_back(next_var_node.get_label().order);
          path.mismatches += 2;

          // Set bitset
          {
            std::bitset<MAX_NUMBER_OF_HAPLOTYPES> num;

            for (long i = 1; i < static_cast<long>(ref_node.out_degree()); ++i)
              num.set(i);

            path.nums.push_back(num);
          }

          break;
        } // is_match_found
      }
    }
  }
}


void
GenotypePaths::remove_paths_within_variant_node()
{
  if (paths.size() <= 1)
    return;

  auto is_path_within_one_variant_node = [&](Path const & path)
  {
    std::vector<Location> s_locs = graph.get_locations_of_a_position(path.start, path);
    std::vector<Location> e_locs = graph.get_locations_of_a_position(path.end, path);

    for (auto const & s : s_locs)
    {
      if (s.node_type != 'V')
        continue;

      for (auto const & e : e_locs)
      {
        if (e.node_type == 'V' && s.node_order == e.node_order && s.offset > 0)
        {
          return true;
        }
      }
    }

    return false;
  };

  // Erase paths which are within a single variant node
  paths.erase(std::remove_if(paths.begin(),
                             paths.end(),
                             is_path_within_one_variant_node
                             ), paths.end()
    );

  // Finally, update the longest path size
  update_longest_path_size();
}


void
GenotypePaths::remove_non_ref_paths_when_read_matches_ref()
{
  if (all_paths_unique()) // paths.size() == 0
    return;

  // Check if there are any paths that support purely reference, and in that case delete every other path.
  auto find_path_it = std::find_if(paths.begin(), paths.end(), [](Path const & p){return p.is_reference();});

  if (find_path_it != paths.end())
  {
    // I found a path that is purely reference... delete all other paths
    paths.erase(std::remove_if(paths.begin(), paths.end(), [](Path const & p){return !p.is_reference();}), paths.end());
  }
}


void
GenotypePaths::remove_fully_special_paths()
{
  auto is_fully_special = [](Path const & p) -> bool
  {
    return p.start_ref_reach_pos() == p.end_ref_reach_pos();
  };

  paths.erase(std::remove_if(paths.begin(),
                             paths.end(),
                             is_fully_special
              )
    , paths.end());
}


void
GenotypePaths::walk_read_ends(seqan::IupacString const & seq, int maximum_mismatches, gyper::Graph const & graph)
{
  if (paths.size() == 0 || paths[0].size() == seqan::length(seq))
    return;

  if (paths.size() > Options::instance()->MAX_SEED_NUMBER_FOR_WALKING)
    return; // Do not walk if we have too many seeds

  if (paths.size() > Options::instance()->MAX_SEED_NUMBER_ALLOWING_MISMATCHES)
    maximum_mismatches = 0; // Only allow exact matches when we have a lot of paths

  std::size_t best_mismatches = 7;
  std::vector<uint32_t> best_end_indexes;
  std::vector<std::vector<KmerLabel> > best_labels;

  for (auto & path : paths)
  {
    // Check if this path can be made longer
    assert(path.read_end_index <= seqan::length(seq) - 1);

    if (path.read_end_index == seqan::length(seq) - 1)
      continue;

    std::vector<Location> s_locs = graph.get_locations_of_a_position(path.end, path);

    if (s_locs.size() == 0 || s_locs.size() > Options::instance()->MAX_NUM_LOCATIONS_PER_PATH)
      continue;

    std::vector<char> kmer;
    kmer.reserve(seqan::length(seq) - path.read_end_index + 1); // It cannot get bigger than this

    for (uint32_t i = path.read_end_index; i < seqan::length(seq); ++i)
      kmer.push_back(seq[i]);

    std::vector<Location> e_locs(1); // Unavailable end
    std::vector<KmerLabel> new_labels;

    // Value less than zero causes default value
    uint32_t mismatches = (maximum_mismatches < 0) ?
      static_cast<uint32_t>(std::min(2 + kmer.size() / 11, best_mismatches)) :
      static_cast<uint32_t>(maximum_mismatches);

    new_labels = graph.iterative_dfs(std::move(s_locs), std::move(e_locs), kmer, mismatches);

    if (new_labels.size() > 0)
    {
      // Add the results
      if (mismatches < best_mismatches)
      {
        best_labels.clear();
        best_labels.push_back(std::move(new_labels));
        best_end_indexes.clear();
        best_end_indexes.push_back(path.read_end_index);
        best_mismatches = mismatches;
      }
      else if (mismatches == best_mismatches)
      {
        best_labels.push_back(std::move(new_labels));
        best_end_indexes.push_back(path.read_end_index);
      }
    }
  }

  if (best_labels.size() > 0)
  {
    for (unsigned i = 0; i < best_labels.size(); ++i)
    {
      add_next_kmer_labels(best_labels[i],
                           best_end_indexes[i],
                           seqan::length(seq) - 1,
                           (int)best_mismatches
        );
    }
  }
}


void
GenotypePaths::walk_read_starts(seqan::IupacString const & seq, int maximum_mismatches, gyper::Graph const & graph)
{
  if (paths.size() == 0 || paths[0].size() == seqan::length(seq))
    return;

  if (paths.size() > Options::instance()->MAX_SEED_NUMBER_FOR_WALKING)
    return; // Do not walk if we have too many seeds

  if (paths.size() > Options::instance()->MAX_SEED_NUMBER_ALLOWING_MISMATCHES)
    maximum_mismatches = 0; // Only allow exact matches when we have a lot of paths

  std::size_t best_mismatches = 7;
  std::vector<std::vector<KmerLabel> > best_labels;
  std::vector<uint32_t> best_start_indexes;

  for (auto & path : paths)
  {
    // Check if the path can be lengthened
    if (path.read_start_index == 0)
      continue;

    assert(path.read_start_index % (K - 1) == 0);
    assert(path.read_start_index >= (K - 1));
    std::vector<char> kmer;
    kmer.reserve(path.read_start_index + 1); // It cannot get bigger than this

    for (uint32_t i = 0; i <= path.read_start_index; ++i)
      kmer.push_back(seq[i]);

    std::vector<Location> e_locs = graph.get_locations_of_a_position(path.start, path);

    if (e_locs.size() == 0 || e_locs.size() > Options::instance()->MAX_NUM_LOCATIONS_PER_PATH)
      continue;

    std::vector<Location> s_locs(1); // Unavailable start
    std::vector<KmerLabel> new_labels;

    // Value less than zero causes default value
    uint32_t mismatches = (maximum_mismatches < 0) ?
      std::min(2 + kmer.size() / 11, best_mismatches) :
      static_cast<uint32_t>(maximum_mismatches);

    new_labels = graph.iterative_dfs(std::move(s_locs), std::move(e_locs), kmer, mismatches);

    if (new_labels.size() > 0)
    {
      // Add the results
      if (mismatches < best_mismatches)
      {
        best_labels.clear();
        best_labels.push_back(std::move(new_labels));
        best_start_indexes.clear();
        best_start_indexes.push_back(path.read_start_index);
        best_mismatches = mismatches;
      }
      else if (mismatches == best_mismatches)
      {
        best_labels.push_back(std::move(new_labels));
        best_start_indexes.push_back(path.read_start_index);
      }
    }
  }

  if (best_labels.size() > 0)
  {
    for (unsigned i = 0; i < best_labels.size(); ++i)
      add_prev_kmer_labels(best_labels[i], 0, best_start_indexes[i], best_mismatches);
  }
}


std::vector<VariantCandidate>
GenotypePaths::find_new_variants() const
{
  std::vector<VariantCandidate> new_variants;

  // Don't try to find variants in perfect reads
  if (paths.size() == 0 || (all_paths_fully_aligned() && paths[0].mismatches == 0))
    return new_variants;

  auto const & path = paths[0];
  // Check if the paths is fully aligned (with mismatches) and is purely on the reference,
  // in this case we assume all the variants are SNPs
  if (all_paths_fully_aligned() && is_purely_reference())
  {
    // Discover SNPs
    uint32_t pos = path.start_ref_reach_pos();
    uint32_t end_pos = path.end_ref_reach_pos() + 1;
    std::vector<char> const reference = graph.get_generated_reference_genome(pos, end_pos);
    assert(pos == path.start);
    assert(end_pos == path.end_pos() + 1);
    assert(reference.size() == read.size());

    {
      int matches = -1; // -1 means we have not found a mismatch yet
      std::vector<char> ref;
      std::vector<char> alt;
      int64_t const MIN_VAR_THRESHOLD = 5;
      uint32_t var_pos = 0;

      for (unsigned i = 0; i < reference.size(); ++i)
      {
        assert(i < seqan::length(read));
        assert(i < seqan::length(qual));

        if (reference[i] != read[i] && read[i] != 'N')
        {
          if (matches == -1)
            var_pos = pos + i;

          matches = 0;
        }
        else if (matches >= 0)
        {
          ++matches;
        }

        if (matches >= 0)
        {
          ref.push_back(reference[i]);
          alt.push_back(read[i]);
        }

        if (matches >= MIN_VAR_THRESHOLD)
        {
          assert(matches == MIN_VAR_THRESHOLD);
          assert(ref.size() >= static_cast<std::size_t>(MIN_VAR_THRESHOLD));
          assert(alt.size() >= MIN_VAR_THRESHOLD);
          assert(var_pos != 0);

          // Create and add the new variant
          {
            VariantCandidate new_var;
            new_var.abs_pos = var_pos;
            std::vector<char> new_ref(ref.begin(), ref.end() - matches);
            std::vector<char> new_alt(alt.begin(), alt.end() - matches);

            // Make sure there are no Ns
            if (std::find(new_ref.begin(), new_ref.end(), 'N') == new_ref.end() &&
                std::find(new_alt.begin(), new_alt.end(), 'N') == new_alt.end()
                )
            {
              new_var.seqs.push_back(std::move(new_ref));
              new_var.seqs.push_back(std::move(new_alt));

              // Determine if the variant is low quality
              {
                unsigned const r = 1u + i - matches - new_var.seqs[0].size();
                assert(new_var.seqs[1].size() > 0);
                unsigned const MAX_QUAL = *std::max_element(qual.begin() + r, qual.begin() + r + new_var.seqs[1].size()) - 33u;
                new_var.is_low_qual = MAX_QUAL < 25u;
              }

              // Determine if it is a proper pair
              new_var.normalize();
              new_variants.push_back(std::move(new_var));
            }
          }

          ref.clear();
          alt.clear();
          matches = -1;
          var_pos = 0;
        }
      }

      // Leftovers
      if (matches >= 0)
      {
        assert(var_pos != 0);

        // Create and add the new variant
        {
          VariantCandidate new_var;
          new_var.abs_pos = var_pos;
          std::vector<char> new_ref(ref.begin(), ref.end() - matches);
          std::vector<char> new_alt(alt.begin(), alt.end() - matches);

          // Make sure there are no Ns
          if (std::find(new_ref.begin(), new_ref.end(), 'N') == new_ref.end() &&
              std::find(new_alt.begin(), new_alt.end(), 'N') == new_alt.end()
              )
          {
            new_var.seqs.push_back(std::move(new_ref));
            new_var.seqs.push_back(std::move(new_alt));

            // Determine if the variant is low quality
            {
              unsigned const r = read.size() - matches - new_var.seqs[0].size();
              assert(new_var.seqs[1].size() > 0);
              unsigned const MAX_QUAL = *std::max_element(qual.begin() + r, qual.begin() + r + new_var.seqs[1].size()) - 33u;
              new_var.is_low_qual = MAX_QUAL < 25u;
            }

            // Make sure the sequences are not empty
            assert(new_var.seqs[0].size() > 0);
            assert(new_var.seqs[1].size() > 0);

            new_var.normalize();
            new_variants.push_back(std::move(new_var));
          }
        }
      }
    }
  }
  else
  {
    // Discover SNPs and indels
    uint32_t pos = path.start_ref_reach_pos() - path.read_start_index;

    // Parameters
    uint32_t constexpr EXTRA_BASES_BEFORE = 50;
    uint32_t constexpr EXTRA_BASES_AFTER = 50;
    uint32_t ref_pos_start;

    // Check if we underflew, and if we didn't prevent an underflow
    if (pos <= path.start_ref_reach_pos() && pos > EXTRA_BASES_BEFORE)
      ref_pos_start = pos - EXTRA_BASES_BEFORE;
    else
      ref_pos_start = 0;

    uint32_t ref_pos_end = static_cast<uint32_t>(pos + read.size() + EXTRA_BASES_AFTER);
    std::vector<char> reference = graph.get_generated_reference_genome(ref_pos_start, ref_pos_end);

    // Disocvery new variants (SNPs and indels)
    new_variants = find_variants_in_alignment(ref_pos_start,
                                              std::move(reference),
                                              static_cast<seqan::Dna5String>(read),
                                              qual
      );
  }

  for (auto & new_var : new_variants)
  {
    // Make sure the sequences are not empty
    assert(new_var.seqs.size() == 2);
    assert(new_var.seqs[0].size() > 0);
    assert(new_var.seqs[1].size() > 0);
    new_var.normalize();
    new_var.is_in_proper_pair = this->is_proper_pair();
    new_var.is_mapq0 = mapq == 0;
    new_var.is_unaligned = is_originally_unaligned;
    new_var.is_clipped = is_originally_clipped;
    new_var.is_first_in_pair = is_first_in_pair;
    new_var.is_seq_reversed = not forward_strand;
    new_var.original_pos = original_pos;
  }

  return new_variants;
}


void
GenotypePaths::remove_short_paths()
{
  if (longest_path_length <= 1)
    return;

  paths.erase(std::remove_if(paths.begin(), paths.end(), [&](Path const & p){
      return p.size() < longest_path_length;
    }), paths.end());

  assert(paths.size() > 0);
}


/*
void
GenotypePaths::remove_paths_with_no_variants()
{
  paths.erase(std::remove_if(paths.begin(), paths.end(), [](Path const & p){
      return p.var_order.size() == 0;
    }), paths.end());
}
*/


bool
GenotypePaths::all_paths_fully_aligned() const
{
  for (auto const & path : paths)
  {
    if (path.size() != this->read.size())
      return false;
  }

  return true;
}


bool
GenotypePaths::is_purely_reference() const
{
  for (auto const & path : paths)
  {
    if (!path.is_purely_reference())
      return false;
  }

  return true;
}


void
GenotypePaths::update_longest_path_size()
{
  longest_path_length = 0;

  for (auto const & path : paths)
    longest_path_length = std::max(path.size(), longest_path_length);
}


/*********************
 * CLASS INFORMATION *
 *********************/

std::size_t
GenotypePaths::longest_path_size() const
{
  return longest_path_length;
}


std::vector<Path>
GenotypePaths::longest_paths() const
{
  std::size_t local_longest_path_length = longest_path_length;
  std::vector<Path> long_paths;

  for (auto path_it = paths.cbegin(); path_it != paths.cend(); ++path_it)
  {
    if (path_it->size() > local_longest_path_length)
    {
      local_longest_path_length = path_it->size();
      long_paths.clear();
      long_paths.push_back(*path_it);
    }
    else if (path_it->size() == local_longest_path_length)
    {
      long_paths.push_back(*path_it);
    }
  }

  assert(paths.size() == 0 || long_paths.size() > 0);
  return long_paths;
}


bool
GenotypePaths::check_no_variant_is_missing() const
{
  for (auto const & path : paths)
  {
    std::vector<uint32_t> expected_orders = graph.get_var_orders(path.start_ref_reach_pos(), path.end_ref_reach_pos());

    if (expected_orders.size() != path.var_order.size())
    {
      BOOST_LOG_TRIVIAL(error) << "[genotype_paths] The number of expected orders did not match, got " << path.var_order.size() << " but expected " << expected_orders.size();
      return false;
    }

    auto all_orders_match = [&](std::vector<uint32_t> const & o1, std::vector<uint32_t> const & o2)
                            {
                              for (auto const & o : o1)
                              {
                                if (std::find(o2.begin(), o2.end(), o) == o2.end())
                                  return false;
                              }

                              return true;
                            };

    if (!all_orders_match(expected_orders, path.var_order) || !all_orders_match(path.var_order, expected_orders))
    {
      BOOST_LOG_TRIVIAL(error) << "[genotype_paths] Orders did not match with the expected orders.";
      return false;
    }
  }

  return true;
}


std::string
GenotypePaths::to_string() const
{
  std::ostringstream ss;
  ss << "read_name=" << details->query_name
     << " first_in_pair=" << (is_first_in_pair ? "Y" : "N")
     << " read=" << std::string(read.begin(), read.end())
     << " paths.size()=" << paths.size()
     << " pos=" << original_pos
     << " mapq=" << static_cast<unsigned>(mapq)
     << "\n";

  for (auto const & path : paths)
  {
    ss << "  path_start=" << path.start_pos()
       << " path_end=" << path.end_pos()
       << " path_start_correct=" << path.start_correct_pos()
       << " path_end_correct=" << path.end_correct_pos()
       << " read_start_index=" << path.read_start_index
       << " read_end_index=" << path.read_end_index
       << " mismatches=" << path.mismatches
       << " path.var_order.size()=" << path.var_order.size()
       << " var_orders=";

    for (auto const var_order : path.var_order)
      ss << var_order << " ";

    ss << "\n";
  }

  return ss.str();
}


bool
GenotypePaths::is_proper_pair() const
{
  return ml_insert_size != INSERT_SIZE_WHEN_NOT_PROPER_PAIR;
}


int
compare_pair_of_genotype_paths(GenotypePaths const & geno1, GenotypePaths const & geno2)
{
  std::size_t const TOTAL_MATCHES_1 = geno1.longest_path_size();
  std::size_t const TOTAL_MATCHES_2 = geno2.longest_path_size();
  std::size_t const MINIMUM_PATH_SIZE = 90;

  if (TOTAL_MATCHES_1 > TOTAL_MATCHES_2 && TOTAL_MATCHES_1 > MINIMUM_PATH_SIZE)
  {
    return 1;
  }
  else if (TOTAL_MATCHES_2 > TOTAL_MATCHES_1 && TOTAL_MATCHES_2 > MINIMUM_PATH_SIZE)
  {
    return 2;
  }
  else if (TOTAL_MATCHES_2 == TOTAL_MATCHES_1 && TOTAL_MATCHES_1 > MINIMUM_PATH_SIZE)
  {
    std::size_t const mismatches1 = geno1.paths[0].mismatches;
    std::size_t const mismatches2 = geno2.paths[0].mismatches;

    if (mismatches1 < mismatches2)
      return 1;
    else if (mismatches2 < mismatches1)
      return 2;
    else
      return 1;
  }

  return 0;
}


int
compare_pair_of_genotype_paths(std::pair<GenotypePaths, GenotypePaths> const & genos1, std::pair<GenotypePaths, GenotypePaths> const & genos2)
{
  std::size_t const TOTAL_MATCHES_1_1 = genos1.first.paths.size() > 0 ? genos1.first.longest_path_size() : 0;
  std::size_t const TOTAL_MATCHES_1_2 = genos1.second.paths.size() > 0 ? genos1.second.longest_path_size() : 0;
  std::size_t const TOTAL_MATCHES_2_1 = genos2.first.paths.size() > 0 ? genos2.first.longest_path_size() : 0;
  std::size_t const TOTAL_MATCHES_2_2 = genos2.second.paths.size() > 0 ? genos2.second.longest_path_size() : 0;
  std::size_t const MAX_SIZE_1 = std::max(TOTAL_MATCHES_1_1, TOTAL_MATCHES_1_2);
  std::size_t const MAX_SIZE_2 = std::max(TOTAL_MATCHES_2_1, TOTAL_MATCHES_2_2);
  std::size_t const PERFECT_PATH_SIZE_1 = genos1.first.read.size();
  assert(PERFECT_PATH_SIZE_1 == genos2.first.read.size());
  std::size_t const PERFECT_PATH_SIZE_2 = genos1.second.read.size();
  assert(PERFECT_PATH_SIZE_2 == genos2.second.read.size());
  std::size_t const MINIMUM_PATH_SIZE = 90;


  // Paths will be chosen in the following order:
  //   1. Alignment with the longest match (which is at least of length MINIMUM_PATH_SIZE)
  //   2. Number of mismatches in the best matches.
  //   3. The length of the worse alignment.
  //   4. The alignment with fewer paths.
  //   5. The alignment with fewer non-ref variant calls.
  // If we haven't already chosen an alignment at this point we throw both of them away.
  // Also, the requirement for an alignment to be considered is that both reads have at least two K-mers aligned or one with 90 bases matched.

  if ((TOTAL_MATCHES_1_1 >= PERFECT_PATH_SIZE_1 && TOTAL_MATCHES_1_2 >= PERFECT_PATH_SIZE_2) ||
      (TOTAL_MATCHES_2_1 >= PERFECT_PATH_SIZE_1 && TOTAL_MATCHES_2_2 >= PERFECT_PATH_SIZE_2)
      )
  {
    // We have a perfect match
    if ((TOTAL_MATCHES_1_1 >= PERFECT_PATH_SIZE_1 && TOTAL_MATCHES_1_2 >= PERFECT_PATH_SIZE_2) &&
        (TOTAL_MATCHES_2_1 >= PERFECT_PATH_SIZE_1 && TOTAL_MATCHES_2_2 >= PERFECT_PATH_SIZE_2)
        )
    {
      assert(genos1.first.paths.size() > 0);
      assert(genos1.second.paths.size() > 0);
      assert(genos2.first.paths.size() > 0);
      assert(genos2.second.paths.size() > 0);
      std::size_t const mismatches1 = genos1.first.paths[0].mismatches + genos1.second.paths[0].mismatches;
      std::size_t const mismatches2 = genos2.first.paths[0].mismatches + genos2.second.paths[0].mismatches;

      if (mismatches1 < mismatches2)
      {
        return 1;
      }
      else if (mismatches2 < mismatches1)
      {
        return 2;
      }
      else
      {
        std::size_t const num_paths1 = genos1.first.paths.size() + genos1.second.paths.size();
        std::size_t const num_paths2 = genos2.first.paths.size() + genos2.second.paths.size();

        if (num_paths1 < num_paths2)
        {
          return 1;
        }
        else if (num_paths2 < num_paths1)
        {
          return 2;
        }
        else
        {
          // Count the number of alternative allele calls
          auto alternative_call_count = [](std::vector<Path> const & paths)
                                        {
                                          std::size_t count = 0;

                                          for (auto const & path : paths)
                                          {
                                            for (auto const & num : path.nums)
                                            {
                                              count += (!num.test(0));
                                            }
                                          }

                                          return count;
                                        };

          std::size_t const COUNT_1 = alternative_call_count(genos1.first.paths) + alternative_call_count(genos1.second.paths);
          std::size_t const COUNT_2 = alternative_call_count(genos2.first.paths) + alternative_call_count(genos2.second.paths);

          if (COUNT_1 >= COUNT_2)
          {
            return 1;
          }
          else
          {
            return 2;
          }
        }
      }
    }
    else if (TOTAL_MATCHES_1_1 >= PERFECT_PATH_SIZE_1 && TOTAL_MATCHES_1_2 >= PERFECT_PATH_SIZE_2)
    {
      return 1;
    }
    else
    {
      // Implies TOTAL_MATCHES_2_1 >= PERFECT_PATH_SIZE_1 && TOTAL_MATCHES_2_2 >= PERFECT_PATH_SIZE_2
      return 2;
    }
  }
  else if (MAX_SIZE_2 >= MINIMUM_PATH_SIZE && MAX_SIZE_2 > MAX_SIZE_1)
  {
    return 2;
  }
  else if (MAX_SIZE_1 >= MINIMUM_PATH_SIZE && MAX_SIZE_1 > MAX_SIZE_2)
  {
    return 1;
  }
  else if (MAX_SIZE_1 >= MINIMUM_PATH_SIZE && MAX_SIZE_2 >= MINIMUM_PATH_SIZE)
  {
    assert(MAX_SIZE_1 == MAX_SIZE_2); // This should be implied
    uint16_t mismatches1 = 10;

    // Get mismatches for the first pair of genotype paths
    {
      if (TOTAL_MATCHES_1_1 == MAX_SIZE_1)
      {
        assert(genos1.first.paths.size() > 0);
        mismatches1 = std::min(mismatches1, genos1.first.paths[0].mismatches);
      }

      if (TOTAL_MATCHES_1_2 == MAX_SIZE_1)
      {
        assert(genos1.second.paths.size() > 0);
        mismatches1 = std::min(mismatches1, genos1.second.paths[0].mismatches);
      }
    }

    uint16_t mismatches2 = 10;

    // Get mismatches for the second pair of genotype paths
    {
      if (TOTAL_MATCHES_2_1 == MAX_SIZE_2)
      {
        assert(genos2.first.paths.size() > 0);
        mismatches2 = std::min(mismatches2, genos2.first.paths[0].mismatches);
      }

      if (TOTAL_MATCHES_2_2 == MAX_SIZE_2)
      {
        assert(genos2.second.paths.size() > 0);
        mismatches2 = std::min(mismatches2, genos2.second.paths[0].mismatches);
      }
    }

    if (mismatches1 < mismatches2)
    {
      return 1;
    }
    else if (mismatches2 < mismatches1)
    {
      return 2;
    }
    else
    {
      // Same number of mismatches
      if (std::min(TOTAL_MATCHES_1_1, TOTAL_MATCHES_1_2) < std::min(TOTAL_MATCHES_2_1, TOTAL_MATCHES_2_2))
        return 1;
      else if (std::min(TOTAL_MATCHES_2_1, TOTAL_MATCHES_2_2) < std::min(TOTAL_MATCHES_1_1, TOTAL_MATCHES_1_2))
        return 2;
      else
        return 0;
    }
  }
  else if (MAX_SIZE_2 == 0u && TOTAL_MATCHES_1_1 >= 63u && TOTAL_MATCHES_1_2 >= 63u)
  {
    return 1;
  }
  else if (MAX_SIZE_1 == 0u && TOTAL_MATCHES_2_1 >= 63u && TOTAL_MATCHES_2_2 >= 63u)
  {
    return 2;
  }

  return 0;
}


} // namespace gyper
