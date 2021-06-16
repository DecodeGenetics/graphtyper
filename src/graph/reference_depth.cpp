#include <algorithm> // std::sort
#include <sstream>   // std::stringstream
#include <string>    // std::string
#include <vector>    // std::vector

#include <graphtyper/graph/graph.hpp>
#include <graphtyper/graph/reference_depth.hpp>
#include <graphtyper/typer/genotype_paths.hpp>
#include <graphtyper/typer/variant_candidate.hpp>
#include <graphtyper/utilities/logging.hpp>

namespace gyper
{
/***********************
 * Global reference
 */

ReferenceDepth::ReferenceDepth()
{
  reference_offset = graph.ref_nodes.size() > 0 ? graph.ref_nodes[0].get_label().order : 0;
}

void ReferenceDepth::set_depth_sizes(long const sample_count, long const reference_size)
{
  std::vector<uint16_t> depth;
  depth.resize(reference_size);
  depths.resize(sample_count, depth);
}

uint16_t ReferenceDepth::get_read_depth(VariantCandidate const & var, long const sample_index) const
{
  assert(depths.size() > 0);
  assert(sample_index < static_cast<long>(depths.size()));
  auto const & depth = depths[sample_index];
  assert(depth.size() > 0);
  assert(var.seqs.size() > 0);
  uint32_t start_pos = var.abs_pos;
  uint32_t end_pos = static_cast<uint32_t>(start_pos + var.seqs[0].size() - 1);

  // We want to avoid getting the coverage at the first pos if possible, since the first positions very often match in
  // more than one path
  if (var.seqs[0].size() > 1)
    ++start_pos;

  long const start_index = start_pos_to_index(start_pos);
  long const end_index = end_pos_to_index(end_pos, depth.size());

  if (start_index < static_cast<long>(depth.size()))
  {
    auto max_depth_it = std::max_element(depth.cbegin() + start_index, depth.cbegin() + end_index);
    assert(max_depth_it != depth.begin() + end_index);
    return *max_depth_it;
  }
  else
  {
    return 0;
  }
}

uint16_t ReferenceDepth::get_read_depth(uint32_t abs_pos, long const sample_index) const
{
  assert(sample_index < static_cast<long>(depths.size()));
  auto const & depth = depths[sample_index];

  if (depth.size() == 0)
  {
    // Happens if there are absolutely no reads in the region
    return 0;
  }
  else
  {
    // std::cerr << "abs_pos = " << abs_pos << " " << start_pos_to_index(abs_pos) << "\n";
    long const index = std::min(start_pos_to_index(abs_pos), static_cast<long>(depth.size() - 1l));
    return depth[index];
  }
}

uint64_t ReferenceDepth::get_total_read_depth_of_samples(VariantCandidate const & var,
                                                         std::vector<uint32_t> const & sample_indexes) const
{
  assert(var.seqs.size() > 0);
  uint32_t start_pos = var.abs_pos;
  uint32_t end_pos = start_pos + var.seqs[0].size() - 1;

  // We want to avoid getting the coverage at the first pos if possible, since the first positions very often
  // match in more than one path
  if (var.seqs[0].size() > 1)
    ++start_pos;

  std::size_t const start_index = start_pos_to_index(start_pos);
  std::size_t const end_index = end_pos_to_index(end_pos, depths[0].size());
  uint64_t total_read_depth = 0;

  for (auto const & pn : sample_indexes)
  {
    assert(pn < depths.size());
    auto const & depth = depths[pn];

    if (depth.size() == 0)
    {
      // This happens only when this pn did not have any reads overlapping the graph,
      // so we can safely continue to the next one
      continue;
    }

    auto max_depth_it = std::max_element(depth.begin() + start_index, depth.begin() + end_index);
    assert(max_depth_it != depth.begin() + end_index);
    total_read_depth += *max_depth_it;
  }

  return total_read_depth;
}

void ReferenceDepth::add_genotype_paths(GenotypePaths const & geno, long const sample_index)
{
  assert(depths.size() > 0);

  if (sample_index >= static_cast<long>(depths.size()))
  {
    print_log(log_severity::error,
              __HERE__,
              " Odd.. sample_index >= expected (",
              sample_index,
              " >= ",
              depths.size(),
              ")");
    return;
    // std::exit(1);
  }

  assert(sample_index < static_cast<long>(depths.size()));

  if (geno.paths.size() == 0 || geno.paths[0].size() < 63)
    return;

  auto & depth = depths[sample_index];
  assert(depth.size() > 0);

  // Simple algorithm for the case there is only one path
  if (geno.paths.size() == 1)
  {
    auto const & path = geno.paths[0];
    long const start_pos = path.start_ref_reach_pos() - path.read_start_index;
    long const end_pos = path.end_ref_reach_pos() + (geno.read_length - 1 - path.read_end_index);

    long const start_index = start_pos_to_index(start_pos);
    auto const end = depth.begin() + end_pos_to_index(end_pos, depth.size());

    if (start_index < static_cast<long>(depth.size()))
    {
      for (auto it = depth.begin() + start_index; it != end && it != depth.end(); ++it)
        ++(*it); // Increase the depth by one
    }
  }
  else
  {
    std::unordered_set<long> local_depth;

    auto increase_local_depth_lambda = [&](long start_pos, long end_pos)
    {
      assert(end_pos >= start_pos);

      // Needed for SV breakpoint indel calling
      if (end_pos < reference_offset)
        return;

      assert(end_pos >= reference_offset);
      assert(depth.size() > 0);
      long const start_index = start_pos_to_index(start_pos);
      auto const end = depth.begin() + end_pos_to_index(end_pos, depth.size());

      // TODO: Super rarely this fails, find out why.
      // Must likely it is related to SV breakpoint indel calling
      if (start_index < static_cast<long>(depth.size()))
      {
        for (auto it = depth.begin() + start_index; it != end && it != depth.end(); ++it)
          local_depth.insert(std::distance(depth.begin(), it)); // Increase the depth by one
      }
    };

    for (auto const & path : geno.paths)
    {
      long const start_pos = path.start_ref_reach_pos() - path.read_start_index;
      long const end_pos = path.end_ref_reach_pos() + (geno.read_length - 1 - path.read_end_index);

      if (end_pos - start_pos >= 50)
        increase_local_depth_lambda(start_pos + 4, end_pos - 4);
      else
        increase_local_depth_lambda(start_pos, end_pos);
    }

    // Commit depth
    for (auto it = local_depth.begin(); it != local_depth.end(); ++it)
    {
      assert(*it < static_cast<long>(depth.size()));
      auto & pn_depth = depth[*it];

      // Check for overflow
      if (pn_depth < 0xFFFFul)
        ++pn_depth;
    }
  }
}

void ReferenceDepth::add_depth(long const start_pos, long const end_pos, long const sample_index)
{
  assert(sample_index < static_cast<long>(depths.size()));
  auto & depth = depths[sample_index];

  long const start_index = start_pos_to_index(start_pos);
  auto const end = depth.begin() + end_pos_to_index(end_pos, depth.size());
  assert(depth.size() > 0);

  if (start_index < static_cast<long>(depth.size()))
  {
    for (auto it = depth.begin() + start_index; it != end && it != depth.end(); ++it)
    {
      ++(*it); // Increase the depth by one
    }
  }
}

long ReferenceDepth::start_pos_to_index(long const start_pos) const
{
  return (start_pos < reference_offset) ? 0 : (start_pos - reference_offset);
}

long ReferenceDepth::end_pos_to_index(long const end_pos, long const depth_size) const
{
  return (end_pos > reference_offset + depth_size) ? depth_size : end_pos + 1 - reference_offset;
}

} // namespace gyper
