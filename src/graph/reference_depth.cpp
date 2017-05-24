#include <algorithm> // std::sort
#include <string> // std::string
#include <sstream> // std::stringstream
#include <vector> // std::vector

#include <graphtyper/graph/graph.hpp>
#include <graphtyper/graph/reference_depth.hpp>
#include <graphtyper/typer/genotype_paths.hpp>
#include <graphtyper/typer/variant_candidate.hpp>


namespace gyper
{

ReferenceDepth::ReferenceDepth()
{
  resize_depth();
}


std::string
ReferenceDepth::print_non_zero_depth(std::size_t const MIN_DEPTH) const
{
  std::stringstream ss;

  for (std::size_t i = 0; i < depth.size(); ++i)
  {
    if (depth[i] >= MIN_DEPTH)
    {
      ss << (i + reference_offset) << "\t" << depth[i] << "\n";
    }
  }

  return ss.str();
}


uint16_t
ReferenceDepth::get_read_depth(VariantCandidate const & var) const
{
  assert(depth.size() > 0);
  assert(var.seqs.size() > 0);
  uint32_t start_pos = var.abs_pos;
  uint32_t end_pos = start_pos + var.seqs[0].size() - 1;

  // We want to avoid getting the coverage at the first pos if possible, since the first positions very often match in more than one path
  if (var.seqs[0].size() > 1)
    ++start_pos;

  std::size_t const start_index = start_pos_to_index(start_pos);
  std::size_t const end_index = end_pos_to_index(end_pos);

  auto max_depth_it = std::max_element(depth.begin() + start_index, depth.begin() + end_index);
  assert(max_depth_it != depth.begin() + end_index);
  return *max_depth_it;
}


void
ReferenceDepth::add_genotype_paths(GenotypePaths const & geno)
{
  if (geno.paths.size() == 0 || geno.paths[0].size() < 63)
    return;

  assert(local_depth.size() == 0);

  for (auto const & path : geno.paths)
  {
    std::size_t const start_pos = path.start_ref_reach_pos() - path.read_start_index;
    std::size_t const end_pos = path.end_ref_reach_pos() + (geno.read.size() - 1 - path.read_end_index);

    if (end_pos - start_pos >= 50)
      increase_local_depth_by_one(start_pos + 4, end_pos - 4);
    else
      increase_local_depth_by_one(start_pos, end_pos);
  }

  commit_local_depth();
}


void
ReferenceDepth::add_reference_depths_from(ReferenceDepth const & ref_depth)
{
  if (ref_depth.depth.size() == 0)
    return;

  assert(ref_depth.reference_offset == this->reference_offset);
  assert(ref_depth.depth.size() == this->depth.size());

  for (std::size_t i = 0; i < depth.size(); ++i)
    this->depth[i] += ref_depth.depth[i];
}


void
ReferenceDepth::clear()
{
  depth.clear();
  depth.shrink_to_fit();
  local_depth.clear();
  local_depth.shrink_to_fit();
}


std::size_t
ReferenceDepth::start_pos_to_index(uint32_t const start_pos) const
{
  return (start_pos < reference_offset) ? 0 : (start_pos - reference_offset);
}


std::size_t
ReferenceDepth::end_pos_to_index(uint32_t const end_pos) const
{
  return (end_pos > reference_offset + depth.size()) ? depth.size() : end_pos + 1 - reference_offset;
}


void
ReferenceDepth::increase_local_depth_by_one(std::size_t const start_pos, std::size_t const end_pos)
{
  assert(end_pos >= start_pos);
  assert(end_pos >= reference_offset);
  assert(depth.size() > 0);
  std::size_t const start_index = start_pos_to_index(start_pos);
  std::size_t const end_index = end_pos_to_index(end_pos);
  assert(start_index < depth.size());

  for (auto it = depth.begin() + start_index; it != depth.begin() + end_index && it != depth.end(); ++it)
    local_depth.push_back(std::distance(depth.begin(), it)); // Increase the depth by one
}


void
ReferenceDepth::commit_local_depth()
{
  std::sort(local_depth.begin(), local_depth.end());
  auto last = std::unique(local_depth.begin(), local_depth.end());

  for (auto it = local_depth.begin(); it != last; ++it)
  {
    assert(*it < depth.size());

    // Make sure we will not overflow
    if (depth[*it] < 0xFFFFul)
      ++depth[*it];
  }

  local_depth.clear();
}


void
ReferenceDepth::resize_depth()
{
  depth.resize(graph.reference.size(), 0u);
  reference_offset = graph.ref_nodes.size() > 0 ? graph.ref_nodes[0].get_label().order : 0;
}


void
GlobalReferenceDepth::set_pn_count(std::size_t const pn_count)
{
  reference_depth_mutexes = std::vector<std::mutex>(pn_count);
  depths.resize(pn_count);
}


void
GlobalReferenceDepth::add_reference_depths_from(ReferenceDepth const & ref_depth, std::size_t const pn_index)
{
  if (ref_depth.depth.size() == 0)
    return;

  if (reference_offset == 0)
    reference_offset = ref_depth.reference_offset;

  assert(reference_offset == ref_depth.reference_offset);
  assert(pn_index < reference_depth_mutexes.size());
  std::lock_guard<std::mutex> lock(reference_depth_mutexes[pn_index]);

  assert(pn_index < this->depths.size());

  if (depths[pn_index].size() == 0)
    depths[pn_index].resize(ref_depth.depth.size());

  assert(ref_depth.depth.size() == depths[pn_index].size());

  for (std::size_t i = 0; i < depths[pn_index].size(); ++i)
  {
    // Check for overflow
    if (static_cast<uint32_t>(depths[pn_index][i]) + static_cast<int32_t>(ref_depth.depth[i]) < 0xFFFFul)
      depths[pn_index][i] += ref_depth.depth[i];
    else
      depths[pn_index][i] = 0xFFFFul;
  }
}


uint16_t
GlobalReferenceDepth::get_read_depth(VariantCandidate const & var, std::size_t const pn_index) const
{
  assert(pn_index < reference_depth_mutexes.size());
  std::lock_guard<std::mutex> lock(reference_depth_mutexes[pn_index]); // Create lock just to be sure, we should in general not need it
  auto & depth = depths[pn_index];
  assert(depth.size() > 0);
  assert(var.seqs.size() > 0);
  uint32_t start_pos = var.abs_pos;
  uint32_t end_pos = start_pos + var.seqs[0].size() - 1;

  // We want to avoid getting the coverage at the first pos if possible, since the first positions very often match in more than one path
  if (var.seqs[0].size() > 1)
    ++start_pos;

  std::size_t const start_index = start_pos_to_index(start_pos);
  std::size_t const end_index = end_pos_to_index(end_pos, depth.size());

  auto max_depth_it = std::max_element(depth.begin() + start_index, depth.begin() + end_index);
  assert(max_depth_it != depth.begin() + end_index);
  return *max_depth_it;
}


uint64_t
GlobalReferenceDepth::get_total_read_depth_of_samples(VariantCandidate const & var, std::vector<uint32_t> const & pn_indexes) const
{
  assert(var.seqs.size() > 0);
  uint32_t start_pos = var.abs_pos;
  uint32_t end_pos = start_pos + var.seqs[0].size() - 1;

  // We want to avoid getting the coverage at the first pos if possible, since the first positions very often match in more than one path
  if (var.seqs[0].size() > 1)
    ++start_pos;

  std::size_t const start_index = start_pos_to_index(start_pos);
  std::size_t const end_index = end_pos_to_index(end_pos, depths[0].size());
  uint64_t total_read_depth = 0;

  for (auto const & pn : pn_indexes)
  {
    assert(pn < depths.size());

    if (depths[pn].size() == 0)
      continue; // This happens only when this pn did not have any reads overlapping the graph, so we can safely continue to the next one

    auto max_depth_it = std::max_element(depths[pn].begin() + start_index, depths[pn].begin() + end_index);
    assert(max_depth_it != depths[pn].begin() + end_index);
    total_read_depth += *max_depth_it;
  }

  return total_read_depth;
}


std::size_t
GlobalReferenceDepth::start_pos_to_index(uint32_t const start_pos) const
{
  return (start_pos < reference_offset) ? 0 : (start_pos - reference_offset);
}


std::size_t
GlobalReferenceDepth::end_pos_to_index(uint32_t const end_pos, std::size_t const depth_size) const
{
  return (end_pos > reference_offset + depth_size) ? depth_size : end_pos + 1 - reference_offset;
}


GlobalReferenceDepth global_reference_depth;

} // namespace gyper
