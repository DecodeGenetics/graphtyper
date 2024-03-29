#include <algorithm>
#include <cassert> // assert
#include <set>
#include <vector> // std::vector<Type>

#include <parallel_hashmap/phmap.h>

#include <graphtyper/graph/graph.hpp>
#include <graphtyper/typer/path.hpp>

namespace gyper
{
Path::Path(Graph const & graph,
           KmerLabel const & l,
           uint16_t const _read_start_index,
           uint16_t const _read_end_index,
           uint16_t const _mismatches) noexcept :
  start(l.start_index),
  end(l.end_index),
  read_start_index(_read_start_index),
  read_end_index(_read_end_index),
  mismatches(_mismatches)
{
  assert(_read_end_index > _read_start_index);

  if (l.variant_id != gyper::INVALID_ID)
  {
    assert(l.variant_id < graph.var_nodes.size());
    // var_order.push_back(graph.var_nodes[l.variant_id].get_label().order);
    var_order.push_back(graph.get_variant_order(l.variant_id));

    phmap::flat_hash_set<uint16_t> new_num;
    new_num.insert(graph.get_variant_num(l.variant_id));
    nums.push_back(std::move(new_num));
  }
}

Path::Path(Path const & p1, Path const & p2) : Path(p2) // Take everything from the latter path
{
  for (long i = 0; i < static_cast<long>(p1.var_order.size()); ++i)
  {
    bool id_found = false;

    for (long j = 0; j < static_cast<long>(var_order.size()); ++j)
    {
      if (p1.var_order[i] == var_order[j])
      {
#if 0 // first implementation seems faster; copies more but has better asymptotic complexity
        nums[j].clear();
        std::set_intersection(p1.nums[i].begin(), p1.nums[i].end(),
                              p2.nums[j].begin(), p2.nums[j].end(),
                              std::insert_iterator<std::set<uint64_t> >{nums[j], nums[j].end()});

#else
        for (auto it = nums[j].begin(), e = nums[j].end(); it != e; /**/)
        {
          if (p1.nums[i].count(*it) == 0)
            it = nums[j].erase(it);
          else
            ++it;
        }
#endif

        if (nums[j].empty())
          return;

        id_found = true;
        break;
      }
    }

    if (!id_found)
    {
      var_order.push_back(p1.var_order[i]);
      nums.push_back(p1.nums[i]);
    }
  }

  read_start_index = p1.read_start_index;
  start = p1.start;
  mismatches += p1.mismatches;
}

/**********************
 * PATH MODIFICATIONS *
 **********************/

void Path::erase_var_order(long const index)
{
  assert(index < static_cast<long>(var_order.size()));
  assert(index < static_cast<long>(nums.size()));
  var_order.erase(var_order.begin() + index);
  nums.erase(nums.begin() + index);
}

void Path::erase_ref_support(long const index)
{
  assert(index < static_cast<long>(var_order.size()));
  assert(index < static_cast<long>(nums.size()));

  if (nums[index].count(0))
    erase_var_order(index); // delete if the sequence supports the reference
}

void Path::merge_with_current(KmerLabel const & l)
{
  assert(l.end_index == end);
  assert(l.start_index == start);

  if (l.variant_id == INVALID_ID)
    return;

  auto const variant_order = graph.get_variant_order(l.variant_id);
  auto const variant_num = graph.get_variant_num(l.variant_id);

  for (int i = 0; i < static_cast<int>(var_order.size()); ++i)
  {
    if (var_order[i] == variant_order)
    {
      nums[i].insert(variant_num);
      return;
    }
  }

  var_order.push_back(variant_order);
  phmap::flat_hash_set<uint16_t> new_num;
  new_num.insert(variant_num);
  nums.push_back(std::move(new_num));
}

/********************
 * PATH INFORMATION *
 ********************/

uint32_t Path::start_pos() const
{
  return start;
}

uint32_t Path::end_pos() const
{
  return end;
}

uint32_t Path::start_correct_pos() const
{
  return graph.get_actual_pos(start);
}

uint32_t Path::start_ref_reach_pos() const
{
  return graph.get_ref_reach_pos(start);
}

uint32_t Path::end_correct_pos() const
{
  return graph.get_actual_pos(end);
}

uint32_t Path::end_ref_reach_pos() const
{
  return graph.get_ref_reach_pos(end);
}

uint32_t Path::size() const
{
  assert(read_end_index != read_start_index);
  return read_end_index - read_start_index + 1u;
}

uint32_t Path::get_read_end_index(uint32_t const read_length) const
{
  return std::min(read_start_index + size() * (K - 1), read_length - 1);
}

bool Path::is_reference() const
{
  for (auto const & num : nums)
  {
    if (num.count(0) == 0)
      return false;
  }

  return true;
}

bool Path::is_purely_reference() const
{
  for (auto const & num : nums)
  {
    if (num.count(0) == 0 || num.size() > 1)
      return false;
  }

  return true;
}

bool Path::is_empty() const
{
  return start == end;
}

} // namespace gyper
