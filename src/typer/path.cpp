#include <bitset> // std::bitset<Size>
#include <vector> // std::vector<Type>

#include <graphtyper/typer/path.hpp>
#include <graphtyper/graph/graph.hpp>


namespace gyper
{

Path::Path(KmerLabel const & l, uint16_t const _read_start_index, uint16_t const _read_end_index, uint16_t const _mismatches) noexcept
  : start(l.start_index)
  , end(l.end_index)
  , read_start_index(_read_start_index)
  , read_end_index(_read_end_index)
  , mismatches(_mismatches)
{
  assert(_read_end_index > _read_start_index);

  if (l.variant_id != INVALID_ID)
  {
    var_order.push_back(l.variant_order);
    std::bitset<MAX_NUMBER_OF_HAPLOTYPES> new_bitset(1);
    new_bitset <<= l.variant_num;
    nums.push_back(new_bitset);
  }
}


Path::Path(Path const & p1, Path const & p2) noexcept
{
  *this = p2; // Take everything from the latter path

  for (unsigned i = 0; i < p1.var_order.size(); ++i)
  {
    bool id_found = false;

    for (unsigned j = 0; j < var_order.size(); ++j)
    {
      if (p1.var_order[i] == var_order[j])
      {
        nums[j] &= p1.nums[i];

        if (nums[j].none())
          return;

        id_found = true;
        break;
      }
    }

    if (not id_found)
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

void
Path::merge_with_current(KmerLabel const & l)
{
  assert(l.end_index == end);
  assert(l.start_index == start);

  if (l.variant_order == INVALID_ID)
    return;

  for (unsigned i = 0; i < var_order.size(); ++i)
  {
    if (var_order[i] == l.variant_order)
    {
      std::bitset<MAX_NUMBER_OF_HAPLOTYPES> new_bitset(1);
      new_bitset <<= l.variant_num;
      nums[i] |= new_bitset;
      return;
    }
  }

  var_order.push_back(l.variant_order);
  std::bitset<MAX_NUMBER_OF_HAPLOTYPES> new_bitset(1);
  new_bitset <<= l.variant_num;
  nums.push_back(std::move(new_bitset));
}


/********************
 * PATH INFORMATION *
 ********************/

uint32_t
Path::start_pos() const
{
  return start;
}


uint32_t
Path::end_pos() const
{
  return end;
}


uint32_t
Path::start_correct_pos() const
{
  return graph.get_actual_pos(start);
}


uint32_t
Path::start_ref_reach_pos() const
{
  return graph.get_ref_reach_pos(start);
}


uint32_t
Path::end_correct_pos() const
{
  return graph.get_actual_pos(end);
}


uint32_t
Path::end_ref_reach_pos() const
{
  return graph.get_ref_reach_pos(end);
}


uint32_t
Path::size() const
{
  assert(read_end_index != read_start_index);
  return read_end_index - read_start_index + 1u;
}


uint32_t
Path::get_read_end_index(uint32_t const read_length) const
{
  return std::min(read_start_index + size() * (K - 1), read_length - 1);
}


bool
Path::is_reference() const
{
  for (auto const & num : nums)
  {
    if (not num.test(0))
      return false;
  }

  return true;
}


bool
Path::is_purely_reference() const
{
  for (auto const & num : nums)
  {
    if (not num.test(0) or num.count() > 1)
      return false;
  }

  return true;
}


bool
Path::is_empty() const
{
  return start == end;
}


} // namespace gyper
