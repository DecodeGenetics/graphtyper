#pragma once

#include <cstdint> // uint32_t
#include <iterator>
#include <vector> // std::vector

uint32_t inline
count_mismatches(std::vector<char> const & read,
                 uint32_t read_offset,
                 std::vector<char> const & graph_dna,
                 uint32_t dna_index,
                 uint32_t max_mismatches
                 )
{
  uint32_t mismatches = 0;
  auto read_it = read.begin() + read_offset;
  auto graph_it = graph_dna.begin() + dna_index;
  // If we find a '<' we allow anything to align there without mismatches

  while (graph_it != graph_dna.end() and read_it != read.end())
  {
    if (*graph_it == '>' || *graph_it == '<')
    {
      return max_mismatches + 1; // Do not allow paths with >
    }
    else if (*graph_it != *read_it && *read_it != 'N' && *graph_it != 'N')
    {
      ++mismatches;

      if (mismatches > max_mismatches)
        return mismatches; // Stop at this point
    }

    ++read_it;
    ++graph_it;
  }

  return mismatches;
}


uint32_t inline
count_mismatches_backward(std::vector<char> const & read,
                          uint32_t read_offset,
                          std::vector<char> const & graph_dna,
                          uint32_t dna_index,
                          uint32_t max_mismatches
                          )
{
  uint32_t mismatches = 0;
  auto read_it = read.rbegin() + read_offset;
  auto graph_it = graph_dna.rbegin() + dna_index;
  // If we find a '>' we allow anything to align there without mismatches

  while (graph_it != graph_dna.rend() && read_it != read.rend())
  {
    if (*graph_it == '>' || *graph_it == '<')
    {
      return max_mismatches + 1; // Do not allow paths with <
    }
    else if (*graph_it != *read_it && *read_it != 'N' && *graph_it != 'N')
    {
      ++mismatches;

      if (mismatches > max_mismatches)
        return mismatches; // Stop at this point
    }

    ++read_it;
    ++graph_it;
  }

  return mismatches;
}


void inline
add_node_dna_to_sequence(std::vector<char> & seq, gyper::Label const & label, uint32_t const from, uint32_t const to)
{
  uint32_t const dna_end = label.order + label.dna.size();

  if (to <= from or dna_end <= from or label.order >= to)
    return;

  auto it1 = label.dna.begin();
  auto it2 = label.dna.end();

  if (label.order >= from)
  {
    if (to < dna_end)
      it2 -= dna_end - to;
  }
  else if (dna_end >= to)
  {
    it1 += from - label.order;
    it2 -= dna_end - to;
  }
  else
  {
    it1 += from - label.order;
  }

  std::copy(it1, it2, std::back_inserter(seq));
}
