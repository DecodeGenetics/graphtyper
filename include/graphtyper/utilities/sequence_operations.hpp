#pragma once

#include <assert.h>
#include <vector>


namespace gyper
{


void inline
remove_common_prefix(uint32_t & pos, std::vector<char> & ref, std::vector<std::vector<char> > & alts, bool keep_one_match = false)
{
  while (ref.size() > 1)
  {
    for (auto alt_it = alts.begin(); alt_it != alts.end(); ++alt_it)
    {
      assert(alt_it->size() > 0);
      if (alt_it->size() == 1 || (*alt_it)[0] != ref[0] || (keep_one_match && (*alt_it)[1] != ref[1])) // Never remove too much
        return;
    }

    // Remove one base from front
    ++pos;   // Increase position of the haplotypes

    // Erase first base from ref
    ref.erase(ref.begin());

    // Erase first base from alts
    for (auto alt_erase_it = alts.begin(); alt_erase_it != alts.end(); ++alt_erase_it)
    {
      alt_erase_it->erase(alt_erase_it->begin());
    }
  }
}


void inline
remove_common_prefix(uint32_t & pos, std::vector<std::vector<char> > & seqs, bool keep_one_match = false)
{
  if (seqs.size() <= 1 or seqs[0].size() == 0)
    return;

  std::vector<char> ref = seqs[0];
  std::vector<std::vector<char> > alts(seqs.begin() + 1, seqs.end());
  remove_common_prefix(pos, ref, alts, keep_one_match);
  seqs[0] = ref;

  for (std::size_t a = 1; a < seqs.size(); ++a)
  {
    seqs[a] = alts[a - 1];
  }
}


void inline
remove_common_suffix(std::vector<char> & ref, std::vector<std::vector<char> > & alts)
{
  while (ref.size() > 1)
  {
    for (auto alt_it = alts.begin(); alt_it != alts.end(); ++alt_it)
    {
      assert(alt_it->size() > 0);

      if (alt_it->size() == 1 or alt_it->back() != ref.back())
        return;
    }

    // Remove one base from back of reference
    ref.pop_back();

    for (auto alt_erase_it = alts.begin(); alt_erase_it != alts.end(); ++alt_erase_it)
      alt_erase_it->pop_back();
  }
}


void inline
remove_common_suffix(std::vector<std::vector<char> > & seqs)
{
  if (seqs.size() <= 1 or seqs[0].size() == 0)
    return;

  std::vector<char> ref = seqs[0];
  std::vector<std::vector<char> > alts(seqs.begin() + 1, seqs.end());
  remove_common_suffix(ref, alts);

  seqs[0] = ref;

  for (std::size_t a = 1; a < seqs.size(); ++a)
  {
    seqs[a] = alts[a - 1];
  }
}


} // namespace gyper
