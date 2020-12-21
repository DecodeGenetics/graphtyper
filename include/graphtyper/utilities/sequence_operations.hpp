#pragma once

#include <cassert>
#include <cstdint>
#include <vector>


namespace gyper
{


void inline
remove_common_prefix(long & pos,
                     std::vector<char> & ref,
                     std::vector<std::vector<char> > & alts,
                     bool keep_one_match = false)
{
  while (ref.size() > 1)
  {
    for (auto const & alt : alts)
    {
      assert(alt.size() > 0);

      if (alt.size() <= 1 || alt[0] != ref[0] || (keep_one_match && alt[1] != ref[1])) // Never remove too much
        return;
    }

    // Remove one base from front
    ++pos; // Increment position for each removed base

    // Erase first base from ref
    ref.erase(ref.begin());

    // Erase first base from each alternative allele
    for (auto & alt : alts) //alt_erase_it = alts.begin(); alt_erase_it != alts.end(); ++alt_erase_it)
      alt.erase(alt.begin());
  }
}


void inline
remove_common_prefix(long & pos, std::vector<std::vector<char> > & seqs, bool keep_one_match = false)
{
  if (seqs.size() <= 1 || seqs[0].size() <= 1)
    return;

  std::vector<char> ref = seqs[0];
  std::vector<std::vector<char> > alts(seqs.begin() + 1, seqs.end());
  remove_common_prefix(pos, ref, alts, keep_one_match);
  seqs[0] = ref;

  for (long a = 1; a < static_cast<long>(seqs.size()); ++a)
    seqs[a] = alts[a - 1];
}


void inline
remove_common_suffix(std::vector<char> & ref, std::vector<std::vector<char> > & alts)
{
  while (ref.size() > 1)
  {
    for (auto const & alt : alts)
    {
      assert(alt.size() > 0);

      if (alt.size() <= 1 || alt[alt.size() - 1] != ref[ref.size() - 1])
        return;
    }

    // Remove one base from back of reference
    ref.pop_back();

    for (auto & alt : alts)
      alt.pop_back();
  }
}


void inline
remove_common_suffix(std::vector<std::vector<char> > & seqs)
{
  if (seqs.size() <= 1 or seqs[0].size() <= 1)
    return;

  std::vector<char> ref = seqs[0];
  std::vector<std::vector<char> > alts(seqs.begin() + 1, seqs.end());
  remove_common_suffix(ref, alts);
  seqs[0] = ref;

  for (long a = 1; a < static_cast<long>(seqs.size()); ++a)
    seqs[a] = alts[a - 1];
}


} // namespace gyper
