#include <algorithm>
#include <cassert>
#include <iostream>
#include <iterator>
#include <string>
#include <sstream>

#include <boost/log/trivial.hpp>

#include <graphtyper/graph/var_record.hpp>


namespace
{

bool
allele_uniqueness(std::vector<char> const & lhs, std::vector<char> const & rhs)
{
  return lhs == rhs;
}


bool
compare_two_alleles(std::vector<char> const & lhs, std::vector<char> const & rhs)
{
  return lhs < rhs;
}


void
insert_prior_sequence(gyper::VarRecord & current, gyper::VarRecord const & previous)
{
  current.ref.insert(current.ref.begin(), previous.ref.begin(), previous.ref.begin() + current.pos - previous.pos);

  for (long i = 0; i < static_cast<long>(current.alts.size()); ++i)
  {
    current.alts[i].insert(current.alts[i].begin(),
                           previous.ref.begin(),
                           previous.ref.begin() + current.pos - previous.pos);
  }

  current.pos = previous.pos;
}


void
extend_previous_records_alts(gyper::VarRecord const & current, gyper::VarRecord & previous)
{
  for (unsigned i = 0; i < previous.alts.size(); ++i)
    previous.alts[i].insert(previous.alts[i].end(),
                            current.ref.end() - (current.ref.size() - previous.ref.size()),
                            current.ref.end());
}


void
extend_current_record(gyper::VarRecord & current, gyper::VarRecord const & previous)
{
  for (unsigned i = 0; i < current.alts.size(); ++i)
    current.alts[i].insert(current.alts[i].end(),
                           previous.ref.end() - (previous.ref.size() - current.ref.size()),
                           previous.ref.end());

  current.ref.insert(current.ref.end(),
                     previous.ref.end() - (previous.ref.size() - current.ref.size()),
                     previous.ref.end());
}


} // anon namespace


std::ostream &
operator<<(std::ostream & os, std::vector<char> const & dt)
{
  for (auto it = dt.begin(); it != dt.end(); ++it)
    os << *it;

  return os;
}


std::string
to_string(std::vector<char> const & dt)
{
  std::stringstream ss;
  ss << dt;
  return ss.str();
}


namespace gyper
{

VarRecord::VarRecord()
  :  pos(0u), ref(0), alts(0)
{}

VarRecord::VarRecord(uint32_t const p)
  : pos(p)
{}

VarRecord::VarRecord(uint32_t const p, std::vector<char> && r, std::vector<std::vector<char> > && a)
  : pos(std::move(p))
  , ref(std::move(r))
  , alts(std::move(a))
{}

void
VarRecord::clear()
{
  pos = 0;
  ref.clear();
  alts.clear();
}


std::string
VarRecord::to_string() const
{
  std::ostringstream ss;
  ss << ref << " @ " << pos << ":";

  for (auto const & alt : alts)
    ss << " " << alt;

  return ss.str();
}


void
VarRecord::merge_one_path(VarRecord && prev_record)
{
  assert(prev_record.pos <= pos);
  std::vector<char> const oref = ref; // Original reference allele

  // Insert new sequences prior to the reference sequence and all alternative sequences (if any)
  insert_prior_sequence(*this, prev_record);

  // Check if the previous reference reacher further than then current one
  if (ref.size() < prev_record.ref.size())
  {
    // In this case, we need to lengthen the current reference and all its variants
    extend_current_record(*this, prev_record);
  }
  else if (ref.size() > prev_record.ref.size())
  {
    // If not, then we need to extend the variants of the previous record
    extend_previous_records_alts(*this, prev_record);
    prev_record.ref = ref;
  }

  assert(ref.size() == prev_record.ref.size());
  std::move(prev_record.alts.begin(), prev_record.alts.end(), std::back_inserter(alts));
  std::sort(alts.begin(), alts.end(), compare_two_alleles);
  alts.erase(std::unique(alts.begin(), alts.end(), allele_uniqueness), alts.end());
}


void
VarRecord::merge(VarRecord && prev_record, long const EXTRA_SUFFIX)
{
  assert(prev_record.pos <= pos);

  // Check if the previous records ends where this one starts
  if (prev_record.pos + prev_record.ref.size() == pos)
  {
    if (EXTRA_SUFFIX == 0)
    {
      // std::cout << "SPECIAL CASE: WE START AND AT THE SAME LOCATION" << std::endl;
      // std::cout << "Current record = "; print_var_record(*this);
      // std::cout << "Prev record = "; print_var_record(prev_record);

      // R A,B
      // S C,D,E
      // =>
      // RS RC,RD,RE,AC,AD,AE,BC,BD,BE
      std::size_t const ORIGINAL_NUMBER_OF_CURR_ALTS = alts.size();

      // Gather new alts here and add them afterwards
      std::vector<std::vector<char> > new_alts;

      // Add AC,AD,AE,BC,BD,BE
      for (auto const & prev_alt : prev_record.alts)
      {
        for (auto const & curr_alt : alts)
        {
          assert(curr_alt.size() > 0);
          assert(prev_alt.size() > 0);
          std::vector<char> seq(prev_alt);
          std::copy(curr_alt.begin(), curr_alt.end(), std::back_inserter(seq));
          assert(seq.size() > 0);
          new_alts.push_back(std::move(seq));
        }

        std::vector<char> seq(prev_alt);
        std::copy(ref.begin(), ref.end(), std::back_inserter(seq));
        assert(seq.size() > 0);
        new_alts.push_back(std::move(seq));
      }

      std::move(new_alts.begin(), new_alts.end(), std::back_inserter(alts));

      // Update C,D,E to RC,RD,RE
      for (std::size_t i = 0; i < ORIGINAL_NUMBER_OF_CURR_ALTS; ++i)
      {
        alts[i].insert(alts[i].begin(), prev_record.ref.begin(), prev_record.ref.end());
      }

      // std::cout << "Current record3 = "; print_var_record(*this);
      // std::cout << "Prev record3 = "; print_var_record(prev_record);

      // Extend the reference to RS
      ref.insert(ref.begin(), prev_record.ref.begin(), prev_record.ref.end());
      // std::cout << "Prev record4 = "; print_var_record(prev_record);
      pos = prev_record.pos;
    }
    else
    {
      merge_one_path(std::move(prev_record));
    }
  }
  else
  {
    long const jump_size = pos - prev_record.pos; // Original reference allele pos
    long const oref_size = ref.size(); // Original reference allele size

    // Insert new sequences prior to the reference sequence and all alternative sequences (if any)
    insert_prior_sequence(*this, prev_record);
    long const oref_size_pre = ref.size(); // reference allele after adding prefix

    // std::cout << "Current record2 = "; print_var_record(*this);
    // std::cout << "Prev record2 = "; print_var_record(prev_record);

    // Check if the previous reference reacher further than then current one
    if (ref.size() < prev_record.ref.size())
    {
      // In this case, we need to lengthen the current reference and all its variants
      extend_current_record(*this, prev_record);
    }
    else if (ref.size() > prev_record.ref.size())
    {
      // If not, then we need to extend the variants of the previous record
      extend_previous_records_alts(*this, prev_record);
      prev_record.ref = ref;
    }

    long const extension_size = ref.size() - oref_size_pre;
    //std::cerr << "Current record3 = " << to_string() << "\n";
    //std::cerr << "Prev record3 = " << prev_record.to_string() << "\n";

    assert(ref.size() == prev_record.ref.size());

    // Create new variants on the previous record if some this variant does not overlap with any of the previous ones
    long const original_size = prev_record.alts.size();

    for (long i = 0; i < original_size; ++i)
    {
      auto & alt = prev_record.alts[i];

      if (static_cast<long>(alt.size()) <= oref_size)
        continue;

      long suffix_matches = 0;
      long const offset = ref.size() - alt.size();
      long const smaller_allele_size = std::min(ref.size(), alt.size());

      //std::cerr << oref_size << " " << alt.size() << "  " << ref.size() << " " << offset << "\n";
      //std::cerr << ref << " " << alt << "\n";

      for (long k = 0; k < smaller_allele_size; ++k)
      {
        if (*(ref.rbegin() + k) == *(alt.rbegin() + k))
        {
          ++suffix_matches;
        }
        else
        {
          break;
        }
      }

      //std::cerr << "suffix matches, jump_size, oref_size, es, extension_size = " << suffix_matches << ","
      //          << jump_size << "," << oref_size << "," << EXTRA_SUFFIX << "," << extension_size << "\n";

      //long const sd = jump_size;
      //bool suffix_matches = true;

      /*for (long k = 0; k < static_cast<long>(oref.size()); ++k)
      {
        if (k + sd >= static_cast<long>(prev_record.alts[i].size()) || oref[k] != *(alt[i].rbegin() + k + sd))
        {
          suffix_matches = false;
          break;
        }
      }*/

      if (suffix_matches >= extension_size + oref_size + EXTRA_SUFFIX) //(oref_size - 1 + jump_size + EXTRA_SUFFIX))
      {
        std::vector<char> prefix_alt(alt.begin(), alt.begin() + jump_size - offset);
        //std::cerr << "prefix alt = " << std::string(prefix_alt.begin(), prefix_alt.end()) << "\n";

        for (long j = 0; j < static_cast<long>(alts.size()); ++j)
        {
          std::vector<char> new_alt(prefix_alt);
          new_alt.insert(new_alt.end(), alts[j].cbegin() + jump_size, alts[j].cend());
          prev_record.alts.push_back(std::move(new_alt));
        }
      }
    }

    std::move(prev_record.alts.begin(), prev_record.alts.end(), std::back_inserter(alts));
  }

  std::sort(alts.begin(), alts.end(), compare_two_alleles);
  alts.erase(std::unique(alts.begin(), alts.end(), allele_uniqueness), alts.end());
}


std::vector<char>
VarRecord::get_common_suffix()
{
  uint64_t suffix_size = 0;
  auto ref_it = ref.rbegin();

  while (ref_it != ref.rend() && suffix_size < ref.size() - 1 &&
         std::all_of(alts.begin(),
                     alts.end(),
                     [&](std::vector<char> const & alt)
    {
      return suffix_size < (alt.size() - 1) && *(alt.rbegin() + suffix_size) == *ref_it;
    }
                     )
         )
  {
    ++ref_it;
    ++suffix_size;
  }

  return std::vector<char>(ref.end() - suffix_size, ref.end());
}


} // namespace gyper
