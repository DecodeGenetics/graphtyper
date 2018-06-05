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


} // anon namespace


std::ostream &
operator<<(std::ostream & os, std::vector<char> const & dt)
{
  for (auto it = dt.begin(); it != dt.end(); ++it)
  {
    os << *it;
  }
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
  :  pos(0u), ref(0), alts(0) {}

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
insert_prior_sequence(VarRecord & current, VarRecord const & previous)
{
  current.ref.insert(current.ref.begin(), previous.ref.begin(), previous.ref.begin() + current.pos - previous.pos);

  for (unsigned i = 0; i < current.alts.size(); ++i)
    current.alts[i].insert(current.alts[i].begin(), previous.ref.begin(), previous.ref.begin() + current.pos - previous.pos);

  current.pos = previous.pos;
}


void
extend_previous_records_alts(VarRecord const & current, VarRecord & previous)
{
  for (unsigned i = 0; i < previous.alts.size(); ++i)
    previous.alts[i].insert(previous.alts[i].end(), current.ref.end() - (current.ref.size() - previous.ref.size()), current.ref.end());
}


void
extend_current_record(VarRecord & current, VarRecord const & previous)
{
  for (unsigned i = 0; i < current.alts.size(); ++i)
    current.alts[i].insert(current.alts[i].end(), previous.ref.end() - (previous.ref.size() - current.ref.size()), previous.ref.end());

  current.ref.insert(current.ref.end(), previous.ref.end() - (previous.ref.size() - current.ref.size()), previous.ref.end());
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
  // std::cout << "Current record final = "; print_var_record(*this);
}


void
VarRecord::merge(VarRecord && prev_record)
{
  assert(prev_record.pos <= pos);

  // Check if the previous records ends where this one starts
  if (prev_record.pos + prev_record.ref.size() == pos)
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
        assert (prev_alt.size() > 0);
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
    std::size_t const opos = pos - prev_record.pos; // Original reference allele pos
    std::vector<char> const oref = ref; // Original reference allele

    // Insert new sequences prior to the reference sequence and all alternative sequences (if any)
    insert_prior_sequence(*this, prev_record);

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

    // std::cout << "Current record3 = "; print_var_record(*this);
    // std::cout << "Prev record3 = "; print_var_record(prev_record);

    assert(ref.size() == prev_record.ref.size());

    // Create new variants on the previous record if some this variant does not overlap with any of the previous ones
    std::size_t const original_size = prev_record.alts.size();

    for (unsigned i = 0; i < original_size; ++i)
    {
      if (prev_record.alts[i].size() < oref.size())
        continue;

      std::size_t const sd = opos;
      bool suffix_matches = true;
      assert(sd == opos);

      for (unsigned k = 0; k < oref.size(); ++k)
      {
        if (k + sd >= prev_record.alts[i].size() or oref[k] != prev_record.alts[i][k + sd])
        {
          suffix_matches = false;
          break;
        }
      }

      if (suffix_matches)
      {
        std::vector<char> prefix_alt(prev_record.alts[i].begin(), prev_record.alts[i].begin() + sd);

        for (unsigned j = 0; j < alts.size(); ++j)
        {
          std::vector<char> new_alt(prefix_alt);
          new_alt.insert(new_alt.end(), alts[j].cbegin() + sd, alts[j].cend());
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

  while(ref_it != ref.rend() && suffix_size < ref.size() - 1 &&
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
