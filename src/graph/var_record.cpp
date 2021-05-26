#include <algorithm>
#include <cassert>
#include <iostream>
#include <iterator>
#include <string>
#include <sstream>
#include <vector>

#include <graphtyper/utilities/logging.hpp>

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/var_record.hpp>


std::ostream &
operator<<(std::ostream & os, std::vector<char> const & dt)
{
  for (auto it = dt.begin(); it != dt.end(); ++it)
    os << *it;

  return os;
}


namespace
{

template <typename T, typename U>
void
copy_events_and_anti_events(T & out, U const & in)
{
  std::copy(in.events.begin(),
            in.events.end(),
            std::inserter(out.events, out.events.begin()));

  std::copy(in.anti_events.begin(),
            in.anti_events.end(),
            std::inserter(out.anti_events, out.anti_events.begin()));
}


void
insert_prior_sequence(gyper::VarRecord & current, gyper::VarRecord const & previous)
{
  assert(current.pos > previous.pos);
  current.ref.seq.insert(current.ref.seq.begin(),
                         previous.ref.seq.begin(),
                         previous.ref.seq.begin() + current.pos - previous.pos);

  for (long i = 0; i < static_cast<long>(current.alts.size()); ++i)
  {
    current.alts[i].seq.insert(current.alts[i].seq.begin(),
                               previous.ref.seq.begin(),
                               previous.ref.seq.begin() + current.pos - previous.pos);
  }

  current.pos = previous.pos;
}


void
extend_record(gyper::VarRecord & current, gyper::VarRecord const & previous)
{
  assert(previous.ref.seq.size() > current.ref.seq.size());

  long const size_diff = previous.ref.seq.size() - current.ref.seq.size();
  assert(size_diff < static_cast<long>(previous.ref.seq.size()));

  auto begin_it = previous.ref.seq.end() - size_diff;
  auto end_it = previous.ref.seq.end();

  for (auto & alt : current.alts)
    std::copy(begin_it, end_it, std::back_inserter(alt.seq));

  std::copy(begin_it, end_it, std::back_inserter(current.ref.seq));
}


void
extend_smaller_record(gyper::VarRecord & current, gyper::VarRecord & previous)
{
  // Check if the previous reference reacher further than then current one
  if (current.ref.seq.size() < previous.ref.seq.size())
  {
    // In this case, we need to lengthen the current reference and all its variants
    extend_record(current, previous);
  }
  else if (current.ref.seq.size() > previous.ref.seq.size())
  {
    // If not, then we need to extend the variants of the previous record
    extend_record(previous, current);
  }
}


void
move_alts(gyper::VarRecord & current, gyper::VarRecord && prev_record)
{
  long const alts_size_original = current.alts.size();

  // Check if any of the inserted alts are already there
  for (long pa{0}; pa < static_cast<long>(prev_record.alts.size()); ++pa)
  {
    auto && prev_alt = prev_record.alts[pa];
    bool is_good{true};

    for (long a{0}; a < alts_size_original; ++a)
    {
      if (current.alts[a].seq == prev_alt.seq)
      {
        is_good = false; // is not unique
        break;
      }
    }

    if (is_good)
      current.alts.push_back(std::move(prev_alt));
  }
}


} // anon namespace


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
  :  pos(0u)
  , ref()
  , alts(0)
  , is_sv(false)
{}

VarRecord::VarRecord(uint32_t const p)
  : pos(p)
  , ref()
  , alts(0)
  , is_sv(false)
{}

VarRecord::VarRecord(uint32_t const p, Ref && r, std::vector<Alt> && a)
  : pos(std::move(p))
  , ref(std::move(r))
  , alts(std::move(a))
  , is_sv(false)
{}


VarRecord::VarRecord(VarRecord const & o)
  : pos(o.pos)
  , ref(o.ref)
  , alts(o.alts)
  , is_sv(o.is_sv)
{}


VarRecord::VarRecord(VarRecord && o)
  : pos(o.pos)
  , ref(std::forward<Ref>(o.ref))
  , alts(std::forward<std::vector<Alt> >(o.alts))
  , is_sv(o.is_sv)
{}


VarRecord &
VarRecord::operator=(VarRecord const & o)
{
  pos = o.pos;
  ref = o.ref;
  alts = o.alts;
  is_sv = o.is_sv;
  return *this;
}


VarRecord &
VarRecord::operator=(VarRecord && o)
{
  pos = o.pos;
  ref = std::move(o.ref);
  alts = std::move(o.alts);
  is_sv = o.is_sv;
  return *this;
}


void
VarRecord::clear()
{
  pos = 0;
  ref.clear();
  alts.clear();
  ref.events.clear();
  ref.anti_events.clear();

  for (Alt & alt : alts)
    alt.clear();

  is_sv = false;
}


std::string
VarRecord::to_string() const
{
  std::ostringstream ss;
  ss << "pos,ref:alts=" << pos << ", " << ref.seq << ":";

  for (auto const & alt : alts)
    ss << " " << alt.seq;

  return ss.str();
}


void
VarRecord::merge_one_path(VarRecord && prev_record)
{
  assert(pos >= prev_record.pos);

  // Insert new sequences prior to the reference sequence and all alternative sequences (if any)
  if (prev_record.pos < pos)
    insert_prior_sequence(*this, prev_record);

  extend_smaller_record(*this, prev_record);
  assert(ref.seq == prev_record.ref.seq); // the records have the same ref now

  // copy events
  copy_events_and_anti_events(ref, prev_record.ref);

  // add prev ref events to current alts
  for (Alt & alt : alts)
    copy_events_and_anti_events(alt, prev_record.ref);

  move_alts(*this, std::move(prev_record));
  // TODO remove this as it messes up haplotype things
  //std::sort(alts.begin(), alts.end(), compare_two_alts);
}


void
VarRecord::merge_all(VarRecord && prev_record)
{
  // Check if the previous records ends where this one starts
  assert(prev_record.pos + prev_record.ref.seq.size() >= pos);

  if (prev_record.pos + prev_record.ref.seq.size() == pos)
  {
    print_log(log_severity::debug, __HERE__, " merging all records");
    // R A,B
    // S C,D,E
    // =>
    // RS RC,RD,RE,AC,AD,AE,BC,BD,BE
    //long const ORIGINAL_NUMBER_OF_CURR_ALTS = alts.size();

    // Gather new alts here and add them afterwards
    VarRecord new_record(prev_record.pos); // its ok that ref is just set to empty string - it is not used

    //std::vector<Alt> new_alts;

    // Add AC,AD,AE,BC,BD,BE
    for (auto const & prev_alt : prev_record.alts)
    {
      for (auto const & curr_alt : alts)
      {
        assert(curr_alt.seq.size() > 0);
        assert(prev_alt.seq.size() > 0);

        if (is_ok_to_merge_alts(prev_alt, curr_alt))
          new_record.alts.push_back(make_alt(prev_alt, curr_alt, 0));
      }

      // Add reference events from current ref to a new alt
      // AS, BS
      Alt new_alt(prev_alt);
      std::copy(ref.seq.begin(), ref.seq.end(), std::back_inserter(new_alt.seq));
      copy_events_and_anti_events(new_alt, ref);
      new_record.alts.push_back(std::move(new_alt));
    }

    // Update C,D,E to RC,RD,RE
    for (long i{0}; i < static_cast<long>(alts.size()); ++i)
    {
      auto & alt = alts[i];
      alt.seq.insert(alt.seq.begin(), prev_record.ref.seq.begin(), prev_record.ref.seq.end());
      copy_events_and_anti_events(alt, prev_record.ref);
    }

    // Extend the reference to RS
    pos = prev_record.pos;
    ref.seq.insert(ref.seq.begin(), prev_record.ref.seq.begin(), prev_record.ref.seq.end());
    copy_events_and_anti_events(ref, prev_record.ref);

    move_alts(*this, std::move(new_record));

    // remove extended alts which are in parity to the added event
    alts.erase(std::remove_if(alts.begin(),
                              alts.end(),
                              [](Alt const & alt)
      {
        for (auto anti_event : alt.anti_events)
        {
          if (alt.events.count(anti_event) == 1)
            return true;
        }

        return false;
      }), alts.end());
  }
  else
  {
    // the records overlap each other
    merge(std::move(prev_record), 0);
  }
}


void
VarRecord::merge(VarRecord && prev_record, long const EXTRA_SUFFIX)
{
  assert(pos >= prev_record.pos);

  long const jump_size = pos - prev_record.pos; // Original reference allele pos
  long const oref_size = ref.seq.size(); // Original reference allele size

  // Insert new sequences prior to the reference sequence and all alternative sequences (if any)
  if (jump_size > 0)
    insert_prior_sequence(*this, prev_record);

  long const oref_size_pre = ref.seq.size(); // reference allele after adding prefix
  assert(oref_size + jump_size == oref_size_pre);

  // Check if the previous reference reacher further than then current one
  extend_smaller_record(*this, prev_record);

  long const extension_size = ref.seq.size() - oref_size_pre;
  assert(prev_record.ref.seq == ref.seq);

  // Create new variants on the previous record if some this variant does not overlap with any of the previous ones
  VarRecord new_record(prev_record.pos); // its ok that ref is just set to empty string - it is not used

  for (long i{0}; i < static_cast<long>(prev_record.alts.size()); ++i)
  {
    auto & prev_alt = prev_record.alts[i];

    if (static_cast<long>(prev_alt.seq.size()) <= oref_size)
      continue;

    long suffix_matches{0};
    long const offset = ref.seq.size() - prev_alt.seq.size();

    // TODO figure out why this is needed
    if (jump_size - offset < 0)
      continue;

    {
      long const smaller_allele_size = std::min(ref.seq.size(), prev_alt.seq.size());

      for (long k{0}; k < smaller_allele_size; ++k)
      {
        if (*(ref.seq.rbegin() + k) == *(prev_alt.seq.rbegin() + k))
          ++suffix_matches;
        else
          break;
      }
    }

    if (suffix_matches >= extension_size + EXTRA_SUFFIX)
    {
      assert(jump_size - offset >= 0);
      assert(jump_size - offset < static_cast<long>(prev_alt.seq.size()));

      Alt prefix_alt(prev_alt);
      prefix_alt.seq = std::vector<char>(prev_alt.seq.begin(), prev_alt.seq.begin() + jump_size - offset);

      for (long j{0}; j < static_cast<long>(alts.size()); ++j)
      {
        Alt const & curr_alt = alts[j];

        if (is_ok_to_merge_alts(prefix_alt, curr_alt))
        {
          Alt new_alt = make_alt(prefix_alt, curr_alt, jump_size);
          new_record.alts.push_back(std::move(new_alt));
        }
      }
    }
  }

  copy_events_and_anti_events(ref, prev_record.ref);

  // add prev ref events to current alts
  for (Alt & alt : alts)
    copy_events_and_anti_events(alt, prev_record.ref);

  // remove extended alts which are in parity to the added event
  prev_record.alts.erase(
    std::remove_if(prev_record.alts.begin(),
                   prev_record.alts.end(),
                   [this](Alt const & prev_alt)
    {
      for (auto anti_event : prev_alt.anti_events)
      {
        if (ref.events.count(anti_event))
          return true;
      }

      return false;
    }), prev_record.alts.end());

  move_alts(*this, std::move(prev_record));
  move_alts(*this, std::move(new_record));
}


void
VarRecord::add_suffix(std::vector<char> && suffix)
{
  for (Alt & alt : alts)
    std::copy(suffix.begin(), suffix.end(), std::back_inserter(alt.seq));

  std::move(suffix.begin(), suffix.end(), std::back_inserter(ref.seq));
}


std::vector<char>
VarRecord::get_common_suffix() const
{
  if (ref.seq.size() == 0 ||
      std::any_of(alts.begin(),
                  alts.end(),
                  is_empty_seq))
  {
    std::vector<char> empty;
    return empty;
  }

  long suffix_size{0};
  auto ref_it = ref.seq.rbegin();

  while (ref_it != ref.seq.rend() &&
         suffix_size < (static_cast<long>(ref.seq.size()) - 1l) &&
         std::all_of(alts.begin(),
                     alts.end(),
                     [&](Alt const & alt)
    {
      return suffix_size < (static_cast<long>(alt.seq.size()) - 1l) && *(alt.seq.rbegin() + suffix_size) == *ref_it;
    }
                     )
         )
  {
    ++ref_it;
    ++suffix_size;
  }

  assert(suffix_size < static_cast<long>(ref.seq.size()));
  return std::vector<char>(ref.seq.end() - suffix_size, ref.seq.end());
}


bool
VarRecord::is_any_seq_larger_than(long val) const
{
  return static_cast<long>(ref.seq.size()) > val ||
         std::any_of(alts.begin(),
                     alts.end(),
                     [val](Alt const & alt)
    {
      return static_cast<long>(alt.seq.size()) > val;
    });
}


bool
VarRecord::is_snp_or_snps() const
{
  return std::all_of(alts.begin(),
                     alts.end(),
                     [this](Alt const & alt){
      return alt.seq.size() == ref.seq.size();
    });
}


bool
VarRecord::operator<(VarRecord const & b)
{
  if (pos < b.pos)
    return true;
  else
    return false;

  if (alts.size() == 1 && b.alts.size() == 1)
  {
    std::vector<char> const & a_ref = ref.seq;
    std::vector<char> const & b_ref = b.ref.seq;
    std::vector<char> const & a_alt = alts[0].seq;
    std::vector<char> const & b_alt = b.alts[0].seq;

    // order is snp, insertion, deletion
    bool const a_is_del = a_ref.size() > a_alt.size();
    bool const a_is_ins = a_alt.size() > a_ref.size();

    bool const b_is_del = b_ref.size() > b_alt.size();
    bool const b_is_ins = b_alt.size() > b_ref.size();

    auto const order_a = a_is_ins + 2 * a_is_del;
    auto const order_b = b_is_ins + 2 * b_is_del;

    return order_a < order_b;
  }
}


} // namespace gyper
