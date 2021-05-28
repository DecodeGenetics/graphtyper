#include <array>
#include <cassert>
#include <initializer_list>
#include <utility>
#include <vector>

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/alt.hpp>

#include <graphtyper/utilities/logging.hpp>


namespace gyper
{

Alt::Alt(std::vector<char> && _seq)
  : seq(std::forward<std::vector<char> >(_seq))
{}


Alt::Alt(std::vector<char> const & _seq)
  : seq(_seq)
{}


Alt::Alt(std::initializer_list<char> l)
  : seq(l)
{}


void
Alt::clear()
{
  //non_merge_start = -1;
  //non_merge_end = -1;
  seq.clear();
  events.clear();
  //parity_events.clear();
  anti_events.clear();
}


bool
Alt::operator==(Alt const & o)
{
  return seq == o.seq;
}


bool
Alt::operator!=(Alt const & o)
{
  return !(seq == o.seq);
}


bool
Alt::operator<(Alt const & o)
{
  return seq < o.seq;
}


bool
is_empty_seq(Alt const & alt)
{
  return alt.seq.size() == 0;
}


bool
compare_two_alts(Alt const & lhs, Alt const & rhs)
{
  return lhs.seq < rhs.seq;
}


Alt
make_alt(Alt const & prev, Alt const & curr, long const jump_size)
{
  Alt new_alt(prev);
  assert(new_alt.events.size() == prev.events.size());
  assert(new_alt.anti_events.size() == prev.anti_events.size());
  assert(jump_size < static_cast<long>(curr.seq.size()));

  std::copy(curr.seq.cbegin() + jump_size, curr.seq.cend(), std::back_inserter(new_alt.seq));

  std::copy(curr.events.begin(),
            curr.events.end(),
            std::inserter(new_alt.events, new_alt.events.begin()));

  //std::copy(curr.parity_events.begin(),
  //          curr.parity_events.end(),
  //          std::inserter(new_alt.parity_events, new_alt.parity_events.begin()));

  std::copy(curr.anti_events.begin(),
            curr.anti_events.end(),
            std::inserter(new_alt.anti_events, new_alt.anti_events.begin()));

  //if (new_alt.non_merge_start == -1)
  //{
  //  new_alt.non_merge_start = curr.non_merge_start;
  //  new_alt.non_merge_end = curr.non_merge_end;
  //}
  //else if (curr.non_merge_start != -1)
  //{
  //  // neither is -1
  //  new_alt.non_merge_start = std::min(new_alt.non_merge_start, curr.non_merge_start);
  //  new_alt.non_merge_end = std::max(new_alt.non_merge_end, curr.non_merge_end);
  //}

  return new_alt;
}


bool
is_ok_to_merge_alts(Alt const & prev_alt, Alt const & curr_alt)
{
  //return true;

  //if (curr_alt.events.size() > 0)
  //{
  //  BOOST_LOG_TRIVIAL(info) << __HERE__ << " " << curr_alt.events.size();
  //}

  //if (prev_alt.non_merge_start != -1 &&
  //    curr_alt.non_merge_start != -1 &&
  //    prev_alt.non_merge_end >= curr_alt.non_merge_start)
  //{
  //  //BOOST_LOG_TRIVIAL(info) << __HERE__ << " Cannot merge " << std::string(curr_alt.seq.begin(), curr_alt.seq.end())
  //  //                        << " and "
  //  //                        << std::string(prev_alt.seq.begin(), prev_alt.seq.end())
  //  //                        << " due to too close indels";
  //  return false;
  //}

  // Do not merge if the current alt is an anti event
  for (auto const & curr_event : curr_alt.events)
  {
    //BOOST_LOG_TRIVIAL(info) << __HERE__ << " " << curr_event.to_string();
    if (curr_event < 0)
      continue;

    if (prev_alt.anti_events.count(curr_event) == 1)
    {
      //BOOST_LOG_TRIVIAL(info) << __HERE__ << " Cannot merge "
      //                        << std::string(curr_alt.seq.begin(), curr_alt.seq.end())
      //                        << " and "
      //                        << std::string(prev_alt.seq.begin(), prev_alt.seq.end())
      //                        << " due to anti event "
      //                        << curr_event;

      return false;
    }


  }

#ifndef NDEBUG
  // make sure it is not the other way around either
  for (auto const & prev_event : prev_alt.events)
  {
    assert(curr_alt.anti_events.count(prev_event) == 0);
  }
#endif

  return true;
}


} // namespace gyper
