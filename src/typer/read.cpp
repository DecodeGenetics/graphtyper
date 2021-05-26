#include <cstdint>
#include <string>
#include <vector>

#include <graphtyper/utilities/logging.hpp>

#include <parallel_hashmap/phmap.h>

#include <graphtyper/constants.hpp>
#include <graphtyper/typer/read.hpp>


namespace gyper
{

bool
Alignment::has_indel_event(Tindel_events::iterator indel_event) const
{
  for (auto const & e : indel_events)
  {
    if (e.event_it == indel_event)
      return e.read_pos != READ_ANTI_SUPPORT;
  }

  return false;
}


bool
Alignment::is_clipped() const
{
  return num_clipped_begin != 0 || num_clipped_end != 0;
}


void
Alignment::add_indel_event(int32_t pos, uint16_t flags, uint8_t mapq, Tindel_events::iterator indel_it)
{
  if (pos == READ_ANTI_SUPPORT)
  {
    ++indel_it->second.anti_count;
  }
  else if (pos == READ_MULTI_SUPPORT)
  {
    ++indel_it->second.multi_count;
  }
  else
  {
    auto & indel_info = indel_it->second;
    ++indel_info.hq_count;

    if ((flags & IS_SEQ_REVERSED) != 0u)
      ++indel_info.sequence_reversed;

    if ((flags & IS_PROPER_PAIR) != 0u)
      ++indel_info.proper_pairs;

    if (mapq < 255 && mapq > indel_info.max_mapq)
      indel_info.max_mapq = mapq;
  }

  indel_events.push_back({pos, indel_it});
}


void
Alignment::replace_indel_events(uint16_t flags, uint8_t mapq, std::vector<ReadIndelEvent> && new_indel_events)
{
  // Lower count on old events
  for (auto & event_and_info : indel_events)
  {
    if (event_and_info.read_pos == READ_ANTI_SUPPORT)
    {
      --event_and_info.event_it->second.anti_count;
    }
    else if (event_and_info.read_pos == READ_MULTI_SUPPORT)
    {
      --event_and_info.event_it->second.multi_count;
    }
    else
    {
      --event_and_info.event_it->second.hq_count;

      if ((flags & IS_SEQ_REVERSED) != 0u)
      {
        if (event_and_info.event_it->second.sequence_reversed > 0)
          --event_and_info.event_it->second.sequence_reversed;
      }

      if ((flags & IS_PROPER_PAIR) != 0u)
      {
        //assert(event_and_info.event_it->second.proper_pairs > 0);
        if (event_and_info.event_it->second.proper_pairs > 0)
          --event_and_info.event_it->second.proper_pairs;
      }
    }
  }

  // Increase count on new events
  for (auto & event_and_info : new_indel_events)
  {
    if (event_and_info.read_pos == READ_ANTI_SUPPORT)
    {
      ++event_and_info.event_it->second.anti_count;
    }
    else if (event_and_info.read_pos == READ_MULTI_SUPPORT)
    {
      ++event_and_info.event_it->second.multi_count;
    }
    else
    {
      ++event_and_info.event_it->second.hq_count;

      if ((flags & IS_SEQ_REVERSED) != 0u)
        ++event_and_info.event_it->second.sequence_reversed;

      if ((flags & IS_PROPER_PAIR) != 0u)
        ++event_and_info.event_it->second.proper_pairs;

      if (mapq < 255 && mapq > event_and_info.event_it->second.max_mapq)
        event_and_info.event_it->second.max_mapq = mapq;
    }
  }

  indel_events = std::move(new_indel_events);
}


std::string
Read::to_string() const
{
  std::ostringstream ss;
  ss << name << " " << alignment.pos << " " << flags << " " << static_cast<long>(mapq) << "\n";
  ss << std::string(sequence.begin(), sequence.end());

  //for (auto const & q : qual)
  //  ss << static_cast<char>(q + 33);

  {
    ss << "\nAlignment: score=" << alignment.score << " num_clipped_begin,_end=" << alignment.num_clipped_begin
       << "," << alignment.num_clipped_end << " num_ins_begin=" << alignment.num_ins_begin;

    for (auto const & support_event : alignment.indel_events)
      ss << " pos=" << static_cast<long>(support_event.read_pos) << " " << &(*(support_event.event_it));

    ss << " ";
  }

  return ss.str();
}


} // namespace gyper
