#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <boost/log/trivial.hpp>

#include <parallel_hashmap/phmap.h>

#include <graphtyper/constants.hpp>
#include <graphtyper/typer/bucket.hpp>
#include <graphtyper/typer/event.hpp>
#include <graphtyper/typer/read.hpp>


namespace gyper
{

std::string
Bucket::to_string() const
{
  std::ostringstream ss;
  ss << " indel events " << "," << events.size();
  ss << " n_reads=" << reads.size();

  for (auto const & read : reads)
    ss << "\n" << read.to_string();

  return ss.str();
}


bool
is_indel_in_bucket(std::vector<Bucket> const & buckets,
                   Event const & indel_event,
                   long const region_begin,
                   long const BUCKET_SIZE)
{
  long const event_bucket_index = (indel_event.pos - region_begin) / BUCKET_SIZE;
  assert(event_bucket_index < static_cast<long>(buckets.size()));
  auto & indel_events = buckets[event_bucket_index].events;
  auto find_it = indel_events.find(indel_event);
  return find_it != indel_events.end();
}


template <typename TBucket>
std::map<Event, EventSupport>::iterator
add_indel_event_to_bucket(std::vector<TBucket> & buckets,
                          Event && new_indel_event,
                          long const region_begin,
                          long const BUCKET_SIZE,
                          std::vector<char> const & reference_sequence,
                          long ref_offset)
{
  long const REF_SIZE = reference_sequence.size();

  assert(new_indel_event.pos >= region_begin);
  long const event_bucket_index = (new_indel_event.pos - region_begin) / BUCKET_SIZE;
  assert(event_bucket_index >= 0);

  if (event_bucket_index >= static_cast<long>(buckets.size()))
    buckets.resize(event_bucket_index + 1);

  assert(event_bucket_index < static_cast<long>(buckets.size()));

  EventSupport new_info;
  std::pair<std::map<Event, EventSupport>::iterator, bool> it_pair =
    buckets[event_bucket_index].events.insert({std::move(new_indel_event), std::move(new_info)});

  assert(it_pair.first != buckets[event_bucket_index].events.end());

  if (it_pair.second)
  {
    Event const & new_event = it_pair.first->first;
    EventSupport & new_event_info = it_pair.first->second;

    long span{0};
    long const count{static_cast<long>(new_event.sequence.size())};

    // Calculate span
    if (new_event.type == 'I' || new_event.type == 'i')
    {
      // Check insterted sequence
      while (span < count)
      {
        if ((ref_offset + span) >= REF_SIZE ||
            new_event.sequence[span] != reference_sequence[ref_offset + span])
        {
          break;
        }

        ++span;
      }

      // Check further if we have matched the entire inserted sequence
      if (span == count)
      {
        while ((ref_offset + span) < REF_SIZE)
        {
          if (reference_sequence[ref_offset + span - count] != reference_sequence[ref_offset + span])
            break;

          ++span;
        }
      }

      // Check for overflow
      if ((span + 1) >= std::numeric_limits<uint16_t>::max())
        span = std::numeric_limits<uint16_t>::max() - 1;

      new_event_info.span = span + 1; // to 1-based system
    }
    else
    {
      assert(new_event.type == 'D' || new_event.type == 'd');

      // Check after deleted sequence
      while ((ref_offset + span) < REF_SIZE)
      {
        if (reference_sequence[ref_offset + span] != reference_sequence[ref_offset + span + count])
          break;

        ++span;
      }

      // Check for overflow
      if ((span + 1) >= std::numeric_limits<uint16_t>::max())
        span = std::numeric_limits<uint16_t>::max() - 1;

      new_event_info.span = span + 1; // to 1-based system
    }
  }

  return it_pair.first;
}


std::map<Event, EventSupport>::iterator
add_snp_event_to_bucket(std::vector<BucketFirstPass> & buckets,
                        Event && event,
                        long const region_begin,
                        long const BUCKET_SIZE)
{
  assert(event.pos >= region_begin);
  long const event_bucket_index = (event.pos - region_begin) / BUCKET_SIZE;
  assert(event_bucket_index >= 0);

  if (event_bucket_index >= static_cast<long>(buckets.size()))
    buckets.resize(event_bucket_index + 1);

  assert(event_bucket_index < static_cast<long>(buckets.size()));

  EventSupport new_info;
  std::pair<std::map<Event, EventSupport>::iterator, bool> it_pair =
    buckets[event_bucket_index].events.insert({std::move(event), std::move(new_info)});

  return it_pair.first;
}


/*
void
add_base_to_bucket(std::vector<Bucket> & buckets,
                   int32_t pos,
                   char seq,
                   char qual,
                   long const region_begin,
                   long const BUCKET_SIZE)
{
  // Any sensible compiler will optimize this
  long const event_bucket_index = (pos - region_begin) / BUCKET_SIZE;
  long const local_pos = (pos - region_begin) % BUCKET_SIZE;

  if (event_bucket_index >= static_cast<long>(buckets.size()))
    buckets.resize(event_bucket_index + 1);

  auto & bucket = buckets[event_bucket_index];

  if (bucket.pileup.size() == 0)
    bucket.pileup.resize(BUCKET_SIZE);

  assert(local_pos < static_cast<long>(bucket.pileup.size()));
  bucket.pileup[local_pos].add_base(seq, qual);
}
*/


// explicit instantiation
template std::map<Event, EventSupport>::iterator
add_indel_event_to_bucket(std::vector<BucketFirstPass> & buckets,
                          Event && new_indel_event,
                          long const region_begin,
                          long const BUCKET_SIZE,
                          std::vector<char> const & reference_sequence,
                          long ref_offset);

template std::map<Event, EventSupport>::iterator
add_indel_event_to_bucket(std::vector<Bucket> & buckets,
                          Event && new_indel_event,
                          long const region_begin,
                          long const BUCKET_SIZE,
                          std::vector<char> const & reference_sequence,
                          long ref_offset);


} // namespace gyper
