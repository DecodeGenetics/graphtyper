#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <parallel_hashmap/phmap.h>

#include <graphtyper/constants.hpp>
#include <graphtyper/typer/bucket.hpp>
#include <graphtyper/typer/event.hpp>
#include <graphtyper/typer/read.hpp>
#include <graphtyper/utilities/logging.hpp>
#include <graphtyper/utilities/options.hpp>

namespace gyper
{
std::string Bucket::to_string() const
{
  std::ostringstream ss;
  ss << " indel events "
     << "," << events.size();
  ss << " n_reads=" << reads.size();

  for (auto const & read : reads)
    ss << "\n" << read.to_string();

  return ss.str();
}

void merge_bucket_lr(std::vector<BucketLR> & old_buckets, std::vector<BucketLR> const & new_buckets)
{
  long const max_bucket_size = std::max(old_buckets.size(), new_buckets.size());

  for (long b{0}; b < max_bucket_size; ++b)
  {
    if (b >= static_cast<long>(new_buckets.size()) || new_buckets[b].pileup.size() == 0)
      return;
    else if (b >= static_cast<long>(old_buckets[b].pileup.size()) || old_buckets[b].pileup.size() == 0)
      old_buckets[b].pileup = new_buckets[b].pileup;

    auto & old_bucket = old_buckets[b];
    auto const & new_bucket = new_buckets[b];
    assert(old_bucket.pileup.size() == new_bucket.pileup.size());

    for (long p{0}; p < static_cast<long>(old_bucket.pileup.size()); ++p)
    {
      auto & old_bc = old_bucket.pileup[p];
      auto const & new_bc = new_bucket.pileup[p];

      old_bc.acgt[0] += new_bc.acgt[0];
      old_bc.acgt[2] += new_bc.acgt[1];
      old_bc.acgt[3] += new_bc.acgt[2];
      old_bc.acgt[3] += new_bc.acgt[3];

      old_bc.acgt_qualsum[0] += new_bc.acgt_qualsum[0];
      old_bc.acgt_qualsum[1] += new_bc.acgt_qualsum[1];
      old_bc.acgt_qualsum[2] += new_bc.acgt_qualsum[2];
      old_bc.acgt_qualsum[3] += new_bc.acgt_qualsum[3];

      old_bc.deleted += new_bc.deleted;
      old_bc.unknown += new_bc.unknown;
    }
  }
}

bool is_indel_in_bucket(std::vector<Bucket> const & buckets,
                        Event const & indel_event,
                        long const region_begin,
                        long const BUCKET_SIZE)
{
  long const event_bucket_index = (indel_event.pos - region_begin) / BUCKET_SIZE;

  if (event_bucket_index >= static_cast<long>(buckets.size()))
    return false;

  auto & indel_events = buckets[event_bucket_index].events;
  auto find_it = indel_events.find(indel_event);
  return find_it != indel_events.end();
}

template <typename TBucket>
std::map<Event, EventSupport>::iterator add_indel_event_to_bucket(std::vector<TBucket> & buckets,
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
        if ((ref_offset + span) >= REF_SIZE || new_event.sequence[span] != reference_sequence[ref_offset + span])
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

std::map<Event, EventSupport>::iterator add_snp_event_to_bucket(std::vector<BucketFirstPass> & buckets,
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

bool add_base_to_bucket(
  std::vector<BucketLR> & buckets, int32_t pos, char seq, char qual, long const region_begin, long const BUCKET_SIZE)
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
  BaseCount & local_pileup = bucket.pileup[local_pos];
  local_pileup.add_base(seq, qual);
  int const c = Options::const_instance()->lr_coverage_filter;
  return c > 0 && local_pileup.get_depth_with_deleted() >= c;
}

// explicit instantiation
template std::map<Event, EventSupport>::iterator add_indel_event_to_bucket(std::vector<BucketFirstPass> & buckets,
                                                                           Event && new_indel_event,
                                                                           long const region_begin,
                                                                           long const BUCKET_SIZE,
                                                                           std::vector<char> const & reference_sequence,
                                                                           long ref_offset);

template std::map<Event, EventSupport>::iterator add_indel_event_to_bucket(std::vector<Bucket> & buckets,
                                                                           Event && new_indel_event,
                                                                           long const region_begin,
                                                                           long const BUCKET_SIZE,
                                                                           std::vector<char> const & reference_sequence,
                                                                           long ref_offset);

} // namespace gyper
