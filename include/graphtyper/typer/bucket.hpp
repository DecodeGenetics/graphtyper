#pragma once

#include <parallel_hashmap/phmap_fwd_decl.h>

#include <graphtyper/typer/event.hpp>
#include <graphtyper/typer/read.hpp>


namespace gyper
{

class BucketFirstPass
{
public:
  int32_t global_max_pos_end{-1};// Max pos end of alignments in this bucket and all previous buckets
  int32_t max_pos_end{-1}; // Max pos end of alignments in this bucket

  std::map<SnpEvent, SnpEventInfo> snps;
  Tindel_events indel_events; // type is std::map<IndelEvent, EventInfo>
};


class Bucket
{
public:
  int32_t global_max_pos_end{-1};// Max pos end of alignments in this bucket and all previous buckets
  int32_t max_pos_end{-1}; // Max pos end of alignments in this bucket

  Tindel_events indel_events; // type is std::map<IndelEvent, EventInfo>
  std::vector<Read> reads;

  std::string to_string() const;
};


bool
is_indel_in_bucket(std::vector<Bucket> const & buckets,
                   IndelEvent const & indel_event,
                   long const region_begin,
                   long const BUCKET_SIZE);


std::map<SnpEvent, SnpEventInfo>::iterator
add_snp_event_to_bucket(std::vector<BucketFirstPass> & buckets,
                        SnpEvent && event,
                        long const region_begin,
                        long const BUCKET_SIZE);


template <typename TBucket>
Tindel_events::iterator
add_indel_event_to_bucket(std::vector<TBucket> & buckets,
                          IndelEvent && event,
                          long const region_begin,
                          long const BUCKET_SIZE,
                          std::vector<char> const & reference_sequence,
                          long ref_offset);

//void
//add_base_to_bucket(std::vector<Bucket> & buckets,
//                   int32_t pos,
//                   char seq,
//                   char qual,
//                   long const region_begin,
//                   long const BUCKET_SIZE);

} // namespace gyper
