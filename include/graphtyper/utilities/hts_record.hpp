#pragma once

#include <graphtyper/utilities/hts_store.hpp>
#include <graphtyper/utilities/hts_utils.hpp>

#include <htslib/sam.h>

namespace gyper
{
struct HtsRecord
{
  bam1_t * record = nullptr;
  int file_index = -1;

  HtsRecord() : record(nullptr), file_index(-1)
  {
  }

  HtsRecord(bam1_t * r, int i) : record(r), file_index(i)
  {
  }

  HtsRecord(HtsRecord const &) = delete;
  HtsRecord & operator=(HtsRecord const &) = delete;

  HtsRecord & operator=(HtsRecord && o)
  {
    if (record)
      bam_destroy1(record);

    record = o.record;
    file_index = o.file_index;
    o.record = nullptr;
    return *this;
  }

  HtsRecord(HtsRecord && o)
  {
    if (record)
      bam_destroy1(record);

    record = o.record;
    file_index = o.file_index;
    o.record = nullptr;
  }

  ~HtsRecord()
  {
    if (record)
      bam_destroy1(record);
  }
};

inline bool Cmp_gt_pair_bam1_t_fun(HtsRecord const & a, HtsRecord const & b)
{
  return gt_pos_seq(a.record, b.record);
}

} // namespace gyper
