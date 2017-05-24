#pragma once

#include <seqan/vcf_io.h>

namespace gyper
{

seqan::VcfRecord inline
create_empty_record(uint32_t const pos = 0)
{
  seqan::VcfRecord record;
  record.rID = 0;
  record.beginPos = pos;
  record.id = ".";
  record.filter = ".";
  record.info = ".";
  record.format = ".";
  return record;
}


} // namespace gyper
