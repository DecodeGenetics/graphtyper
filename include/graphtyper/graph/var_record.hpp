#pragma once

#include <iostream>
#include <vector>

#include <htslib/kseq.h>

std::ostream & operator<<(std::ostream & os, std::vector<char> const & dt);
std::string to_string(std::vector<char> const & dt);

namespace gyper
{

class VarRecord
{
public:
  uint32_t pos;
  std::vector<char> ref;
  std::vector<std::vector<char> > alts;
  bool is_sv = false;

  VarRecord();
  VarRecord(uint32_t const pos);
  VarRecord(uint32_t const pos, std::vector<char> && ref, std::vector<std::vector<char> > && alts);

  std::string to_string() const;

  void clear();
  void merge_one_path(VarRecord && prev_record);
  void merge(VarRecord && prev_record, long EXTRA_SUFFIX = 0);

  std::vector<char> get_common_suffix();
};

} // namespace gyper
