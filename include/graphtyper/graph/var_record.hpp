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

  VarRecord();
  // VarRecord(uint32_t pos, std::vector<char> ref, std::vector<std::vector<char> > alts);
  VarRecord(uint32_t && pos, std::vector<char> && ref, std::vector<std::vector<char> > && alts);
  void clear();

  void merge_one_path(VarRecord && prev_record);
  void merge(VarRecord && prev_record);

  std::vector<char> get_common_suffix();
};

} // namespace gyper

void print_var_record(gyper::VarRecord const & var_record);
