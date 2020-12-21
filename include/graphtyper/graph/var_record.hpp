#pragma once

#include <cstdint>
#include <vector>

#include <graphtyper/graph/alt.hpp>
#include <graphtyper/graph/ref.hpp>
#include <graphtyper/graph/snp_event.hpp>


std::ostream &
operator<<(std::ostream & os, std::vector<char> const & dt);

namespace gyper
{

class VarRecord
{
public:
  uint32_t pos;
  Ref ref;
  std::vector<Alt> alts;
  bool is_sv{false};

  VarRecord();
  VarRecord(uint32_t const pos);
  VarRecord(uint32_t const pos, Ref && ref, std::vector<Alt> && alts);
  VarRecord(VarRecord const &);
  VarRecord(VarRecord &&);
  ~VarRecord() = default;

  VarRecord & operator=(VarRecord const & o);
  VarRecord & operator=(VarRecord && o);

  std::string to_string() const;

  void clear();
  void merge_one_path(VarRecord && prev_record);
  void merge_all(VarRecord && prev_record);
  void merge(VarRecord && prev_record, long EXTRA_SUFFIX = 0);
  void add_suffix(std::vector<char> && suffix);

  std::vector<char> get_common_suffix() const;
  bool is_any_seq_larger_than(long val) const;
  bool is_snp_or_snps() const;

  bool operator<(VarRecord const & b);
};


} // namespace gyper
