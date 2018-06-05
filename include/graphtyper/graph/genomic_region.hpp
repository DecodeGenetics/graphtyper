#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <boost/serialization/access.hpp>

#include <graphtyper/constants.hpp>


namespace gyper
{

static const uint32_t AS_LONG_AS_POSSIBLE = 0xFFFFFFFFULL;

class VarRecord;

class GenomicRegion
{
  friend class boost::serialization::access;

public:
  uint16_t rID;
  std::string chr;
  uint32_t begin;
  uint32_t end;
  uint32_t region_to_refnode;

  // Not serialized
  std::array<uint32_t, 26u> offsets;

  /****************
   * CONSTRUCTORS *
   ****************/
  explicit GenomicRegion(uint32_t _region_to_refnode = 0);
  GenomicRegion(uint16_t && r, std::string && c, uint32_t && b, uint32_t && e, uint32_t _region_to_refnode = 0);
  explicit GenomicRegion(std::string region, uint32_t _region_to_refnode = 0);

  uint32_t get_absolute_begin_position() const;
  uint32_t get_absolute_end_position() const;

  uint32_t get_absolute_position(std::string const & chromosome, uint32_t contig_position) const;
  uint32_t get_absolute_position(uint32_t contig_position) const;
  std::pair<std::string, uint32_t> get_contig_position(uint32_t absolute_position) const;

  std::string to_string() const;

  void check_if_var_records_match_reference_genome(std::vector<VarRecord> const & var_records, std::vector<char> const & reference);
  void add_reference_to_record_if_they_have_a_matching_prefix(VarRecord & var_record, std::vector<char> const & reference);

private:
  template <typename Archive>
  void serialize(Archive & ar, unsigned int);
};

} // namespace gyper
