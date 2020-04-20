#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <boost/serialization/access.hpp>


namespace gyper
{

class VarRecord;
class Graph;

class GenomicRegion
{
  friend class boost::serialization::access;

public:
  std::string chr;
  uint32_t begin{0};
  uint32_t end{0};

  /****************
   * CONSTRUCTORS *
   ****************/
  explicit GenomicRegion();
  explicit GenomicRegion(std::string const & region);
  explicit GenomicRegion(std::string const & chrom, long begin, long end);

  void clear();
  void pad(long N_bases); // pad region by N_bases
  void pad_end(long N_bases); // pad end of region by N_bases
  uint32_t get_absolute_begin_position() const;
  uint32_t get_absolute_end_position() const;
  uint32_t get_absolute_position(std::string const & chromosome, uint32_t contig_position) const;
  uint32_t get_absolute_position(uint32_t contig_position) const;
  std::pair<std::string, uint32_t> get_contig_position(uint32_t absolute_position, Graph const & graph) const;
  std::string to_string() const;

  void check_if_var_records_match_reference_genome(std::vector<VarRecord> const & var_records,
                                                   std::vector<char> const & reference);
  void add_reference_to_record_if_they_have_a_matching_prefix(VarRecord & var_record,
                                                              std::vector<char> const & reference);

private:
  template <typename Archive>
  void serialize(Archive & ar, unsigned int);
};

} // namespace gyper
