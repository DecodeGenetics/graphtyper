#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include <cereal/access.hpp>


namespace gyper
{

class VarRecord;
class Graph;

class GenomicRegion
{
  friend class cereal::access;

public:
  std::string chr;
  long begin{0};
  long end{0};

  /****************
   * CONSTRUCTORS *
   ****************/
  explicit GenomicRegion();
  explicit GenomicRegion(std::string_view region);
  explicit GenomicRegion(std::string_view chrom, long begin, long end);
  explicit GenomicRegion(std::string && chrom, long begin, long end);

  void clear();
  void pad(long N_bases); // pad region by N_bases
  void pad_end(long N_bases); // pad end of region by N_bases
  long get_absolute_begin_position() const;
  long get_absolute_end_position() const;
  long get_absolute_position(std::string const & chromosome, long contig_position) const;
  long get_absolute_position(long contig_position) const;
  std::pair<std::string, long> get_contig_position(long absolute_position, Graph const & graph) const;
  std::string to_string() const;
  std::string to_file_string() const;

  void check_if_var_records_match_reference_genome(std::vector<VarRecord> const & var_records,
                                                   std::vector<char> const & reference);
  void add_reference_to_record_if_they_have_a_matching_prefix(VarRecord & var_record,
                                                              std::vector<char> const & reference);

private:
  template <typename Archive>
  void serialize(Archive & ar, unsigned int);
};

} // namespace gyper
