#pragma once

#include <cstdint>
#include <string>
#include <unordered_set>
#include <vector>

#include <graphtyper/graph/genomic_region.hpp>

namespace gyper
{
class GenotypePaths;
class Graph;

class Primers
{
  std::vector<GenomicRegion> left;
  std::vector<GenomicRegion> right;

  std::unordered_set<uint32_t> left_var_orders;
  std::unordered_set<uint32_t> right_var_orders;

public:
  Primers() = delete;
  Primers(std::string const & primer_bedpe);

  //* non-const */
  // Read the primer BEDPE
  void read(std::string const & primer_bedpe);
  void reset_var_orders_to_check(Graph const & graph);

  //* const */
  void check(GenotypePaths & genos) const;
  void check_left(GenotypePaths & genos) const;
  void check_right(GenotypePaths & genos) const;
};

} // namespace gyper
