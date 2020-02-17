#pragma once

#include <cstdint>
#include <vector>

#include <graphtyper/graph/graph.hpp>


namespace gyper
{

class GenotypePaths;
class VariantCandidate;


class ReferenceDepth
{
public:
  ReferenceDepth();

  uint32_t reference_offset = 0;
  std::vector<std::vector<uint16_t> > depths{};

  /***************
   * INFORMATION *
   ***************/
  long start_pos_to_index(long start_pos) const;
  long end_pos_to_index(long end_pos, long depth_size) const;

  /**
   * \brief Retrieve the graph aligned read depth at a VariantCandidate site
   * \param var Variant candidate to check read depth for.
   * \param sample_index the index of the PN to check for.
   */
  uint16_t get_read_depth(VariantCandidate const & var, long sample_index) const;
  uint16_t get_read_depth(uint32_t abs_pos, long sample_index) const;
  uint64_t get_total_read_depth_of_samples(VariantCandidate const & var,
                                           std::vector<uint32_t> const & sample_indexes
                                           ) const;

  /****************
   * MODIFICATION *
   ****************/
  void set_depth_sizes(long sample_count, long reference_size = graph.reference.size());
  void add_depth(long start_pos, long end_pos, long sample_index);
  void add_genotype_paths(GenotypePaths const & geno, long sample_index);
};

} // namespace gyper
