#pragma once

#include <string>
#include <unordered_map>
#include <vector>


namespace gyper
{

class Graph;

class AbsolutePosition
{
public:
  std::vector<uint32_t> offsets;
  std::unordered_map<std::string, uint32_t> chromosome_to_offset;

  AbsolutePosition() = default;
  AbsolutePosition(Graph const & graph);

  ///* non-const methods */
  void calculate_offsets(Graph const & graph);

  ///* const methods */
  bool is_contig_available(std::string const & chromosome) const;
  uint32_t get_absolute_position(std::string const & chromosome, uint32_t contig_position) const;

  std::pair<std::string, uint32_t>
  get_contig_position(uint32_t absolute_position,
                      Graph const & graph) const;


};

extern AbsolutePosition absolute_pos;

} // namespace gyper
