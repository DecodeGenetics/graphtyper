#pragma once

#include <string>
#include <unordered_map>
#include <vector>


namespace gyper
{

class AbsolutePosition
{
public:
  AbsolutePosition();

  void calculate_offsets();
  bool is_contig_available(std::string const & chromosome) const;
  uint32_t get_absolute_position(std::string const & chromosome, uint32_t contig_position) const;
  std::pair<std::string, uint32_t> get_contig_position(uint32_t absolute_position) const;

  std::vector<uint32_t> offsets;
  std::unordered_map<std::string, uint32_t> chromosome_to_offset;

};

extern AbsolutePosition absolute_pos;

} // namespace gyper
