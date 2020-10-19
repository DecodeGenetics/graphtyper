#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include <boost/serialization/access.hpp>


namespace gyper
{

struct Contig;

class AbsolutePosition
{
public:
  std::vector<uint32_t> offsets;
  std::unordered_map<std::string, uint32_t> chromosome_to_offset;

  AbsolutePosition() = default;
  AbsolutePosition(std::vector<Contig> const & contigs);

  ///* non-const methods */
  void calculate_offsets(std::vector<Contig> const & contigs);

  ///* const methods */
  bool is_contig_available(std::string const & chromosome) const;
  uint32_t get_absolute_position(std::string const & chromosome, uint32_t contig_position) const;

  std::pair<std::string, uint32_t>
  get_contig_position(uint32_t absolute_position,
                      std::vector<Contig> const & contigs) const;


  template <typename Archive>
  void
  serialize(Archive & ar, unsigned int)
  {
    ar & offsets;
    ar & chromosome_to_offset;
  }


};

extern AbsolutePosition absolute_pos;

} // namespace gyper
