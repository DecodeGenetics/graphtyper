#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include <cereal/access.hpp>


namespace gyper
{

struct Contig;

class AbsolutePosition
{
public:
  std::vector<long> offsets;
  std::unordered_map<std::string, long> chromosome_to_offset;

  AbsolutePosition() = default;
  AbsolutePosition(std::vector<Contig> const & contigs);

  ///* non-const methods */
  void calculate_offsets(std::vector<Contig> const & contigs);

  ///* const methods */
  bool is_contig_available(std::string const & chromosome) const;
  long get_absolute_position(std::string const & chromosome, long contig_position) const;

  std::pair<std::string, long>
  get_contig_position(long absolute_position,
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
