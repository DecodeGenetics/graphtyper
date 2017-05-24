#include <cassert>

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/absolute_position.hpp>


namespace gyper
{

AbsolutePosition::AbsolutePosition()
{
  offsets[0] = 0;
  chromosome_to_offset[chromosome_names[0]] = 0;

  for (std::size_t i = 1; i < 24ul; ++i)
  {
    offsets[i] = offsets[i - 1] + chromosome_lengths[i - 1];
    chromosome_to_offset[chromosome_names[i]] = offsets[i];
  }
}


uint32_t
AbsolutePosition::get_absolute_position(std::string const & chromosome, uint32_t const contig_position) const
{
  return chromosome_to_offset.at(chromosome) + contig_position;
}


std::pair<std::string, uint32_t>
AbsolutePosition::get_contig_position(uint32_t const absolute_position) const
{
  // std::cout << "absolute_position = " << absolute_position << std::endl;
  std::size_t i = std::distance(offsets.begin(),
                                std::lower_bound(offsets.begin(),
                                                 offsets.end(),
                                                 absolute_position
                                                 )
                                );

  assert(i > 0);
  assert(static_cast<std::size_t>(i - 1ul) < chromosome_names.size());
  return {
           chromosome_names[i - 1], absolute_position - offsets[i - 1]
  };
}


AbsolutePosition absolute_pos;

} // namespace gyper
