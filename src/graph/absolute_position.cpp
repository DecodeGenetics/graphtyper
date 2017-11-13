#include <cassert>

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/absolute_position.hpp>

#include <boost/log/trivial.hpp>

namespace gyper
{

AbsolutePosition::AbsolutePosition()
{
  offsets[0] = 0;
  chromosome_to_offset[chromosome_names[0]] = 0;

  for (std::size_t i = 1; i < 26ul; ++i)
  {
    offsets[i] = offsets[i - 1] + chromosome_lengths[i - 1];
    chromosome_to_offset[chromosome_names[i]] = offsets[i];
  }
}


uint32_t
AbsolutePosition::get_absolute_position(std::string const & chromosome, uint32_t const contig_position) const
{
  std::string chr;

  if (chromosome.substr(0, 3) != std::string("chr"))
    chr = "chr" + chromosome;
  else
    chr = chromosome;

  uint32_t abs_pos;

  try
  {
    abs_pos = chromosome_to_offset.at(chr) + contig_position;
  }
  catch (std::out_of_range const &)
  {
    abs_pos = chromosome_to_offset.at("chrUn") + contig_position;
    BOOST_LOG_TRIVIAL(warning) << "[gyper::graph::absolute_position] No chromosome \"" << chr
                               << "\" available, using \"chrUn\" instead";
  }

  return abs_pos;
}


std::pair<std::string, uint32_t>
AbsolutePosition::get_contig_position(uint32_t const absolute_position) const
{
  auto offset_it = std::lower_bound(offsets.begin(), offsets.end(), absolute_position);
  std::size_t const i = std::distance(offsets.begin(), offset_it);
  assert(i > 0);
  assert(i <= chromosome_names.size());
  return std::make_pair<std::string, uint32_t>(std::string(chromosome_names[i - 1]), absolute_position - offsets[i - 1]);
}


AbsolutePosition absolute_pos;

} // namespace gyper
