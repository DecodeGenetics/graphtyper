#include <cassert>
#include <cstdlib>
#include <stdlib.h>
#include <string>

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/absolute_position.hpp>
#include <graphtyper/graph/graph.hpp>
#include <graphtyper/utilities/logging.hpp>

namespace gyper
{
AbsolutePosition::AbsolutePosition(std::vector<Contig> const & contigs)
{
  this->calculate_offsets(contigs);
}

void AbsolutePosition::calculate_offsets(std::vector<Contig> const & contigs)
{
  if (contigs.size() == 0 || contigs.size() == offsets.size())
    return;

  offsets.clear();
  chromosome_to_offset.clear();

  offsets.resize(contigs.size());
  offsets[0] = 0;
  chromosome_to_offset[contigs[0].name] = 0;

  for (long i = 1; i < static_cast<long>(offsets.size()); ++i)
  {
    offsets[i] = offsets[i - 1] + contigs[i - 1].length;
    chromosome_to_offset[contigs[i].name] = offsets[i];
  }
}

bool AbsolutePosition::is_contig_available(std::string const & contig) const
{
  return chromosome_to_offset.count(contig) > 0;
}

long AbsolutePosition::get_absolute_position(std::string const & chromosome, long const contig_position) const
{
  long abs_pos{0};

  try
  {
    abs_pos = chromosome_to_offset.at(chromosome) + contig_position;
  }
  catch (std::out_of_range const &)
  {
    std::cerr << "[gyper::absolute_position] ERROR: No chromosome \"" << chromosome
              << "\" available. Available chromosomes are:\n";

    for (auto it = chromosome_to_offset.begin(); it != chromosome_to_offset.end(); ++it)
      std::cerr << it->first << "\n";

    std::cerr << std::endl;
    assert(false);
    std::exit(113);
  }

  return abs_pos;
}

std::pair<std::string, long> AbsolutePosition::get_contig_position(long const absolute_position,
                                                                   std::vector<Contig> const & contigs) const
{
  auto offset_it = std::lower_bound(offsets.begin(), offsets.end(), absolute_position);
  long const i = std::distance(offsets.begin(), offset_it);
  assert(i > 0);
  assert(i <= static_cast<long>(contigs.size()));
  return std::make_pair<std::string, long>(std::string(contigs[i - 1].name), absolute_position - offsets[i - 1]);
}

AbsolutePosition absolute_pos;

} // namespace gyper
