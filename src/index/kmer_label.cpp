#include <cstdint>

#include <graphtyper/graph/graph.hpp>
#include <graphtyper/index/kmer_label.hpp>

namespace gyper
{

bool
KmerLabel::operator==(KmerLabel const & c2) const
{
  return start_index == c2.start_index &&
         end_index == c2.end_index &&
         variant_id == c2.variant_id;
}


bool
KmerLabel::operator!=(KmerLabel const & c2) const
{
  return !(*this == c2);
}


} // namespace gyper
