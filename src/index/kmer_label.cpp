#include <cstdint>
#include <sstream>
#include <string>

#include <cereal/archives/binary.hpp>

#include <graphtyper/graph/graph.hpp>
#include <graphtyper/index/kmer_label.hpp>

namespace gyper
{

KmerLabel::KmerLabel(uint32_t const s, uint32_t const e) noexcept
  : start_index(s)
  , end_index(e)
  , variant_id(INVALID_ID)
{}

KmerLabel::KmerLabel(uint32_t const s, uint32_t const e, uint32_t const i) noexcept
  : start_index(s)
  , end_index(e)
  , variant_id(i)
{}

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


template <class Archive>
void
KmerLabel::serialize(Archive & ar, const unsigned int)
{
  ar & start_index;
  ar & end_index;
  ar & variant_id;
}


std::string
KmerLabel::to_string() const
{
  std::ostringstream ss;
  ss << "start,end,ID=" << start_index << "," << end_index << "," << variant_id;
  return ss.str();
}


/***************************
 * EXPLICIT INSTANTIATIONS *
 ***************************/
template void KmerLabel::serialize<cereal::BinaryInputArchive>(cereal::BinaryInputArchive &,
                                                                    const unsigned int);
template void KmerLabel::serialize<cereal::BinaryOutputArchive>(cereal::BinaryOutputArchive &,
                                                                    const unsigned int);

} // namespace gyper
