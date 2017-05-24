#include <assert.h>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include <graphtyper/graph/label.hpp>

namespace gyper
{

Label::Label() noexcept
  : order(INVALID_ID)
  , dna(0)
  , variant_num(INVALID_NUM)
{}


Label::Label(Label const & l) noexcept
  : order(l.order)
  , dna(l.dna)
  , variant_num(l.variant_num)
{}


Label::Label(Label && l) noexcept
  : order(std::forward<uint32_t>(l.order))
  , dna(std::forward<std::vector<char> >(l.dna))
  , variant_num(std::forward<uint16_t>(l.variant_num))
{}


Label::Label(uint32_t const & _order, std::vector<char> && _dna, uint16_t const & _variant_num) noexcept
  : order(_order)
  , dna(std::move(_dna))
  , variant_num(_variant_num)
{}


uint32_t
Label::reach() const
{
  return order + dna.size() - 1;
}


template <typename Archive>
void
Label::serialize(Archive & ar, const unsigned int)
{
  ar & order;
  ar & dna;
  ar & variant_num;
}


/***************************
 * EXPLICIT INSTANTIATIONS *
 ***************************/

template void Label::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive &, const unsigned int);
template void Label::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive &, const unsigned int);

} // namespace gyper
