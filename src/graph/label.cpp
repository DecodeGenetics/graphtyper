#include <assert.h>

#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>

#include <graphtyper/graph/label.hpp>
#include <graphtyper/graph/var_record.hpp>

namespace gyper
{
Label::Label() noexcept : order(INVALID_ID), dna(0), variant_num(INVALID_NUM)
{
}

Label::Label(Label const & l) noexcept : order(l.order), dna(l.dna), variant_num(l.variant_num)
{
}

Label::Label(Label && l) noexcept :
  order(std::forward<uint32_t>(l.order)),
  dna(std::forward<std::vector<char>>(l.dna)),
  variant_num(std::forward<uint16_t>(l.variant_num))
{
}

Label::Label(uint32_t const _order, std::vector<char> && _dna, uint32_t const _variant_num) noexcept :
  order(_order), dna(std::move(_dna)), variant_num(_variant_num)
{
}

Label::Label(uint32_t const _order, Alt && alt, uint32_t const _variant_num) :
  order(_order), dna(std::move(alt.seq)), variant_num(_variant_num)
{
  // TODO do stuff with events and anti_events
}

uint32_t Label::reach() const
{
  return order + dna.size() - 1;
}

template <typename Archive>
void Label::serialize(Archive & ar, const unsigned int)
{
  ar & order;
  ar & dna;
  ar & variant_num;
}

/***************************
 * EXPLICIT INSTANTIATIONS *
 ***************************/

template void Label::serialize<cereal::BinaryInputArchive>(cereal::BinaryInputArchive &, const unsigned int);
template void Label::serialize<cereal::BinaryOutputArchive>(cereal::BinaryOutputArchive &, const unsigned int);

} // namespace gyper
