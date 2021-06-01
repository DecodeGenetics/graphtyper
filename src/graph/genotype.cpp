#include <utility>

#include <cereal/archives/binary.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/vector.hpp>

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/genotype.hpp>

namespace gyper
{
Genotype::Genotype() : id(0), num(1), first_variant_node(0)
{
}

Genotype::Genotype(uint32_t const i, uint16_t const n, uint32_t const fvn) : id(i), num(n), first_variant_node(fvn)
{
}

template <typename Archive>
void Genotype::serialize(Archive & ar, const unsigned int /*version*/)
{
  ar & id;
  ar & num;
  ar & first_variant_node;
}

/***************************
 * EXPLICIT INSTANTIATIONS *
 ***************************/

template void Genotype::serialize<cereal::BinaryInputArchive>(cereal::BinaryInputArchive &, const unsigned int);
template void Genotype::serialize<cereal::BinaryOutputArchive>(cereal::BinaryOutputArchive &, const unsigned int);

} // namespace gyper
