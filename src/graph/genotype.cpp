#include <utility>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/genotype.hpp>


namespace gyper
{

Genotype::Genotype()
  : id(0)
  , num(1)
  , first_variant_node(0)
{}


Genotype::Genotype(uint32_t const i, uint16_t const n, uint32_t const fvn)
  : id(i)
  , num(n)
  , first_variant_node(fvn)
{}


template <typename Archive>
void
Genotype::serialize(Archive & ar, const unsigned int /*version*/)
{
  ar & id;
  ar & num;
  ar & first_variant_node;
}


/***************************
 * EXPLICIT INSTANTIATIONS *
 ***************************/

template void Genotype::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive &, const unsigned int);
template void Genotype::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive &, const unsigned int);

} // namespace gyper
