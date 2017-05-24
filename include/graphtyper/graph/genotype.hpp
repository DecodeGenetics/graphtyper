#pragma once

#include <cstdint>

#include <boost/serialization/access.hpp>


namespace gyper
{

class Genotype
{
  friend class boost::serialization::access;

public:
  uint32_t id;
  uint16_t num;
  uint32_t first_variant_node;

  Genotype();
  Genotype(uint32_t const i, uint16_t const n, uint32_t const fvn);

private:
  template <class Archive>
  void serialize(Archive & ar, const unsigned int version);
};

} // namespace gyper
