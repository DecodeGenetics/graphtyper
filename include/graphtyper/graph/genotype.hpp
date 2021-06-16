#pragma once

#include <cstdint>

#include <cereal/access.hpp>

namespace gyper
{
class Genotype
{
  friend class cereal::access;

public:
  uint32_t id;
  uint16_t num;
  uint32_t first_variant_node;

  Genotype();
  Genotype(uint32_t i, uint16_t n, uint32_t fvn);

private:
  template <class Archive>
  void serialize(Archive & ar, unsigned int version);
};

} // namespace gyper
