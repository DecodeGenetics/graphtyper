#pragma once

#include <cassert>
#include <cstdint> // uint32_t
#include <utility> // std::pair


namespace gyper
{

std::pair<uint16_t, uint16_t> inline
to_pair(uint32_t const index)
{
  // index = x + (y+1)*y/2 = x + (y^2 + y)/2
  uint32_t y = 1u;

  while ((y * y + y) / 2u <= index)
    ++y;

  --y;
  // uint32_t x = index - (y*y + y)/2u;
  return std::make_pair<uint16_t, uint16_t>(index - (y * y + y) / 2u, std::move(y));
}


uint32_t inline
to_index(uint32_t const x, uint32_t const y)
{
  assert(x <= y);
  return x + (y + 1u) * y / 2u;
}


} // namespace gyper
