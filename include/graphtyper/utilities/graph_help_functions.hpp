#pragma once

#include <cassert>
#include <cstdint> // uint32_t
#include <utility> // std::pair


namespace gyper
{

std::pair<uint16_t, uint16_t> inline
to_pair(long const index)
{
  // index = x + (y+1)*y/2 = x + (y^2 + y)/2
  long y = 1;

  while ((y * y + y) / 2 <= index)
    ++y;

  --y;
  // uint32_t x = index - (y*y + y)/2u;
  return std::make_pair<uint16_t, uint16_t>(index - (y * y + y) / 2u, static_cast<uint16_t>(y));
}


long inline
to_index(long const x, long const y)
{
  assert(x <= y);
  return x + (y + 1) * y / 2;
}


} // namespace gyper
