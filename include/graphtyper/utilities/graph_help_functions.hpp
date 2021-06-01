#pragma once

#include <cassert>
#include <cstdint> // uint32_t
#include <utility> // std::pair

namespace gyper
{
std::pair<uint16_t, uint16_t> inline to_pair(long const index)
{
  // index = x + (y+1)*y/2 = x + (y^2 + y)/2
  long y = 1;

  while ((y * y + y) / 2 <= index)
    ++y;

  --y;
  return std::make_pair<uint16_t, uint16_t>(static_cast<uint16_t>(index - (y * y + y) / 2u), static_cast<uint16_t>(y));
}

long inline to_index(long const x, long const y)
{
  assert(x <= y);
  return x + (y + 1) * y / 2;
}

long inline to_index_safe(long const x, long const y)
{
  if (x <= y)
    return x + (y + 1) * y / 2;
  else
    return y + (x + 1) * x / 2;
}

} // namespace gyper
