#pragma once

#include <cstdint>

namespace gyper
{
class Location
{
public:
  char node_type;
  uint32_t node_index;
  uint32_t node_order;
  uint32_t offset;

  bool inline is_available() const
  {
    return node_type != 'U';
  }

  bool inline is_unavailable() const
  {
    return node_type == 'U';
  }

  Location() : node_type('U'), node_index(0), node_order(0), offset(0)
  {
  }

  Location(char && _node_type, uint32_t _node_index, uint32_t _node_order, uint32_t _offset) :
    node_type(std::move(_node_type)), node_index(_node_index), node_order(_node_order), offset(_offset)
  {
  }
};

} // namespace gyper
