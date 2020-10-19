#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include <graphtyper/graph/read_strand.hpp>
#include <graphtyper/typer/event.hpp>


namespace gyper
{

class EventCall
{
public:
  IndelEvent event; // includes pos, type and sequence

  uint32_t count{0};
  uint32_t anti_count{0};
  uint32_t multi_count{0};

  uint16_t span{0};
};

//class EventCalls
//{
//public:
//
//};

} // namespace gyper
