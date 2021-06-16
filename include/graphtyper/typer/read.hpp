#pragma once

#include <cstdint>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include <parallel_hashmap/phmap_fwd_decl.h>

#include <graphtyper/typer/event.hpp>

namespace gyper
{
int32_t static constexpr READ_ANTI_SUPPORT{-1};
int32_t static constexpr READ_MULTI_SUPPORT{-2};

struct ReadIndelEvent
{
public:
  int32_t read_pos{0};
  Tindel_events::iterator event_it;

  ReadIndelEvent(long _read_pos, Tindel_events::iterator _event_it) : read_pos(_read_pos), event_it(_event_it)
  {
  }
};

struct ReadSnpEvent
{
public:
  int32_t read_pos{0};
  int32_t ref_pos{0};

  ReadSnpEvent(int32_t _read_pos, int32_t ref_pos) : read_pos(_read_pos), ref_pos(ref_pos)
  {
  }
};

class Alignment
{
public:
  int32_t pos{-1};
  int32_t pos_end{-1};
  int32_t score{std::numeric_limits<int32_t>::min()}; // Alignment score
  int16_t num_clipped_begin{0};                       // Number of clipped bases at the begin of the read
  int16_t num_clipped_end{0};                         // Number of clipped bases at the end of the read
  int16_t num_ins_begin{0};                           // Number of bases in insertion at the begin of the read
  // int16_t num_ins_end{0}; // Number of bases in insertion at the end of the read
  // std::vector<ReadSnpEvent> snp_events; // read pos containing snp events
  // std::vector<std::map<SnpEvent, SnpEventInfo>::iterator> snp_events;
  std::vector<ReadIndelEvent> indel_events; // pair of read_offset and event

  bool has_indel_event(Tindel_events::iterator) const;
  bool is_clipped() const;
  void add_indel_event(int32_t pos, uint16_t flags, uint8_t mapq, Tindel_events::iterator indel_it);
  void replace_indel_events(uint16_t flags, uint8_t mapq, std::vector<ReadIndelEvent> && new_indel_events);
};

class Read
{
public:
  std::string name;
  int32_t mate_pos{-1};
  uint16_t flags{0};
  uint8_t mapq{255};

  std::vector<char> sequence{};
  std::vector<char> qual{};
  Alignment alignment{};

  std::string to_string() const;
};

} // namespace gyper
