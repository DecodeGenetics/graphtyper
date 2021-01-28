#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <sstream>
#include <map>
#include <unordered_map>
#include <vector>


namespace gyper
{

std::array<char, 4> constexpr index2base = {'A', 'C', 'G', 'T'};

long
base2index(char const base);


struct BaseCount
{
  std::array<int32_t, 4> acgt;
  std::array<int64_t, 4> acgt_qualsum;
  int32_t deleted{0};
  int32_t unknown{0};

  long get_depth_without_deleted() const;
  long get_depth_with_deleted() const;
  std::string to_string() const;

  void add_base(char seq, char qual);

};


class Event
{
public:
  // Event type is the same as the CIGAR mappings, that is:
  //   I for insertion
  //   D for deletion
  //   X for mismatches
  // other cigar mappings are not used in this context
  uint32_t pos{0}; // event is between pos-1 and pos in a 1-based system
  char type{'\0'};
  std::vector<char> sequence{};

  Event() = default;

  Event(uint32_t _pos, char _type)
    : pos(_pos)
    , type(_type)
  {}

  Event(uint32_t _pos, char _type, std::vector<char> && _sequence)
    : pos{_pos}
    , type{_type}
    , sequence{std::forward<std::vector<char> >(_sequence)}
  {}

  inline std::string
  to_string() const
  {
    std::ostringstream ss;
    ss << pos << " " << type << ' ' << std::string(sequence.begin(), sequence.end());
    return ss.str();
  }


  /*********************
   * OPERATOR OVERLOAD *
   *********************/
  bool operator==(Event const & b) const;
  bool operator!=(Event const & b) const;
  bool operator<(Event const & b) const;
};


class EventSupport
{
public:
  uint16_t hq_count{0};
  uint16_t lq_count{0};
  uint16_t proper_pairs{0};
  uint16_t first_in_pairs{0};
  uint16_t sequence_reversed{0u};
  uint16_t clipped{0};
  uint8_t max_mapq{0};
  uint8_t max_distance{0};
  int32_t uniq_pos1{-1};
  int32_t uniq_pos2{-1};
  int32_t uniq_pos3{-1};
  std::map<Event, uint16_t> phase;

  void clear();
  int get_raw_support() const;
  double corrected_support() const;
  bool has_good_support(long const cov) const;
  std::string to_string() const;
  uint32_t log_qual(uint32_t eps = 7) const;
  bool is_good_indel(uint32_t eps = 7) const;

  // indel things
  uint16_t multi_count{0};
  uint16_t anti_count{0};
  uint16_t span{1};

  bool has_realignment_support{false};
  bool has_indel_good_support{false};

  uint32_t max_log_qual{0};
  int max_log_qual_file_i{-1};


};


//std::vector<uint8_t> get_phred_biallelic(uint32_t count, uint32_t anti_count, uint32_t eps);
uint32_t get_log_qual(uint32_t count, uint32_t anti_count, uint32_t eps = 7);
uint32_t get_log_qual_double(double count, double anti_count, double eps = 7.0);

/*
class EventInfo
{
public:
  uint32_t count{0};
  uint32_t multi_count{0};
  uint32_t anti_count{0};
  uint16_t span{1};

  // indels with realignment support are good enough to try to realign reads to them
  bool has_realignment_support{false};
  bool has_good_support{false};

  uint32_t max_log_qual{0};
  int max_log_qual_file_i{-1};

  uint32_t log_qual(uint32_t eps = 7) const;
  void clean_counts();
};
*/


bool
apply_indel_event(std::vector<char> & sequence,
                  std::vector<int32_t> & ref_pos,
                  Event const & indel_event,
                  long const offset,
                  bool const is_debug = false);

Event
make_deletion_event(std::vector<char> const & reference_sequence, long ref_offset, int32_t pos, long count);

Event
make_insertion_event(int32_t pos, std::vector<char> && event_sequence);

} // namespace gyper


namespace std
{

template <>
struct hash<gyper::Event>
{
  size_t
  operator()(gyper::Event const & event)
  {
    size_t h1 = std::hash<uint32_t>()(event.pos);
    size_t h2 = std::hash<char>()(event.type);
    size_t h3 = 42 ^ event.sequence.size();

    for (auto const & seq : event.sequence)
      h3 ^= std::hash<char>()(seq);

    return h1 ^ (h2 << 1) ^ (h3 + 0x9e3779b9);
  }


};

} // namespace std


namespace gyper
{

//using Tevents = phmap::node_hash_map<Event, uint32_t>;
using Tindel_events = std::map<Event, EventSupport>; // maps events to count

} // namespace gyper
