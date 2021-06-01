#pragma once

#include <initializer_list>
//#include <set>
#include <unordered_set>
#include <vector>

namespace gyper
{
class Alt
{
public:
  // int32_t non_merge_start{-1};
  // int32_t non_merge_end{-1};
  std::vector<char> seq;
  std::unordered_set<long> events;
  std::unordered_set<long> anti_events;

  Alt() = default;
  Alt(Alt const &) = default;
  Alt(Alt &&) = default;
  ~Alt() = default;
  Alt(std::vector<char> && _seq);
  Alt(std::vector<char> const & _seq);
  Alt(std::initializer_list<char> l);

  void clear();

  bool operator==(Alt const & o);
  bool operator!=(Alt const & o);
  bool operator<(Alt const & o);

  Alt & operator=(Alt const &) = default;
  Alt & operator=(Alt &&) = default;
};

bool compare_two_alts(Alt const & lhs, Alt const & rhs);
bool is_empty_seq(Alt const & alt);
Alt make_alt(Alt const & prev, Alt const & curr, long const jump_size);
bool is_ok_to_merge_alts(Alt const & prev_alt, Alt const & curr_alt);
} // namespace gyper
