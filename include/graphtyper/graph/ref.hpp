#pragma once

#include <initializer_list>
#include <unordered_set>
#include <vector>

namespace gyper
{

class Ref
{
public:
  std::vector<char> seq;
  std::unordered_set<long> events;
  std::unordered_set<long> anti_events;

  Ref() = default;
  Ref(std::vector<char> && _seq);
  Ref(std::vector<char> const & _seq);
  Ref(std::initializer_list<char> l);

  void clear();

  bool operator==(Ref const & o);
  bool operator!=(Ref const & o);
  bool operator<(Ref const & o);

};

} // namespace gyper
