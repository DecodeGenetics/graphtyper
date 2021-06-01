#include <vector>

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/ref.hpp>

namespace gyper
{
Ref::Ref(std::vector<char> && _seq) : seq(std::forward<std::vector<char>>(_seq))
{
}

Ref::Ref(std::vector<char> const & _seq) : seq(_seq)
{
}

Ref::Ref(std::initializer_list<char> l) : seq(l)
{
}

void Ref::clear()
{
  seq.clear();
  events.clear();
  // anti_events.clear();
}

bool Ref::operator==(Ref const & o)
{
  return seq == o.seq;
}

bool Ref::operator!=(Ref const & o)
{
  return !(seq == o.seq);
}

bool Ref::operator<(Ref const & o)
{
  return seq < o.seq;
}

} // namespace gyper
