#pragma once
#include <algorithm>

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <boost/functional/hash.hpp>


namespace seqan
{

inline
std::size_t
hash_value(String<char> const & s)
{
  return boost::hash_range(begin(s), end(s));
}


inline
std::size_t
hash_value(String<Dna> const & s)
{
  std::size_t hash_val = 0;
  std::size_t const MAX_VAL = std::max(seqan::length(s), static_cast<std::size_t>(32));

  for (unsigned i = 0; i < MAX_VAL; ++i)
  {
    hash_val *= 4;
    hash_val += ordValue(s[i]);
  }

  return hash_val;
}

} // namespace seqan


namespace std
{

template <>
struct hash<seqan::String<seqan::Dna> >
{
  std::size_t
  operator()(seqan::String<seqan::Dna> const & s) const
  {
    std::size_t hash_val = 0;

    for (seqan::Iterator<seqan::String<seqan::Dna> const>::Type it = seqan::begin(s); it != seqan::end(s); ++it)
    {
      hash_val = hash_val * 4 + seqan::ordValue(*it);
    }

    return hash_val;
  }


};

template <>
struct hash<seqan::String<char> >
{
  std::size_t
  operator()(seqan::String<char> const & s) const
  {
    return boost::hash_range(seqan::begin(s), seqan::end(s));
  }


};

} // namespace std
