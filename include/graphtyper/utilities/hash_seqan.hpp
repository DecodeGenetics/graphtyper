#pragma once

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

  for (unsigned i = 0; i < length(s); ++i)
  {
    hash_val *= 4;
    hash_val += ordValue(s[i]);

    if (i == 32)
      return hash_val;
  }

  return hash_val;
}


inline
std::size_t
hash_value_optional1(String<Dna> const & s)
{
  std::size_t hash_val = 0;

  for (unsigned i = 0; i < length(s); ++i)
  {
    hash_val = hash_val * 4 + ordValue(s[i]);

    if (i == 32)
      return hash_val;
  }

  return hash_val;
}


inline
std::size_t
hash_value_optional2(String<Dna> const & s)
{
  std::size_t hash_val = 0;

  for (Iterator<String<Dna> const>::Type it = begin(s); it != end(s); ++it)
  {
    hash_val = hash_val * 4 + ordValue(*it);
  }

  return hash_val;
}


inline
std::size_t
hash_value_optional3(String<Dna> const & s)
{
  std::size_t hash_val = 0;
  unsigned i = 0;

  for (Iterator<String<Dna> const>::Type it = begin(s); it != end(s); ++it, ++i)
  {
    hash_val = hash_val * 4 + ordValue(*it);

    if (i == 32)
    {
      break;
    }
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
