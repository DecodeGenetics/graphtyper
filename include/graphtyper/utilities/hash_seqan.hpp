#pragma once

#include <seqan/sequence.h>

namespace seqan
{
std::size_t hash_value(String<char> const & s);
std::size_t hash_value(String<Dna> const & s);

} // namespace seqan

namespace std
{
template <>
struct hash<seqan::String<seqan::Dna>>
{
  std::size_t operator()(seqan::String<seqan::Dna> const & s) const;
};

template <>
struct hash<seqan::String<char>>
{
  std::size_t operator()(seqan::String<char> const & s) const;
};

} // namespace std
