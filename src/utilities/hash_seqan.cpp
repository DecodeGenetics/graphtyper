#include <algorithm>
#include <string_view>

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <graphtyper/utilities/hash_seqan.hpp>

namespace std
{

std::size_t
hash<seqan::String<seqan::Dna> >::operator()(seqan::String<seqan::Dna> const & s) const
{
  std::size_t hash_val = 0;

  for (seqan::Iterator<seqan::String<seqan::Dna> const>::Type it = seqan::begin(s); it != seqan::end(s); ++it)
  {
    hash_val = hash_val * 4 + seqan::ordValue(*it);
  }

  return hash_val;
}


std::size_t
hash<seqan::String<char> >::operator()(seqan::String<char> const & s) const
{
  return std::hash<std::string_view>{}(std::string_view{&*seqan::begin(s), seqan::length(s)});
}


} // namespace std
