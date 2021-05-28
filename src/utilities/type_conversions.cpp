#include <bitset> // std::bitset
#include <iomanip> // std::hex

#include <graphtyper/utilities/type_conversions.hpp>

#include <graphtyper/utilities/logging.hpp>

#include <seqan/basic.h>
#include <seqan/sequence.h>

namespace gyper
{

/**
 * @brief Converts a string/list of DNA bases to a unsigned 64 bit integer.
 *
 * @param s String of DNA bases.
 * @return The new unsigned 64 bit integer.
 */
uint64_t
to_uint64(char const c)
{
  switch (c)
  {
  case 'A':
    return 0ul;

  case 'C':
    return 1ul;

  case 'G':
    return 2ul;

  case 'T':
    return 3ul;

  default:
    print_log(log_severity::warning, "[graphtyper::type_conversions] Invalid character ", c);
    std::exit(1);
    return 0ul;
  }
}


/*
uint64_t
to_uint64(seqan::DnaString const & s, std::size_t i)
{
  SEQAN_ASSERT_MSG(seqan::length(s) - i >= 32, "Cannot read 32 bases from read!");
  uint64_t d = 0;

  for (std::size_t const j = i + 32; i < j; ++i)
  {
    d <<= 2;
    d += seqan::ordValue(s[i]);
  }

  return d;
}
*/


uint64_t
to_uint64(std::vector<char> const & s, int i)
{
  assert(s.size() >= 32u);
  uint64_t d = 0ull;

  for (int const j = i + 32; i < j; ++i)
  {
    d <<= 2;
    d += to_uint64(s[i]);
  }

  return d;
}


uint64_t
to_uint64(std::string const & s, int i)
{
  assert(s.size() >= 32u);
  uint64_t d = 0ull;

  for (int const j = i + 32; i < j; ++i)
  {
    d <<= 2;
    d += to_uint64(s[i]);
  }

  return d;
}


/** WIP
uint16_t
to_uint16(char const c)
{
  switch (c)
  {
  case 'a':
  case 'A':
    return 0u;

  case 'c':
  case 'C':
    return 1u;

  case 'g':
  case 'G':
    return 2u;

  case 't':
  case 'T':
    return 3u;

  default:
    print_log(log_severity::warning, "[graphtyper::type_conversions] Invalid character ", c);
    return 4ul;
  }
}


uint16_t
to_uint16(std::vector<char> const & s, std::size_t i)
{
  assert(s.size() >= 8 + i);
  uint16_t d = 0u;

  for (std::size_t const j = i + 8; i < j; ++i)
  {
    d <<= 2;
    d += to_uint16(s[i]);
  }

  return d;
}



int32_t
to_uint16_with_check(std::vector<char> const & s, long i)
{
  assert(static_cast<long>(s.size()) >= 8 + i);
  int32_t d{0};

  for (long const j = i + 8; i < j; ++i)
  {
    d <<= 2;
    uint16_t ret = to_uint16(s[i]);

    if (ret == 4)
      return -1;
    else
      d += ret;
  }

  return d;
}


std::vector<uint16_t>
to_uint16_vec(std::string const & s)
{
  assert(s.size() >= 8); // Otherwise its impossible to read 8 bases from read!
  uint16_t current{0u}
  std::vector<uint16_t> uints;
  long const sequence_size = s.size();
  long const max_8mers{sequence_size - 7};
  uints.reserve(max_8mers);
  int bp_since_N{0};

  for (long i = 0; i < max_8mers; ++i)
  {
    char const c = s[i];
    ++bp_since_N;

    switch (c)
    {
    case 'a':
    case 'A':
      current *= 4;
      //current += 0;
      break;

    case 'c':
    case 'C':
      current *= 4;
      current += 1;

    case 'g':
    case 'G':
      current *= 4;
      current += 2;
      break;

    case 't':
    case 'T':
      current *= 4;
      current += 3;
      break;

    default:
      bp_since_N = 0;
      break;
    }

    if (bp_since_N >= 8)
      uints.push_back(current);
  }

  return uints;
}
*/


template <typename TSeq>
std::vector<uint64_t>
to_uint64_vec(TSeq const & s, std::size_t i)
{
  assert(seqan::length(s) >= 32 + i);  // Cannot read 32 bases from read!"
  std::vector<uint64_t> uints(1, 0u);

  for (unsigned const j = i + 32; i < j; ++i)
  {
    std::size_t const origin_size = uints.size();

    if (origin_size > 97)
      return std::vector<uint64_t>();

    for (std::size_t u = 0; u < origin_size; ++u)
    {
      std::bitset<4> const iupac(seqan::ordValue(s[i]));

      if (iupac.all() || iupac.none())
      {
        uints.push_back(uints[u] * 4 + 0); // A
        uints.push_back(uints[u] * 4 + 1); // C
        uints.push_back(uints[u] * 4 + 2); // G
        uints[u] <<= 2; uints[u] += 3;     // T
      }
      else
      {
        std::size_t set_count = iupac.count();

        auto check_set_count =
          [&uints, u, &set_count](std::size_t const to_add)
          {
            assert(set_count != 0);

            if (set_count == 1)
            {
              uints[u] *= 4;
              uints[u] += to_add;
            }
            else
            {
              uints.push_back(uints[u] * 4 + to_add);
            }

            --set_count;
          };

        if (iupac.test(0))
          check_set_count(0); // A
        if (iupac.test(1))
          check_set_count(1); // C
        if (iupac.test(2))
          check_set_count(2); // G
        if (iupac.test(3))
          check_set_count(3); // T
      }
    }
  }

  return uints;
}


// Explicit instantation
template std::vector<uint64_t> to_uint64_vec<seqan::Dna5String>(seqan::Dna5String const & s, std::size_t i);
template std::vector<uint64_t> to_uint64_vec<seqan::IupacString>(seqan::IupacString const & s, std::size_t i);


std::array<uint64_t, 96>
to_uint64_vec_hamming_distance_1(uint64_t const key)
{
  std::array<uint64_t, 96> hamming1;

  for (uint8_t bb = 0; bb < 32; ++bb)
  {
    // We add mask ^ exact kmer
    // E.g. with bb = 2
    // Mask for flipping both bits will be
    // 0x0000000000000030
    hamming1[bb * 3 + 0] = (1ull << (bb * 2)) ^ key; // Flip second bit
    hamming1[bb * 3 + 1] = (2ull << (bb * 2)) ^ key; // Flip first bit
    hamming1[bb * 3 + 2] = (3ull << (bb * 2)) ^ key; // Flip both bits
  }

  return hamming1;
}


/*
seqan::DnaString
to_dna(uint64_t const & d, uint8_t k)
{
  seqan::String<seqan::Dna> new_dna_string = "";

  while (k > 0)
  {
    --k;

    switch ((d & (0x0000000000000003ULL << 2 * k)) >> 2 * k)
    {
    case A_VALUE:
      seqan::append(new_dna_string, seqan::Dna('A'));
      break;

    case C_VALUE:
      seqan::append(new_dna_string, seqan::Dna('C'));
      break;

    case G_VALUE:
      seqan::append(new_dna_string, seqan::Dna('G'));
      break;

    case T_VALUE:
      seqan::append(new_dna_string, seqan::Dna('T'));
      break;
    }
  }

  return new_dna_string;
}
*/

std::string
to_dna_str(uint64_t const d, uint8_t k)
{
  std::string new_dna;

  while (k > 0)
  {
    --k;

    switch ((d & (0x0000000000000003ULL << 2 * k)) >> 2 * k)
    {
    case A_VALUE:
      new_dna.push_back('A');
      break;

    case C_VALUE:
      new_dna.push_back('C');
      break;

    case G_VALUE:
      new_dna.push_back('G');
      break;

    case T_VALUE:
      new_dna.push_back('T');
      break;
    }
  }

  return new_dna;
}


std::array<uint64_t, 3>
get_mismatches_of_last_base(uint64_t const d)
{
  uint64_t const d2 = ((d >> 2) << 2);

  switch (d & 0x0000000000000003ULL)
  {
  case A_VALUE: return {{
      d2 | C_VALUE, d2 | G_VALUE, d2 | T_VALUE
    }};

  case C_VALUE: return {{
      d2 | A_VALUE, d2 | G_VALUE, d2 | T_VALUE
    }};

  case G_VALUE: return {{
      d2 | A_VALUE, d2 | C_VALUE, d2 | T_VALUE
    }};

  default: /*T*/ return {{
      d2 | A_VALUE, d2 | C_VALUE, d2 | G_VALUE
    }};
  }
}


std::array<uint64_t, 3>
get_mismatches_of_first_base(uint64_t const d)
{
  uint64_t const d2 = ((d << 2) >> 2);

  switch (d & 0xC000000000000000ULL)
  {
  case A_LAST_VALUE: return {{
      d2 | C_LAST_VALUE, d2 | G_LAST_VALUE, d2 | T_LAST_VALUE
    }};

  case C_LAST_VALUE: return {{
      d2 | A_LAST_VALUE, d2 | G_LAST_VALUE, d2 | T_LAST_VALUE
    }};

  case G_LAST_VALUE: return {{
      d2 | A_LAST_VALUE, d2 | C_LAST_VALUE, d2 | T_LAST_VALUE
    }};

  default: /*T*/ return {{
      d2 | A_LAST_VALUE, d2 | C_LAST_VALUE, d2 | G_LAST_VALUE
    }};
  }
}


/**
 * @brief Inserts all elements from one map to another and optionally deletes the elements from the old map.
 */
template <typename Map1, typename Map2>
void
map_swap(Map1 & map1, Map2 & map2)
{
  for (auto it = map1.begin(); it != map1.end(); ++it)
  {
    map2[it->first] = it->second;
    // it = kmers->erase(it);
  }
}


template <typename Map1, typename Map2>
void
map_swap(Map1 & map1, Map2 & map2, bool const & delete_as_i_go)
{
  if (!delete_as_i_go)
  {
    return map_swap(map1, map2);
  }

  for (auto it = map1.begin(); it != map1.end(); ++it)
  {
    map2[it->first] = it->second;
    it = map1.erase(it);
  }
}


template <class TListItems>
std::vector<TListItems>
to_vector(std::forward_list<TListItems> const & q)
{
  std::vector<TListItems> new_list;

  for (auto q_it = q.cbegin(); q_it != q.cend(); ++q_it)
  {
    new_list.push_back(*q_it);
  }

  return new_list;
}


std::vector<char>
to_vec(std::string && str)
{
  std::vector<char> vec(str.size());

  for (unsigned i = 0; i < str.size(); ++i)
  {
    vec[i] = str[i];
  }

  return vec;
}


} // namespace gyper
