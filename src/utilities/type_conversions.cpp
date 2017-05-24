#include <bitset> // std::bitset
#include <iomanip> // std::hex
#include <iostream> // std::cerr, std::endl

#include <graphtyper/utilities/type_conversions.hpp>


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
  switch(c)
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
      std::cerr << "[graphtyper::type_conversions] WARNING: Invalid character " << c << std::endl;
      return 0ul;
  }
}


uint64_t
to_uint64(seqan::DnaString const & s, std::size_t i)
{
  SEQAN_ASSERT_MSG(seqan::length(s) - i >= gyper::K, "Cannot read 32 bases from read!");
  uint64_t d = 0;

  for (std::size_t const j = i + K; i < j; ++i)
  {
    d *= 4;
    d += seqan::ordValue(s[i]);
  }

  return d;
}


uint64_t
to_uint64(std::vector<char> const & s)
{
  assert(s.size() == 32u);
  uint64_t d = 0;

  for (unsigned i = 0; i < 32; ++i)
  {
    d *= 4u;
    d += to_uint64(s[i]);
  }

  return d;
}


template <typename TSeq>
std::vector<uint64_t>
to_uint64_vec(TSeq const & s, std::size_t i)
{
  assert(seqan::length(s) - i >= gyper::K);  // Cannot read 32 bases from read!"
  std::vector<uint64_t> uints(1, 0u);

  for (unsigned const j = i + K; i < j; ++i)
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
        uints[u] *= 4;     uints[u] += 3;  // T
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
std::array<uint64_t, 2160>
to_uint64_vec_hamming_distance_2(uint64_t const key)
{
  std::array<uint64_t, 2160> hamming2;
  uint16_t i = 0;

  for (uint8_t bb = 0; bb < 16; ++bb)
  {
    // We add mask1 ^ mask2 ^ exact kmer where mask1 and mask2 do not change the same DNA base
    for (uint8_t cc = 0; cc < 16; ++cc)
    {
      if (cc == bb)
        continue;

      hamming2[i++] = (1 << (bb * 2)) ^ (1 << cc * 2) ^ key;
      hamming2[i++] = (1 << (bb * 2)) ^ (2 << cc * 2) ^ key;
      hamming2[i++] = (1 << (bb * 2)) ^ (3 << cc * 2) ^ key;
      hamming2[i++] = (2 << (bb * 2)) ^ (1 << cc * 2) ^ key;
      hamming2[i++] = (2 << (bb * 2)) ^ (2 << cc * 2) ^ key;
      hamming2[i++] = (2 << (bb * 2)) ^ (3 << cc * 2) ^ key;
      hamming2[i++] = (3 << (bb * 2)) ^ (1 << cc * 2) ^ key;
      hamming2[i++] = (3 << (bb * 2)) ^ (2 << cc * 2) ^ key;
      hamming2[i++] = (3 << (bb * 2)) ^ (3 << cc * 2) ^ key;
    }
  }

  return hamming2;
}
*/

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
