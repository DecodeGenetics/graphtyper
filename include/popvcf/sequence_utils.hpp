#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <charconv>
#include <string>
#include <string_view>
#include <vector>

namespace popvcf
{
uint32_t constexpr CHAR_SET_SIZE = 69;
uint32_t constexpr CHAR_SET_SIZE_2BYTES = CHAR_SET_SIZE * CHAR_SET_SIZE;
char constexpr CHAR_SET_MIN = ':';

long constexpr ENC_BUFFER_SIZE{4 * 65536}; //!< Buffer size of arrays when encoding
long constexpr DEC_BUFFER_SIZE{8 * 65536}; //!< Buffer size of arrays when decoding

//! Data type of an encoding array buffer
using Tenc_array_buf = std::array<char, ENC_BUFFER_SIZE>;

//! Data type of a decoding array buffer
using Tdec_array_buf = std::array<char, DEC_BUFFER_SIZE>;

inline char int_to_ascii(uint32_t in)
{
  assert(in < CHAR_SET_SIZE);

  return CHAR_SET_MIN + in;
}

inline uint32_t ascii_to_int(char in)
{
  assert(in >= CHAR_SET_MIN);

  return static_cast<uint32_t>(in) - CHAR_SET_MIN;
}

inline std::string int_to_ascii_string(uint32_t in)
{
  std::string str;

  while (in >= CHAR_SET_SIZE)
  {
    uint32_t rem = in % CHAR_SET_SIZE;
    in = in / CHAR_SET_SIZE;
    str.push_back(int_to_ascii(rem));
  }

  assert(in < CHAR_SET_SIZE);
  str.push_back(int_to_ascii(in));
  return str;
}

inline uint32_t ascii_string_view_to_int(std::string_view in)
{
  uint32_t const in_size = in.size();
  assert(in_size > 0);
  uint32_t out{0};
  uint32_t pow{1};

  for (uint32_t c{0}; c < in_size; ++c)
  {
    out += pow * ascii_to_int(in[c]);
    pow *= CHAR_SET_SIZE;
  }

  return out;
}

inline uint32_t ascii_cstring_to_int(char const * b, char const * e)
{
  uint32_t out{ascii_to_int(*b)};
  ++b;
  uint32_t pow{CHAR_SET_SIZE};

  while (b != e)
  {
    out += pow * ascii_to_int(*b);
    pow *= CHAR_SET_SIZE;
    ++b;
  }

  return out;
}

template <typename Tint, typename Tbuffer_out>
inline void to_chars(Tint char_val, Tbuffer_out & buffer_out)
{
  std::size_t constexpr ARR_SIZE{6};
  std::array<char, ARR_SIZE> a;
  std::size_t i{0};

  while (char_val >= CHAR_SET_SIZE)
  {
    Tint rem = char_val % CHAR_SET_SIZE;
    char_val = char_val / CHAR_SET_SIZE;
    assert(i < ARR_SIZE);
    a[i++] = int_to_ascii(rem);
  }

  assert(char_val < CHAR_SET_SIZE);
  assert(i < ARR_SIZE);
  a[i++] = int_to_ascii(char_val);
  buffer_out.insert(buffer_out.end(), a.data(), a.data() + i);
}

template <typename Tstring>
std::vector<std::string_view> split_string(Tstring const & str, char const delimiter);

template <typename Tit>
long get_vcf_pos(Tit begin, Tit end)
{
  auto find_it1 = std::find(begin, end, '\t');
  auto find_it2 = std::find(find_it1 + 1, end, '\t');
  long vcf_pos{0};
  std::from_chars(find_it1 + 1, find_it2, vcf_pos);
  return vcf_pos;
}

template <typename Tbuffer_in>
inline void resize_input_buffer(Tbuffer_in & buffer_in, std::size_t const new_size)
{
  buffer_in.resize(new_size);
}

template <>
inline void resize_input_buffer(Tdec_array_buf & /*buffer_in*/, std::size_t const /*new_size*/)
{
  // Do nothing. Arrays are not resized
}

template <>
inline void resize_input_buffer(Tenc_array_buf & /*buffer_in*/, std::size_t const /*new_size*/)
{
  // Do nothing. Arrays are not resized
}

} // namespace popvcf
