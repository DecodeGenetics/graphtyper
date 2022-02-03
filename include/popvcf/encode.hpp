#pragma once

#include <charconv>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include <parallel_hashmap/phmap.h>

#include "sequence_utils.hpp"

namespace popvcf
{
//! Buffer size when encoding
long constexpr ENC_BUFFER_SIZE{4 * 65536};

//! Data type of an array buffer
using Tarray_buf = std::array<char, ENC_BUFFER_SIZE>; //!< Buffer type

class EncodeData
{
public:
  std::size_t bytes_read{0};
  std::size_t field{0};         // current vcf field
  std::size_t b{0};             // begin index in buffer_in
  std::size_t i{b};             // index in buffer_in
  bool header_line{true};       //!< True iff in header line
  bool no_previous_line{false}; //! Set to skip using previous line

  std::string prev_contig{};
  int64_t prev_pos{0};
  uint32_t prev_num_unique_fields{};
  std::vector<std::string> prev_unique_fields{};
  std::vector<uint32_t> prev_field2uid{};
  phmap::flat_hash_map<std::string, uint32_t> prev_map_to_unique_fields{};

  std::string contig{};
  int64_t pos{0};
  uint32_t num_unique_fields{};
  std::vector<std::string> unique_fields{};
  std::vector<uint32_t> field2uid{};
  phmap::flat_hash_map<std::string, uint32_t> map_to_unique_fields{};

  inline void clear_line()
  {
    field = 0; // reset field index

    if (not no_previous_line)
    {
      std::swap(prev_contig, contig);
      prev_pos = pos;
      prev_num_unique_fields = num_unique_fields;
      std::swap(prev_unique_fields, unique_fields);
      std::swap(prev_field2uid, field2uid);
      std::swap(prev_map_to_unique_fields, map_to_unique_fields);
    }

    contig.resize(0);
    pos = 0;
    num_unique_fields = 0;
    unique_fields.resize(0);
    field2uid.resize(0);
    map_to_unique_fields.clear(); // clear map every line
  }
};

template <typename Tint, typename Tbuffer_out>
inline void to_chars(Tint char_val, Tbuffer_out & buffer_out)
{
  while (char_val >= CHAR_SET_SIZE)
  {
    auto rem = char_val % CHAR_SET_SIZE;
    char_val = char_val / CHAR_SET_SIZE;
    buffer_out.push_back(int_to_ascii(rem));
  }

  assert(char_val < CHAR_SET_SIZE);
  buffer_out.push_back(int_to_ascii(char_val));
}

//! Encodes an input buffer. Output is written in \a buffer_out.
template <typename Tbuffer_out, typename Tbuffer_in>
inline void encode_buffer(Tbuffer_out & buffer_out, Tbuffer_in & buffer_in, EncodeData & ed)
{
  buffer_out.reserve(ENC_BUFFER_SIZE / 2);
  std::size_t constexpr N_FIELDS_SITE_DATA{9}; // how many fields of the VCF contains site data

  while (ed.i < ed.bytes_read)
  {
    char const b_in = buffer_in[ed.i];

    if (b_in != '\t' && b_in != '\n')
    {
      ++ed.i;
      continue; // we are in a vcf field
    }

    if (ed.field == 0)
    {
      // check if in header line and store contig
      if (buffer_in[ed.b] == '#')
      {
        ed.header_line = true;
      }
      else
      {
        ed.header_line = false;
        ed.contig = std::string(&buffer_in[ed.b], ed.i - ed.b);
      }
    }
    else if (ed.header_line == false && ed.field == 1)
    {
      std::from_chars(&buffer_in[ed.b], &buffer_in[ed.i], ed.pos);

      if (ed.contig != ed.prev_contig || (ed.pos / 10000) != (ed.prev_pos / 10000))
      {
        // previous line is not available
        ed.prev_num_unique_fields = 0;
        ed.prev_map_to_unique_fields.clear();
      }
    }

    if (ed.header_line || ed.field < N_FIELDS_SITE_DATA)
    {
      ++ed.i; // adds '\t' or '\n'
      std::copy(&buffer_in[ed.b], &buffer_in[ed.i], std::back_inserter(buffer_out));
    }
    else
    {
      assert(buffer_in[ed.b] >= '!');
      assert(buffer_in[ed.b] <= '9');

      // check if it is in the current line
      auto insert_it = ed.map_to_unique_fields.insert(
        std::pair<std::string, uint32_t>(std::piecewise_construct,
                                         std::forward_as_tuple(&buffer_in[ed.b], ed.i - ed.b),
                                         std::forward_as_tuple(ed.num_unique_fields)));

      long const field_idx = ed.field - N_FIELDS_SITE_DATA;
      assert(field_idx == static_cast<long>(ed.field2uid.size()));

      if (insert_it.second == true)
      {
        ed.field2uid.push_back(ed.unique_fields.size());
        ++ed.num_unique_fields; // unique field
        ed.unique_fields.emplace_back(&buffer_in[ed.b], ed.i - ed.b);

        assert(ed.num_unique_fields == static_cast<long>(ed.unique_fields.size()));

        if (field_idx < static_cast<long>(ed.prev_field2uid.size()) &&
            ed.prev_unique_fields[ed.prev_field2uid[field_idx]] == ed.unique_fields[insert_it.first->second])
        {
          /* Case 0: unique and same as above. */
          buffer_out.push_back('$');
          buffer_out.push_back(buffer_in[ed.i]);
          ++ed.i;
        }
        else
        {
          // check if it is in the previous line
          auto prev_find_it = ed.prev_map_to_unique_fields.find(insert_it.first->first);

          if (prev_find_it == ed.prev_map_to_unique_fields.end())
          {
            /* Case 1: Field is unique in the current line and is not in the previous line. */
            ++ed.i;                                                                        // adds '\t' or '\n'
            std::copy(&buffer_in[ed.b], &buffer_in[ed.i], std::back_inserter(buffer_out)); // just copy as is
          }
          else
          {
            /* Case 2: Field is unique in the current line but identical to a field in the previous line. */
            buffer_out.push_back('%');
            popvcf::to_chars(prev_find_it->second, buffer_out);
            buffer_out.push_back(buffer_in[ed.i]); // write '\t' or '\n'
            ++ed.i;
          }
        }
      }
      else
      {
        ed.field2uid.push_back(insert_it.first->second);

        if (field_idx < static_cast<long>(ed.prev_field2uid.size()) &&
            ed.prev_unique_fields[ed.prev_field2uid[field_idx]] == ed.unique_fields[insert_it.first->second])
        {
          /* Case 3: Field is not unique and same has the field above. */
          buffer_out.push_back('&');
          buffer_out.push_back(buffer_in[ed.i]);
          ++ed.i;
        }
        else
        {
          /* Case 4: Field is a duplicate in the current line. */
          popvcf::to_chars(insert_it.first->second, buffer_out);
          buffer_out.push_back(buffer_in[ed.i]); // write '\t' or '\n'
          ++ed.i;
        }
      }

      assert((field_idx + 1) == static_cast<long>(ed.field2uid.size()));
      assert(ed.field2uid[0] == 0);
    }

    assert(b_in == buffer_in[ed.i - 1]); // i should have been already incremented here
    ed.b = ed.i;                         // set begin index of next field

    // check if we need to clear line or increment field
    if (b_in == '\n')
      ed.clear_line();
    else
      ++ed.field;
  } // ends inner loop

  // copy the remaining data to the beginning of the input buffer
  std::copy(&buffer_in[ed.b], &buffer_in[ed.i], &buffer_in[0]);
  ed.i = ed.i - ed.b;
  ed.b = 0;
  ed.bytes_read = ed.i;
}

//! Encode a gzipped file and write to stdout
void encode_file(std::string const & input_fn,
                 bool const is_bgzf_input,
                 std::string const & output_fn,
                 std::string const & output_mode,
                 bool const is_bgzf_output,
                 int const compression_threads,
                 bool const no_previous_line);

} // namespace popvcf
