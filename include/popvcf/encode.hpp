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
class EncodeData
{
public:
  std::size_t field{0};   //!< current vcf field.
  std::size_t in_size{0}; //!< Size of inut buffer.
  std::size_t b{0};       //!< begin index in buffer_in
  std::size_t i{b};       //!< index in buffer_in
  bool header_line{true}; //!< True iff in header line

  /* Data fields from previous line. */
  std::vector<std::string> prev_unique_fields{};
  std::vector<uint32_t> prev_field2uid{};
  phmap::flat_hash_map<std::string, uint32_t> prev_map_to_unique_fields{};

  /* Data fields from current line. */
  std::string contig{};
  int64_t pos{0};
  int32_t stored_alt{0};
  int32_t n_alt{-1};
  std::vector<std::string> unique_fields{};
  std::vector<uint32_t> field2uid{};
  phmap::flat_hash_map<std::string, uint32_t> map_to_unique_fields{};

  /* Data fields for the next line. */
  std::string next_contig{};
  int64_t next_pos{0};

  inline void clear_line(int64_t next_pos, int32_t next_n_alt)
  {
    next_n_alt += stored_alt;
    stored_alt = 0;

    if (next_contig != contig || (next_pos / 10000) != (pos / 10000))
    {
      /// Previous line is not available, clear values
      prev_unique_fields.resize(0);
      prev_field2uid.resize(0);
      prev_map_to_unique_fields.clear();
    }
    else if (next_n_alt == n_alt)
    {
      /// Only swap out from this line if we have the same amount of alts
      std::swap(prev_unique_fields, unique_fields);
      std::swap(prev_field2uid, field2uid);
      std::swap(prev_map_to_unique_fields, map_to_unique_fields);
    }

    /// Clear data from this line for the next
    contig = next_contig;
    pos = next_pos;
    n_alt = next_n_alt;
    unique_fields.resize(0);
    field2uid.resize(0);
    map_to_unique_fields.clear();
  }
};

template <typename Tbuffer_in>
inline void set_input_size(Tbuffer_in & buffer_in, EncodeData & ed)
{
  ed.in_size = buffer_in.size();
}

template <>
inline void set_input_size(Tenc_array_buf & /*buffer_in*/, EncodeData & /*ed*/)
{
  // Do nothing.
  // NOTE: dd.in_size must be set prior to calling decode_buffer in arrays
}

//! Encodes an input buffer. Output is written in \a buffer_out.
template <typename Tbuffer_out, typename Tbuffer_in>
inline void encode_buffer(Tbuffer_out & buffer_out, Tbuffer_in & buffer_in, EncodeData & ed)
{
  set_input_size(buffer_in, ed);
  buffer_out.reserve(ENC_BUFFER_SIZE);
  std::size_t constexpr N_FIELDS_SITE_DATA{9}; // how many fields of the VCF contains site data
  int64_t next_pos{0};

  while (ed.i < ed.in_size)
  {
    char const b_in = buffer_in[ed.i];

    if (b_in != '\t' && b_in != '\n')
    {
      ++ed.i;
      continue; // we are in a vcf field
    }

    if (ed.field == 0) /*CHROM field*/
    {
      // check if in header line and store contig
      ed.header_line = buffer_in[ed.b] == '#'; // check if in header line

      if (not ed.header_line)
        ed.next_contig.assign(&buffer_in[ed.b], ed.i - ed.b);
    }
    else if (not ed.header_line)
    {
      if (ed.field == 1) /*POS field*/
      {
        std::from_chars(&buffer_in[ed.b], &buffer_in[ed.i], next_pos);
      }
      else if (ed.field == 4) /*ALT field*/
      {
        int32_t next_n_alt = std::count(&buffer_in[ed.b], &buffer_in[ed.i], ',');
        ed.clear_line(next_pos, next_n_alt);
      }
    }

    if (ed.header_line || ed.field < N_FIELDS_SITE_DATA)
    {
      ++ed.i; // adds '\t' or '\n' and then insert the field to the output buffer
      buffer_out.insert(buffer_out.end(), &buffer_in[ed.b], &buffer_in[ed.i]);
    }
    else
    {
      assert(buffer_in[ed.b] >= '!');
      assert(buffer_in[ed.b] <= '9');

      // check if it is in the current line
      auto insert_it = ed.map_to_unique_fields.insert(
        std::pair<std::string, uint32_t>(std::piecewise_construct,
                                         std::forward_as_tuple(&buffer_in[ed.b], ed.i - ed.b),
                                         std::forward_as_tuple(ed.unique_fields.size())));

      long const field_idx = ed.field - N_FIELDS_SITE_DATA;
      assert(field_idx == static_cast<long>(ed.field2uid.size()));

      if (insert_it.second == true)
      {
        ed.field2uid.push_back(ed.unique_fields.size());
        ed.unique_fields.emplace_back(&buffer_in[ed.b], ed.i - ed.b);

        if (field_idx < static_cast<long>(ed.prev_field2uid.size()) &&
            ed.prev_unique_fields[ed.prev_field2uid[field_idx]] == ed.unique_fields[insert_it.first->second])
        {
          /* Case 0: unique and same as above. */
          buffer_out.push_back('$');

          if (b_in == '\n') /* never skip newline */
            buffer_out.push_back('\n');

          ++ed.i;
        }
        else
        {
          // check if it is in the previous line
          auto prev_find_it = ed.prev_map_to_unique_fields.find(insert_it.first->first);

          if (prev_find_it == ed.prev_map_to_unique_fields.end())
          {
            /* Case 1: Field is unique in the current line and is not in the previous line. */
            ++ed.i; // adds '\t' or '\n'
            buffer_out.insert(buffer_out.end(), &buffer_in[ed.b], &buffer_in[ed.i]);
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

          if (b_in == '\n') /* never skip newline */
            buffer_out.push_back('\n');

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
      ed.field = 0; // reset field index
    else
      ++ed.field;
  } // ends inner loop

  if (ed.field >= 3 && ed.field < N_FIELDS_SITE_DATA)
  {
    // write the data even if the field is not complete
    buffer_out.insert(buffer_out.end(), &buffer_in[ed.b], &buffer_in[ed.i]);

    if (ed.field == 4) /*ALT field*/
      ed.stored_alt = std::count(&buffer_in[ed.b], &buffer_in[ed.i], ',');

    ed.i = 0;
  }
  else
  {
    // copy the remaining data to the beginning of the input buffer
    std::copy(&buffer_in[ed.b], &buffer_in[ed.i], &buffer_in[0]);
    ed.i = ed.i - ed.b;
  }

  ed.b = 0;
  ed.in_size = ed.i;
  resize_input_buffer(buffer_in, ed.i);
}

//! Encode a gzipped file and write to stdout
void encode_file(std::string const & input_fn,
                 bool const is_bgzf_input,
                 std::string const & output_fn,
                 std::string const & output_mode,
                 bool const is_bgzf_output,
                 int const compression_threads);

} // namespace popvcf
