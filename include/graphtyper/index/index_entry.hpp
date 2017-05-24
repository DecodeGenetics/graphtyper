#pragma once
#include <cstdint>

#include <seqan/basic.h>

#include <graphtyper/utilities/type_conversions.hpp> // to_uint64()


namespace gyper
{

/**
 * @brief A datastructure of new entries to be added to the index.
 * @details This class is a helper class to the Indexer class. In it we store already calculated information to improve performance of the
 *          index construction algorithm.
 */
class IndexEntry
{
public:
  uint64_t dna;                            /** \brief A string of DNA bases represented as a 64 bit integer. */
  uint32_t start_index;              /** \brief The index where the variant starts on the reference genome. */
  std::vector<uint32_t> variant_id;
  uint32_t total_var_num = 1u;
  unsigned total_var_count = 0u;

  IndexEntry(uint64_t const d, uint32_t const s)
    : dna(d), start_index(s), variant_id(0), total_var_num(1u), total_var_count(0u) {}

  IndexEntry(uint64_t const d, uint32_t const s, uint32_t const i, bool const is_reference, unsigned const var_num)
    : dna(d), start_index(s), variant_id(1, i), total_var_num(var_num)
  {
    if (!is_reference)
    {
      ++total_var_count;
    }
  }


  IndexEntry &
  operator=(IndexEntry rhs)
  {
    std::swap(dna, rhs.dna);
    std::swap(start_index, rhs.start_index);
    std::swap(variant_id, rhs.variant_id);
    std::swap(total_var_num, rhs.total_var_num);
    return *this;
  }


  void inline
  add_to_dna(seqan::Dna const base)
  {
    dna <<= 2; // This is exactly the same as multiplying by four, but this is way cooler.
    dna += seqan::ordValue(base);
  }


};

} // namespace gyper
