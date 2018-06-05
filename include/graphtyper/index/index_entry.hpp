#pragma once
#include <cstdint>

#include <graphtyper/utilities/type_conversions.hpp> // to_uint64()


namespace gyper
{

/**
 * @brief A datastructure of new entries to be added to the index.
 * @details This class is a helper class to the Indexer class. In it we store already calculated
 *          information to improve performance of the index construction algorithm.
 */
class IndexEntry
{
public:
  uint64_t dna = 0u;                      /** \brief A string of DNA bases represented as a 64 bit integer. */
  uint32_t start_index = 0u;              /** \brief The index where the variant starts on the reference genome. */
  std::vector<uint32_t> variant_id;
  uint32_t total_var_num = 1u;
  uint32_t total_var_count = 0u;
  uint8_t valid = 0u;

  IndexEntry(uint32_t const s)
    : start_index(s)
  {}

  IndexEntry(uint32_t const s, uint32_t const i, bool const is_reference, unsigned const var_num)
    : start_index(s), variant_id(1, i), total_var_num(var_num), total_var_count(static_cast<uint32_t>(!is_reference))
  {}

  void inline
  add_to_dna(char const base)
  {
    dna <<= 2; // This is exactly the same as multiplying by four, but this is way cooler.

    if (valid > 0)
    {
      --valid;
    }
    else
    {
      switch(base)
      {
      case 'A': break;
      case 'C': dna += 1; break;
      case 'G': dna += 2; break;
      case 'T': dna += 3; break;

      default:
      {
        //std::cerr << "[graphtyper::index_entry] Non-valid base: '" << base << "' \n";
        valid = static_cast<uint8_t>(gyper::K);
        break;
      }

      }
    }
  }


};

} // namespace gyper
