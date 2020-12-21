#pragma once

#include <cstdint>
#include <set>
#include <unordered_set>
#include <vector>

#include <graphtyper/graph/snp_event.hpp>


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
  uint64_t dna{0u}; /** \brief A string of DNA bases represented as a 64 bit integer. */
  uint32_t start_index{0u};              /** \brief The index where the variant starts on the reference genome. */
  std::set<uint32_t> variant_id;
  uint32_t total_var_num{1u};
  uint32_t total_var_count{0u};
  uint8_t valid{0u};
  std::unordered_set<long> events;
  std::unordered_set<long> anti_events;

  IndexEntry(uint32_t const s);
  IndexEntry(uint32_t const s, uint32_t const i, bool const is_reference, unsigned const var_num);

  void add_to_dna(char const base);
};

} // namespace gyper
