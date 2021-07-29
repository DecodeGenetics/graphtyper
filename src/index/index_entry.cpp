#include <cstdint>
#include <set>

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/snp_event.hpp>
#include <graphtyper/index/index_entry.hpp>

namespace gyper
{
IndexEntry::IndexEntry(uint32_t const s) : start_index(s)
{
}

IndexEntry::IndexEntry(uint32_t const s, uint32_t const i, bool const is_reference, unsigned const var_num) :
  start_index(s), total_var_num(var_num), total_var_count(static_cast<uint32_t>(!is_reference))
{
  variant_id.insert(i);
}

void IndexEntry::add_to_dna(char const base)
{
  dna <<= 2; // This is exactly the same as multiplying by four, but this is way cooler.

  if (valid > 0)
  {
    --valid;
  }
  else
  {
    switch (base)
    {
    case 'A':
      break;

    case 'C':
      dna += 1;
      break;

    case 'G':
      dna += 2;
      break;

    case 'T':
      dna += 3;
      break;

    default:
    {
      valid = static_cast<uint8_t>(gyper::K);
      break;
    }
    }
  }
}

} // namespace gyper
