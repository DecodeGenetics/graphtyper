#pragma once

#include <cstdint>
#include <vector>

#include <cereal/access.hpp>

#include <graphtyper/constants.hpp>


namespace gyper
{

class Graph;

class KmerLabel
{
private:
  friend class cereal::access;

public:
  uint32_t start_index{0}; /** \brief The index where the variant starts on the reference genome. */
  uint32_t end_index{0};   /** \brief The index where the variant ends on the reference genome. */
  uint32_t variant_id{INVALID_ID};  /** \brief The variant id, i.e. the index of the variant on the graph. */


  KmerLabel() = default;
  KmerLabel(uint32_t const s, uint32_t const e) noexcept;

  // Used when add a KmerLabel to the index.
  KmerLabel(uint32_t const s, uint32_t const e, uint32_t const i) noexcept;

  // Comparison operators
  bool operator==(KmerLabel const & c2) const;
  bool operator!=(KmerLabel const & c2) const;

  std::string to_string() const;


private:
  template <class Archive>
  void
  serialize(Archive & ar, const unsigned int);
};

// Type aliases
using TKmerLabels = std::vector<std::vector<KmerLabel> >;

} // namespace gyper
