#include <graphtyper/graph/node.hpp>

namespace gyper
{

template <typename TType>
std::size_t
Node<TType>::label_order() const
{
  return label.order;
}


template <typename TType>
std::size_t
Node<TType>::label_dna_size() const
{
  return label.dna.size();
}


} // namespace gyper
