#include <cassert>
#include <iostream>

#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>

#include <graphtyper/graph/node.hpp>


namespace gyper
{

RefNode::RefNode(Label && l, std::vector<TNodeIndex> && c) noexcept
  : label(std::move(l))
  , out_var_ids(std::move(c))
{}


RefNode::RefNode(RefNode const & rn) noexcept
  : label(rn.label)
  , out_var_ids(rn.out_var_ids)
{}


RefNode::RefNode(RefNode && o) noexcept
  : label(std::move(o.label))
  , out_var_ids(std::move(o.out_var_ids))
{}


Label const &
RefNode::get_label() const
{
  return label;
}


TNodeIndex
RefNode::get_var_index(unsigned const & index) const
{
  assert(index < out_var_ids.size());
  return out_var_ids[index];
}


std::vector<TNodeIndex> const &
RefNode::get_vars() const
{
  return out_var_ids;
}


std::size_t
RefNode::out_degree() const
{
  return out_var_ids.size();
}


/***********
 * PRIVATE *
 ***********/

template <typename Archive>
void
RefNode::serialize(Archive & ar, const unsigned int)
{
  ar & label;
  ar & out_var_ids;
}


/***************************
 * EXPLICIT INSTANTIATIONS *
 ***************************/
template void RefNode::serialize<cereal::BinaryInputArchive>(cereal::BinaryInputArchive &,
                                                                  const unsigned int);
template void RefNode::serialize<cereal::BinaryOutputArchive>(cereal::BinaryOutputArchive &,
                                                                  const unsigned int);

} // namespace gyper
