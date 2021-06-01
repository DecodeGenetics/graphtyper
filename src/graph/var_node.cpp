#include <iostream>
#include <unordered_set>

#include <cereal/archives/binary.hpp>
#include <cereal/types/unordered_set.hpp>

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/node.hpp>

namespace gyper
{
/**********
 * PUBLIC *
 **********/

// VarNode::VarNode(Label && l, TNodeIndex && ori, std::unordered_set<long> && events) noexcept
//  : label(std::forward<Label>(l))
//  , out_ref_id(ori)
//  , events(std::forward<std::unordered_set<long> >(events))
//{}

VarNode::VarNode(Label && l,
                 TNodeIndex && ori,
                 std::unordered_set<long> && events,
                 std::unordered_set<long> && anti_events) noexcept :
  label(std::forward<Label>(l)),
  out_ref_id(ori),
  events(std::forward<std::unordered_set<long>>(events)),
  anti_events(std::forward<std::unordered_set<long>>(anti_events))
{
}

VarNode::VarNode(VarNode const & vn) noexcept :
  label(vn.label), out_ref_id(vn.out_ref_id), events(vn.events), anti_events(vn.anti_events)
{
}

VarNode::VarNode(VarNode && o) noexcept :
  label(std::move(o.label)),
  out_ref_id(std::move(o.out_ref_id)),
  events(std::move(o.events)),
  anti_events(std::move(o.anti_events))
{
}

Label const & VarNode::get_label() const
{
  return label;
}

TNodeIndex VarNode::get_out_ref_index() const
{
  return out_ref_id;
}

std::size_t VarNode::out_degree() const
{
  return 1;
}

/***********
 * PRIVATE *
 ***********/
template <typename Archive>
void VarNode::serialize(Archive & ar, unsigned int const)
{
  ar & label;
  ar & out_ref_id;
  ar & events;
  ar & anti_events;
}

/***************************
 * EXPLICIT INSTANTIATIONS *
 ***************************/
template void VarNode::serialize<cereal::BinaryInputArchive>(cereal::BinaryInputArchive &, const unsigned int);
template void VarNode::serialize<cereal::BinaryOutputArchive>(cereal::BinaryOutputArchive &, const unsigned int);

} // namespace gyper
