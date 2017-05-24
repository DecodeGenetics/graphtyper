#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <graphtyper/graph/node.hpp>


namespace gyper
{

/**********
 * PUBLIC *
 **********/

VarNode::VarNode(Label && l, TNodeIndex && ori) noexcept
  : label(std::move(l))
  , out_ref_id(std::move(ori))
{}

VarNode::VarNode(VarNode const & vn) noexcept
  : label(vn.label)
  , out_ref_id(vn.out_ref_id)
{}

Label const &
VarNode::get_label() const
{
  return label;
}


TNodeIndex
VarNode::get_out_ref_index() const
{
  return out_ref_id;
}


void
VarNode::change_label_order(uint32_t change)
{
  // Make sure we do not overflow
  assert(change + label.order >= change);
  assert(change + label.order >= label.order);

  label.order += change;
}


std::size_t
VarNode::out_degree() const
{
  return 1;
}


/***********
 * PRIVATE *
 ***********/
VarNode::VarNode()
  : label(), out_ref_id(0) {}

template <typename Archive>
void
VarNode::serialize(Archive & ar, const unsigned int)
{
  ar & label;
  ar & out_ref_id;
}


/***************************
 * EXPLICIT INSTANTIATIONS *
 ***************************/
template void VarNode::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive &, const unsigned int);
template void VarNode::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive &, const unsigned int);

} // namespace gyper
