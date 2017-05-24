#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>

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


std::vector<TNodeIndex>
RefNode::get_vars() const
{
  return out_var_ids;
}


void
RefNode::change_label_order(uint32_t change)
{
  // Make sure we do not overflow
  assert(change + label.order >= change);
  assert(change + label.order >= label.order);

  label.order += change;
}


std::size_t
RefNode::out_degree() const
{
  return out_var_ids.size();
}


/***********
 * PRIVATE *
 ***********/

RefNode::RefNode()
  : label(), out_var_ids(0) {}


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
template void RefNode::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive &, const unsigned int);
template void RefNode::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive &, const unsigned int);

} // namespace gyper
