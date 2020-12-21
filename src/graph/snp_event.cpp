#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <graphtyper/graph/snp_event.hpp>

#include <cstdint>
#include <string>
#include <sstream>


namespace gyper
{

SnpEvent::SnpEvent(uint32_t _pos, char _base)
  : pos(_pos)
  , base(_base)
{}

std::string
SnpEvent::to_string() const
{
  std::ostringstream ss;
  ss << pos << " " << base;
  return ss.str();
}


bool
SnpEvent::operator==(SnpEvent const & b) const
{
  return (pos == b.pos) && (base == b.base);
}


bool
SnpEvent::operator!=(SnpEvent const & b) const
{
  return !(*this == b);
}


bool
SnpEvent::operator<(SnpEvent const & b) const
{
  return pos < b.pos || (pos == b.pos && base < b.base);
}


template <typename Archive>
void
SnpEvent::serialize(Archive & ar, const unsigned int)
{
  ar & pos;
  ar & base;
}


/***************************
 * EXPLICIT INSTANTIATIONS *
 ***************************/
template void SnpEvent::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive &,
                                                                   const unsigned int);
template void SnpEvent::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive &,
                                                                   const unsigned int);

} // namespace gyper
