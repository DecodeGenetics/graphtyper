#include <cstdint>
#include <string>
#include <sstream>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include <graphtyper/graph/sv.hpp>


std::string
SV::get_type() const
{
  switch(type)
  {
  case DEL:
    return "DEL";

  case DEL_ALU:
    return "DEL:ME:ALU";

  case DUP:
    return "DUP";

  case INS:
    return "INS";

  case INS_ALU:
    return "INS:ME:ALU";

  case INV:
    return "INV";

  default:
    return "SV";

  } // switch
}


std::string
SV::get_allele() const
{
  std::ostringstream ss;

  ss << '<'
     << get_type()
     << ":SVSIZE="
     << size
     << '>';

  return ss.str();
}


template <typename Archive>
void
SV::serialize(Archive & ar, const unsigned int /*version*/)
{
  ar & type;
  ar & length;
  ar & size;
  ar & end;
  ar & n_clusters;
  ar & or_end;
  ar & or_start;
  ar & seq;
}


/***************************
 * EXPLICIT INSTANTIATIONS *
 ***************************/

template void SV::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive &, const unsigned int);
template void SV::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive &, const unsigned int);
