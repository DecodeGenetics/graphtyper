#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>

#include <graphtyper/graph/haplotype_calls.hpp>

namespace gyper
{

HaplotypeCalls::HaplotypeCalls()
  : hap_calls(0)
{}

HaplotypeCalls::HaplotypeCalls(THapCalls const & _hap_calls)
  : hap_calls(_hap_calls)
{}

template <typename Archive>
void
HaplotypeCalls::serialize(Archive & ar, const unsigned int)
{
  ar & hap_calls;
}


/***************************
 * EXPLICIT INSTANTIATIONS *
 ***************************/

template void HaplotypeCalls::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive &, const unsigned int);
template void HaplotypeCalls::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive &, const unsigned int);

/*********************************
 * FUNCTIONS TO MAKE LIFE EASIER *
 *********************************/

void
save_calls(HaplotypeCalls & calls, std::string filename)
{
  std::ofstream ofs(filename.c_str(), std::ios::binary);

  if (!ofs.is_open())
  {
    std::cerr << "Could not save calls to location '" << filename << "'" << std::endl;
    std::exit(1);
  }

  boost::archive::binary_oarchive oa(ofs);
  oa << calls;
}


void
load_calls(HaplotypeCalls & calls, std::string filename)
{
  std::ifstream ifs(filename.c_str(), std::ios::binary);
  assert(ifs.is_open());
  boost::archive::binary_iarchive ia(ifs);
  ia >> calls;
}


THapCalls
load_calls(std::string filename)
{
  HaplotypeCalls calls;
  std::ifstream ifs(filename.c_str(), std::ios::binary);
  assert(ifs.is_open());
  boost::archive::binary_iarchive ia(ifs);
  ia >> calls;
  return calls.hap_calls;
}


} // namespace gyper
