#include <fstream>
#include <string>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/log/trivial.hpp>
#include <boost/serialization/vector.hpp>

#include <graphtyper/graph/haplotype_calls.hpp>


namespace gyper
{

HaplotypeCall::HaplotypeCall(std::vector<uint16_t> && _calls,
                             std::vector<Genotype> const & _gts
  )
  : calls(std::move(_calls))
  , gts(_gts)
{}


template<class Archive>
void
HaplotypeCall::serialize(Archive &ar, unsigned const int /*version*/)
{
  ar & calls;
  ar & gts;
}


HaplotypeCalls::HaplotypeCalls(THapCalls const & _hap_calls)
  : hap_calls(_hap_calls)
{}


template <typename Archive>
void
HaplotypeCalls::serialize(Archive & ar, unsigned const int /*version*/)
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
save_calls(HaplotypeCalls & calls, std::string const & filename)
{
  std::ofstream ofs(filename.c_str(), std::ios::binary);

  if (!ofs.is_open())
  {
    BOOST_LOG_TRIVIAL(error) << "[graphtyper::haplotype_calls] Could not save calls to location '"
                             << filename
                             << "'";
    std::exit(1);
  }

  boost::archive::binary_oarchive oa(ofs);
  oa << calls;
}


THapCalls
load_calls(std::string filename)
{
  HaplotypeCalls calls{};
  std::ifstream ifs(filename.c_str(), std::ios::binary);
  assert(ifs.is_open());
  boost::archive::binary_iarchive ia(ifs);
  ia >> calls;
  return calls.get_hap_calls();
}


} // namespace gyper
