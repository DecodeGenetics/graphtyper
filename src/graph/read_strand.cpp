#include <algorithm>
#include <sstream>

#include <cereal/archives/binary.hpp>

#include <graphtyper/graph/read_strand.hpp>


namespace gyper
{

long
ReadStrand::get_weight() const
{
  return r1_forward + r1_reverse + r2_forward + r2_reverse;
}


long
ReadStrand::get_r1_count() const
{
  return r1_forward + r1_reverse;
}


long
ReadStrand::get_reverse_count() const
{
  return r1_reverse + r2_reverse;
}


long
ReadStrand::get_max_bias() const
{
  return std::max({r1_forward, r1_reverse, r2_forward, r2_reverse});
}


void
ReadStrand::merge_with(ReadStrand const & o)
{
  r1_forward += o.r1_forward;
  r1_reverse += o.r1_reverse;
  r2_forward += o.r2_forward;
  r2_reverse += o.r2_reverse;
}


std::string
ReadStrand::str() const
{
  std::ostringstream ss;
  ss << "F1,R1,F2,R2 counts are " << r1_forward << "," << r1_reverse << "," << r2_forward << "," << r2_reverse << '\n';
  return ss.str();
}


template <typename Archive>
void
ReadStrand::serialize(Archive & ar, unsigned const int /*version*/)
{
  ar & r1_forward;
  ar & r1_reverse;
  ar & r2_forward;
  ar & r2_reverse;
}


/***************************
 * EXPLICIT INSTANTIATIONS *
 ***************************/

template void ReadStrand::serialize<cereal::BinaryInputArchive>(cereal::BinaryInputArchive &,
                                                                     const unsigned int);
template void ReadStrand::serialize<cereal::BinaryOutputArchive>(cereal::BinaryOutputArchive &,
                                                                     const unsigned int);

} // namespace gyper
