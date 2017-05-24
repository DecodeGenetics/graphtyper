#pragma once
#include <iostream>
#include <vector>

#include <boost/serialization/access.hpp>

#include <graphtyper/graph/haplotype.hpp>
#include <graphtyper/typer/vcf_writer.hpp>

namespace gyper
{

using THapCalls = std::vector<std::pair<std::vector<uint16_t>, std::vector<Genotype> > >;

class HaplotypeCalls
{
  friend class boost::serialization::access;

public:
  HaplotypeCalls();
  HaplotypeCalls(THapCalls const & hap_calls);
  THapCalls hap_calls;

private:
  template <class Archive>
  void serialize(Archive & ar, const unsigned int version);
};


void save_calls(HaplotypeCalls & calls, std::string filename);
void load_calls(HaplotypeCalls & calls, std::string filename);
THapCalls load_calls(std::string filename);

}
