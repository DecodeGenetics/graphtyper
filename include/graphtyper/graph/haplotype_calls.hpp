#pragma once
#include <iostream>
#include <vector>

#include <boost/serialization/access.hpp>

#include <graphtyper/graph/haplotype.hpp>
#include <graphtyper/typer/vcf_writer.hpp>

namespace gyper
{

class HaplotypeCall
{
  friend class boost::serialization::access;

public:
  std::vector<uint16_t> calls{};
  std::vector<Genotype> gts{};

  HaplotypeCall() = default;
  explicit HaplotypeCall(std::vector<uint16_t> && _calls,
                std::vector<Genotype> const & _gts
    );

private:
  template <class Archive>
  void serialize(Archive & ar, unsigned int version);
};


using THapCalls = std::vector<HaplotypeCall>;

class HaplotypeCalls
{
  friend class boost::serialization::access;

public:
  HaplotypeCalls() = default;
  explicit HaplotypeCalls(THapCalls const & hap_calls);

  /** \brief Return a list of haplotype calls. */
  THapCalls inline get_hap_calls(){return hap_calls;}

private:
  THapCalls hap_calls;

  template <class Archive>
  void serialize(Archive & ar, unsigned int version);
};


void save_calls(HaplotypeCalls & calls, std::string const & filename);
THapCalls load_calls(std::string filename);

}
