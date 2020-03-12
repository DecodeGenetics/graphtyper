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
  std::vector<uint16_t> calls;
  std::vector<Genotype> gts;
  std::vector<std::vector<ReadStrand> > read_strand;
  std::vector<double> haplotype_impurity;
  long num_samples{0};

  HaplotypeCall() = default;
  explicit HaplotypeCall(Haplotype const & hap);

  void merge_with(HaplotypeCall const & other);
  void make_calls_unique();

private:
  template <class Archive>
  void serialize(Archive & ar, unsigned int version);
};


class HaplotypeCalls
{
  friend class boost::serialization::access;

public:
  HaplotypeCalls() = default;
  explicit HaplotypeCalls(std::vector<HaplotypeCall> const & hap_calls);

  /** \brief Return a list of haplotype calls. */
  inline std::vector<HaplotypeCall>
  get_hap_calls()
  {
    return hap_calls;
  }


private:
  std::vector<HaplotypeCall> hap_calls;

  template <class Archive>
  void serialize(Archive & ar, unsigned int version);
};


void save_calls(HaplotypeCalls & calls, std::string const & filename);
std::vector<HaplotypeCall> load_calls(std::string filename);

}
