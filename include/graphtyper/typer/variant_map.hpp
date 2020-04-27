#pragma once

#include <map> // std::map
#include <string> // std::string
#include <vector> // std::vector
#include <unordered_map>

#include <graphtyper/typer/variant_candidate.hpp>
#include <graphtyper/typer/variant_support.hpp>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>


namespace gyper
{

class ReferenceDepth;
class Vcf;

using VarMap = std::unordered_map<VariantCandidate, VariantSupport, VariantCandidateHash>;
using PoolVarMap = std::map<VariantCandidate, std::vector<VariantSupport> >;

class VariantMap
{
  friend class boost::serialization::access;


public:
  void set_samples(std::vector<std::string> const & new_samples);
  void add_variants(std::vector<VariantCandidate> && vars, long sample_index);
  void create_varmap_for_all(ReferenceDepth const & reference_depth);
  void filter_varmap_for_all();
  void clear();

  /**
   * \brief Opens input file and creates one variant map with all variant maps in input file
   * \param path Input file path with many variant maps
   */
  void load_many_variant_maps(std::string const & path);
  void load_many_variant_maps(std::vector<std::string> const & paths);

#ifndef NDEBUG
  void write_stats(std::string const & prefix = "");
#endif // NDEBUG

  void get_vcf(Vcf & new_variant_vcf, std::string const & output_name);


  //* Instance variables */
  std::vector<std::string> samples;
  PoolVarMap pool_varmap; /** \brief A varmap for all samples */

  // Below is not serialized
  long minimum_variant_support{5};
  double minimum_variant_support_ratio{0.25};
  std::vector<VarMap> varmaps; /** \brief List of varmaps, one for each sample */

  template <class Archive>
  void inline
  serialize(Archive & ar, unsigned const int)
  {
    ar & samples;
    ar & pool_varmap;
    ar & minimum_variant_support;
    ar & minimum_variant_support_ratio;
  }


};


/**
 * \brief Saves a serialized version of the variant map.
 * \param path Path to save the variant map.
 * \param map Map to save.
 */
void save_variant_map(std::string const & path, VariantMap const & map);

/**
 * \brief Loads a variant map from a serialized file.
 * \param path to the variant map.
 * \return Variant map.
 */
VariantMap load_variant_map(std::string const & path);

} // namespace gyper
