#pragma once

#include <mutex>
#include <string>
#include <vector>
#include <unordered_map>

#include <graphtyper/typer/variant_candidate.hpp>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>


namespace gyper
{

class ReferenceDepth;
class VariantSupport;

using VarMap = std::unordered_map<VariantCandidate, VariantSupport, VariantCandidateHash>;
using PoolVarMap = std::map<VariantCandidate, std::vector<VariantSupport> >;

class VariantMap
{
  friend class boost::serialization::access;

public:
  void set_samples(std::vector<std::string> const & new_samples);
  void add_variants(std::vector<VariantCandidate> && vars, std::size_t pn_index);
  void create_varmap_for_all(); // Uses the global reference depth
  void filter_varmap_for_all();

  /**
   * \brief Opens input file and creates one variant map with all variant maps in input file
   * \param path Input file path with many variant maps
   */
  void load_many_variant_maps(std::string const & path);
  void write_vcf(std::string const & output_name);

//private:
  void set_pn_count(std::size_t pn_count);

  std::vector<std::string> samples;
  PoolVarMap pool_varmap; /** \brief A varmap for all samples */
  std::vector<VarMap> varmaps; /** \brief List of varmaps, one for each sample */
  std::vector<std::mutex> map_mutexes;

  template <class Archive>
  void inline
  serialize(Archive & ar, unsigned const int)
  {
    ar & samples;
    ar & pool_varmap;
  }
};


extern VariantMap global_varmap;

/**
 * \brief Saves a serialized version of the variant map.
 * \param path Path to save the variant map.
 * \param map Map to save.
 */
void save_variant_map(std::string const & path, VariantMap const & map = global_varmap);

/**
 * \brief Loads a variant map from a serialized file.
 * \param path to the variant map.
 * \return Variant map.
 */
VariantMap load_variant_map(std::string const & path);

} // namespace gyper
