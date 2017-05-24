#pragma once

#include <mutex>
#include <string>
#include <vector>
#include <unordered_map>

#include <graphtyper/typer/variant_candidate.hpp>


namespace gyper
{

class ReferenceDepth;
class VariantSupport;

using VarMap = std::unordered_map<VariantCandidate, std::pair<uint32_t, uint32_t>, VariantCandidateHash>;
using PoolVarMap = std::map<VariantCandidate, std::vector<VariantSupport> >;

class VariantMap
{
public:
  void set_pn_count(std::size_t const pn_count);
  void add_variants(std::vector<VariantCandidate> && vars, std::size_t const pn_index);
  void create_varmap_for_all(); // Uses the global reference depth
  void filter_varmap_for_all();
  void write_vcf(std::string const & output_name);

private:
  PoolVarMap pool_varmap; /** \brief A varmap for all samples */
  std::vector<VarMap> varmaps; /** \brief List of varmaps, one for each sample */
  std::vector<std::mutex> map_mutexes;
};


extern VariantMap global_varmap;

} // namespace gyper
