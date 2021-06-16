#pragma once

#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <parallel_hashmap/phmap_fwd_decl.h>

#include <graphtyper/graph/genotype.hpp>
#include <graphtyper/graph/graph.hpp>
#include <graphtyper/typer/genotype_paths.hpp>
#include <graphtyper/typer/read_stats.hpp>
#include <graphtyper/typer/segment.hpp>
#include <graphtyper/utilities/options.hpp>

namespace gyper
{
class Primers;
class HaplotypeCall;

class VcfWriter
{
public:
  explicit VcfWriter(uint32_t variant_distance = 60);
  //  explicit VcfWriter(std::vector<std::string> const & samples, uint32_t variant_distance = 60);

  /*******************
   * CLASS MODIFIERS *
   *******************/
  void set_samples(std::vector<std::string> const & samples);

  void update_haplotype_scores_geno(GenotypePaths & geno, long pn_index, Primers const * primers);

  phmap::flat_hash_map<std::pair<uint16_t, uint16_t>, std::vector<std::pair<uint16_t, uint16_t>>>
  push_to_haplotype_scores(GenotypePaths & geno, long pn_index);

  void update_haplotype_scores_geno(std::pair<GenotypePaths *, GenotypePaths *> & geno_paths,
                                    long const pn_index,
                                    Primers const * primers);

  /*********************
   * CLASS DATA ACCESS *
   *********************/
  std::vector<HaplotypeCall> get_haplotype_calls() const;

private:
  std::unordered_map<uint32_t, uint32_t> id2hap;

public:
  std::vector<std::string> pns;
  std::vector<Haplotype> haplotypes;

/**********************
 * DEBUG ONLY METHODS *
 **********************/
#ifndef NDEBUG
  // void print_variant_group_details() const;
  void print_statistics_headers() const;
  void print_variant_details() const;
  void print_geno_statistics(std::stringstream & read_ss,
                             std::stringstream & path_ss,
                             GenotypePaths const & geno,
                             long pn_index);

  void update_statistics(std::vector<GenotypePaths> & genos, long pn_index);
  void update_statistics(std::vector<std::pair<GenotypePaths, GenotypePaths>> & genos, long pn_index);
  void update_statistics(GenotypePaths const & geno, long pn_index);
#endif // NDEBUG
};

} // namespace gyper
