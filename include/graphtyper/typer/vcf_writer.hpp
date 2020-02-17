#pragma once

#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <graphtyper/typer/genotype_paths.hpp>
#include <graphtyper/graph/graph.hpp>
#include <graphtyper/graph/genotype.hpp>
#include <graphtyper/typer/read_stats.hpp>
#include <graphtyper/typer/segment.hpp>


namespace gyper
{

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
  void update_haplotype_scores_geno(GenotypePaths & geno, long pn_index);
  void push_to_haplotype_scores(GenotypePaths & geno, long pn_index);

  /*********************
   * CLASS DATA ACCESS *
   *********************/
  std::vector<HaplotypeCall> get_haplotype_calls() const;

private:
  std::unordered_map<uint32_t, std::pair<uint32_t, uint32_t> > id2hap; // first = haplotype, second = local genotype id

public:
  std::vector<std::string> pns;
  std::vector<Haplotype> haplotypes;

/**********************
 * DEBUG ONLY METHODS *
 **********************/
#ifndef NDEBUG
  void print_variant_group_details() const;
  void print_statistics_headers() const;
  void print_variant_details() const;
  void print_geno_statistics(std::stringstream & read_ss,
                             std::stringstream & path_ss,
                             GenotypePaths const & geno,
                             long pn_index
                             );

  void update_statistics(std::vector<GenotypePaths> & genos, long pn_index);
  void update_statistics(std::vector<std::pair<GenotypePaths, GenotypePaths> > & genos,
                         long pn_index
                         );
  void update_statistics(GenotypePaths const & geno, long pn_index);
#endif // NDEBUG

};

} // namespace gyper
