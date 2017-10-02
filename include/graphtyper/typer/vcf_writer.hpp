#pragma once

#include <array>
#include <bitset>
#include <iostream>
#include <map>
#include <vector>
#include <sstream>
#include <string>
#include <mutex>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/vcf_io.h>

#include <graphtyper/typer/genotype_paths.hpp>
#include <graphtyper/graph/graph.hpp>
#include <graphtyper/graph/genotype.hpp>
#include <graphtyper/typer/read_stats.hpp>
#include <graphtyper/typer/segment.hpp>


namespace gyper
{

class VcfWriter
{
  /*******************
   * TYPE DEFINTIONS *
   *******************/
  using ExplainMap = std::map<uint32_t, std::vector<std::bitset<MAX_NUMBER_OF_HAPLOTYPES> > >;

public:
  VcfWriter(std::vector<std::string> const & samples, uint32_t variant_distance = 60);

  std::vector<std::vector<double> > hap_scores;

  /*******************
   * CLASS MODIFIERS *
   *******************/
  // read_pair = 0 is unpaired, 1 is first in pair, 2 is second in pair
  //void update_statistics(GenotypePaths & geno, std::size_t const pn_index, unsigned const read_pair);

  void print_haplotype_details() const;
  void print_statistics_headers() const;
  void print_variant_details() const;
  void print_geno_statistics(std::stringstream & read_ss,
                             std::stringstream & path_ss,
                             GenotypePaths const & geno,
                             std::size_t const pn_index
    );

  void update_statistics(std::vector<GenotypePaths> & genos, std::size_t const pn_index);
  void update_statistics(std::vector<std::pair<GenotypePaths, GenotypePaths> > & genos,
                         std::size_t const pn_index
    );

  void update_haplotype_scores_from_path(GenotypePaths & geno, std::size_t const pn_index, unsigned const read_pair);
  void update_haplotype_scores_from_paths(std::vector<GenotypePaths> & genos, std::size_t const pn_index);
  void update_haplotype_scores_from_paths(std::vector<std::pair<GenotypePaths, GenotypePaths> > & genos, std::size_t const pn_index);

  void find_path_explanation(GenotypePaths const & gt_path,
                             std::vector<std::pair<uint32_t, std::bitset<MAX_NUMBER_OF_HAPLOTYPES> > > & ids_and_path_explain
                             );

  std::vector<uint32_t>
  explain_map_to_haplotype_scores(std::size_t const pn_index,
                                  ExplainMap const & explain_map
                                  );

  uint32_t
  explain_map_specific_indexes_to_haplotype_scores(
    std::size_t const pn_index,
    std::pair<uint32_t, uint32_t> const index,
    ExplainMap const & explain_map
  ) const;


  /*********************
   * CLASS DATA ACCESS *
   *********************/
  std::vector<std::pair<std::vector<uint16_t>, std::vector<Genotype> > > get_haplotype_calls() const;
  std::vector<Genotype> get_gts() const;

private:
  std::mutex mutable haplotype_mutex;
  std::mutex mutable io_mutex;
  std::vector<std::string> pns;
  std::string pn;
  std::unordered_map<uint32_t, std::pair<uint32_t, uint32_t> > id2hap; // first = haplotype, second = local genotype id
  std::unique_ptr<ReadStats> read_stats = nullptr;

  seqan::VcfRecord create_multi_record(Genotype const & gt, bool const WRITE_PL);

public:
  std::vector<Haplotype> haplotypes;
  void generate_statistics();

};

} // namespace gyper
