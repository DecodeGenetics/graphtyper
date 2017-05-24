#pragma once

#include <array>
#include <bitset>
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <mutex>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/vcf_io.h>

#include <graphtyper/typer/genotype_paths.hpp>
#include <graphtyper/graph/graph.hpp>
#include <graphtyper/graph/genotype.hpp>
#include <graphtyper/typer/segment.hpp>


namespace gyper
{

class VcfWriter
{
public:
  VcfWriter(std::vector<std::string> const & samples, uint32_t variant_distance = 60);

  std::vector<std::vector<double> > hap_scores;

  /*******************
   * CLASS MODIFIERS *
   *******************/
  // read_pair = 0 is unpaired, 1 is first in pair, 2 is second in pair
  void update_haplotype_scores_from_path(GenotypePaths & geno, std::size_t const pn_index, unsigned const read_pair);
  void update_haplotype_scores_from_paths(std::vector<GenotypePaths> & genos, std::size_t const pn_index);
  void update_haplotype_scores_from_paths(std::vector<std::pair<GenotypePaths, GenotypePaths> > & genos, std::size_t const pn_index);

  void find_path_explanation(GenotypePaths const & gt_path,
                             std::vector<std::pair<uint32_t, std::bitset<MAX_NUMBER_OF_HAPLOTYPES> > > & ids_and_path_explain
                             );

  std::vector<uint32_t> explain_map_to_haplotype_scores(std::size_t const pn_index,
                                                        std::map<uint32_t, std::vector<std::bitset<MAX_NUMBER_OF_HAPLOTYPES> > > const & explain_map
                                                        );
  uint32_t
  explain_map_specific_indexes_to_haplotype_scores(
    std::size_t const pn_index,
    std::pair<uint32_t, uint32_t> const index,
    std::map<uint32_t, std::vector<std::bitset<MAX_NUMBER_OF_HAPLOTYPES> > > const & explain_map
  ) const;


  /*********************
   * CLASS DATA ACCESS *
   *********************/
  std::vector<std::pair<std::vector<uint16_t>, std::vector<Genotype> > > get_haplotype_calls() const;
  std::vector<Genotype> get_gts() const;

private:
  std::mutex mutable haplotype_mutex;
  std::string pn;
  std::unordered_map<uint32_t, std::pair<uint32_t, uint32_t> > id2hap; // first = haplotype, second = local genotype id

  seqan::VcfRecord create_multi_record(Genotype const & gt, bool const WRITE_PL);

public:
  std::vector<Haplotype> haplotypes;
  void generate_statistics(std::vector<std::string> const & pns);

};

std::string inline
get_contig(std::string const & chromosome)
{
  std::unordered_map<std::string, std::string> contig_map;

  contig_map["chr1"] = "<ID=chr1,length=248956422>";
  contig_map["chr10"] = "<ID=chr10,length=133797422>";
  contig_map["chr11"] = "<ID=chr11,length=135086622>";
  contig_map["chr12"] = "<ID=chr12,length=133275309>";
  contig_map["chr13"] = "<ID=chr13,length=114364328>";
  contig_map["chr14"] = "<ID=chr14,length=107043718>";
  contig_map["chr15"] = "<ID=chr15,length=101991189>";
  contig_map["chr16"] = "<ID=chr16,length=90338345>";
  contig_map["chr17"] = "<ID=chr17,length=83257441>";
  contig_map["chr18"] = "<ID=chr18,length=80373285>";
  contig_map["chr19"] = "<ID=chr19,length=58617616>";
  contig_map["chr2"] = "<ID=chr2,length=242193529>";
  contig_map["chr20"] = "<ID=chr20,length=64444167>";
  contig_map["chr21"] = "<ID=chr21,length=46709983>";
  contig_map["chr22"] = "<ID=chr22,length=50818468>";
  contig_map["chr3"] = "<ID=chr3,length=198295559>";
  contig_map["chr4"] = "<ID=chr4,length=190214555>";
  contig_map["chr5"] = "<ID=chr5,length=181538259>";
  contig_map["chr6"] = "<ID=chr6,length=170805979>";
  contig_map["chr7"] = "<ID=chr7,length=159345973>";
  contig_map["chr8"] = "<ID=chr8,length=145138636>";
  contig_map["chr9"] = "<ID=chr9,length=138394717>";

  if (contig_map.count(chromosome) == 1)
    return contig_map[chromosome];
  else
    return std::string();
}


} // namespace gyper
