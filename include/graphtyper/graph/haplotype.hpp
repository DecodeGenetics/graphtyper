#pragma once

#include <cstdint> // uint64_t
#include <bitset> // std::bitset
#include <unordered_map> // std::unordered_map
#include <vector> // std::vector

#include <graphtyper/constants.hpp> // MAX_NUMBER_OF_HAPLOTYPES
#include <graphtyper/graph/genotype.hpp> // gyper::Genotype
#include <graphtyper/typer/var_stats.hpp> // gyper::MapQ

namespace gyper
{

struct HapStats
{
  std::vector<uint8_t> hap_coverage;
  std::vector<uint8_t> hap_unique_coverage;
  std::vector<std::vector<std::vector<char> > > predecessor;
  std::vector<std::vector<std::vector<char> > > successor;
  std::vector<std::pair<uint32_t, uint32_t> > hap_b;
  std::vector<std::vector<std::pair<uint32_t, uint32_t> > > pair_info;
};

struct HapSample
{
  std::vector<uint16_t> log_score;
  std::vector<std::vector<uint16_t> > gt_coverage;
  std::unique_ptr<HapStats> stats;
  // uint8_t ambiguous_depth = 0; // TODO
  uint16_t max_log_score = 0u;
};


class Haplotype
{
public:
  std::vector<Genotype> gts; /** \brief A list of genotypes this haplotype has. */
  std::vector<HapSample> hap_samples;
  std::vector<std::pair<uint16_t, uint16_t> > calls; /** \brief All calls of each sample. This is used for segment calling only. */
#ifndef NDEBUG
  std::vector<std::bitset<MAX_NUMBER_OF_HAPLOTYPES> > unique_gts; // One per gt
#endif // NDEBUG
  std::vector<VarStats> var_stats;

  Haplotype() noexcept;

  /*******************
   * CLASS MODIFIERS *
   *******************/
  void add_genotype(Genotype && gt);
  void check_for_duplicate_haplotypes();
  void clear_and_resize_samples(std::size_t const new_size);
  void clear();
  void add_coverage(uint32_t const local_genotype_id, uint16_t const c);
  void add_explanation(uint32_t const local_genotype_id, std::bitset<MAX_NUMBER_OF_HAPLOTYPES> const & e);
  void update_max_log_score();

  /*********************
   * CLASS INFORMATION *
   *********************/
  std::vector<uint16_t> get_haplotype_calls() const;
  uint32_t get_genotype_num() const;
  bool has_too_many_genotypes() const;
  std::vector<uint32_t> get_genotype_ids() const;
  uint32_t best_score_of_a_path(std::size_t const pn_index,
                                std::bitset<MAX_NUMBER_OF_HAPLOTYPES> const & e1,
                                std::bitset<MAX_NUMBER_OF_HAPLOTYPES> const & e2
                                /*std::vector<uint16_t> const & max_scores_per_y*/
                                ) const;

  /** Update likelihood and stats */
  void explain_to_score(std::size_t const pn_index, bool const has_low_quality_snp, bool const non_unique_paths, uint8_t const mapq, bool const fully_aligned, std::size_t const mismatches);
  void clipped_reads_to_stats(bool const fully_aligned);
  void graph_complexity_to_stats();
  void mapq_to_stats(uint8_t const mapq);
  void realignment_to_stats(bool const is_unaligned_read, uint32_t const original_pos, uint32_t const new_pos);
  void strand_to_stats(bool const forward_strand, bool const is_first_in_pair);
  void coverage_to_gts(std::size_t const pn_index);
  std::bitset<MAX_NUMBER_OF_HAPLOTYPES> explain_to_path_explain();

private:
  std::vector<uint16_t> coverage; // per gt. Contains the called allele
  std::vector<std::bitset<MAX_NUMBER_OF_HAPLOTYPES> > explains; // per gt

  uint32_t resize_score_vectors(std::size_t const pn_index);
  std::bitset<MAX_NUMBER_OF_HAPLOTYPES> find_which_haplotypes_explain_the_read(uint32_t const cnum) const;
  std::vector<uint16_t> find_with_how_many_errors_haplotypes_explain_the_read(uint32_t const cnum) const;
};

} // namespace gyper
