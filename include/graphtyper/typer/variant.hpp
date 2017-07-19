#pragma once

#include <map> // std::map
#include <string> // std::string
#include <vector> // std::vector

#include <graphtyper/graph/genotype.hpp> // gyper::Genotype
#include <graphtyper/graph/haplotype.hpp> // gyper::Haplotype
#include <graphtyper/typer/sample_call.hpp> // gyper::SampleCall

namespace gyper
{

class VariantCandidate;

class Variant
{
public:
  uint32_t abs_pos;
  std::vector<std::vector<char> > seqs;
  std::vector<SampleCall> calls;
  std::map<std::string, std::string> infos;
  std::vector<uint8_t> phase;

  Variant() noexcept;
  Variant(Variant const & var) noexcept;
  Variant(Variant && var) noexcept;
  Variant(Genotype const & gt);
  Variant(std::vector<Genotype> const & gts, std::vector<uint32_t> const & hap_calls);
  Variant(VariantCandidate const & var_candidate) noexcept;

  /******************
   * CLASS MODIFERS *
   ******************/
  void generate_infos();
  void add_call(Haplotype const & hap);
  bool add_base_in_back(bool const add_N = false);
  bool add_base_in_front(bool const add_N = false);
  void normalize(); /** \brief Defined here: http://genome.sph.umich.edu/wiki/Variant_Normalization */
  void trim_sequences(bool const keep_one_match = false);
  void remove_uncalled_alleles();

  /*********************
   * CLASS INFORMATION *
   ********************/
  std::string print() const;  // for debugging
  std::string determine_variant_type() const;
  bool is_normalized() const;
  bool is_snp() const;
  bool is_snp_or_snps() const;
  bool is_transition_snp() const;
  bool is_transversion_snp() const;
  bool is_indel() const;
  bool is_with_matching_first_bases() const;
  uint64_t get_seq_depth() const;
  uint64_t get_seq_depth_of_allele(uint16_t const allele_id) const;
  std::vector<uint64_t> get_seq_depth_of_all_alleles() const;
  uint64_t get_rooted_mapq() const;
  std::vector<uint64_t> get_rooted_mapq_per_allele() const;
  uint64_t get_qual() const; // Gets quality of record


  /*********************
   * OPERATOR OVERLOAD *
   *********************/
  bool operator==(Variant const & b) const;
  bool operator!=(Variant const & b) const;
  bool operator<(Variant const & b) const;
  Variant & operator=(Variant const & var);
};

// Hash function for Variant
struct VariantHash
{
  std::size_t operator()(Variant const & v) const;
};

std::vector<Variant> break_down_variant(Variant && variant, std::size_t const THRESHOLD);
std::vector<Variant> extract_sequences_from_aligned_variant(Variant const && variant, std::size_t const THRESHOLD);
std::vector<Variant> simplify_complex_haplotype(Variant && variant, std::size_t const THRESHOLD);
std::vector<Variant> break_multi_snps(Variant const && var);


} // namespace gyper
