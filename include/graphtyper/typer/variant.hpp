#pragma once

#include <map> // std::map
#include <string> // std::string
#include <vector> // std::vector

#include <boost/serialization/access.hpp>

#include <graphtyper/graph/genotype.hpp> // gyper::Genotype
#include <graphtyper/graph/haplotype.hpp> // gyper::Haplotype
#include <graphtyper/typer/sample_call.hpp> // gyper::SampleCall
#include <graphtyper/typer/var_stats.hpp>
#include <graphtyper/utilities/options.hpp>

namespace gyper
{

class VariantCandidate;

class Variant
{
  friend class boost::serialization::access;

public:
  uint32_t abs_pos;
  std::vector<std::vector<char> > seqs;
  std::vector<SampleCall> calls;
  VarStats stats;
  std::map<std::string, std::string> infos;
  std::string suffix_id;
  bool is_info_generated{false};
  char type{'.'};

  Variant() noexcept;

  // The rule of 5
  Variant(Variant const & var) noexcept;
  Variant(Variant && var) noexcept;
  Variant & operator=(Variant const & o) noexcept;
  Variant & operator=(Variant && o) noexcept;
  ~Variant() = default;

  Variant(Genotype const & gt);
  Variant(std::vector<Genotype> const & gts, std::vector<uint16_t> const & hap_calls);
  Variant(VariantCandidate const & var_candidate) noexcept;

  /******************
   * CLASS MODIFERS *
   ******************/
  void update_camou_phred(long const ploidy);
  void generate_infos();
  bool add_base_in_back(bool const add_N = false);
  bool add_base_in_front(bool const add_N = false);
  void expanded_normalized();
  void normalize(); /** \brief Defined here: http://genome.sph.umich.edu/wiki/Variant_Normalization */
  void remove_common_prefix(bool const keep_one_match = false);
  void trim_sequences(bool const keep_one_match = false);

  /*********************
   * CLASS INFORMATION *
   ********************/
  std::string print() const;  // for debugging
  std::string determine_variant_type() const;
  bool is_normalized() const;
  bool is_snp_or_snps() const;
  bool is_with_matching_first_bases() const;
  bool is_sv() const;
  uint64_t get_seq_depth() const;
  uint64_t get_seq_depth_of_allele(uint16_t const allele_id) const;
  std::vector<uint64_t> get_seq_depth_of_all_alleles() const;
  uint64_t get_qual() const; // Gets quality of a record
  double get_qual_by_depth() const; // Gets total quality by depth (QD) of a record
  std::vector<double> get_qual_by_depth_per_alt_allele() const;

  /*********************
   * OPERATOR OVERLOAD *
   *********************/
  bool operator==(Variant const & b) const;
  bool operator!=(Variant const & b) const;
  bool operator<(Variant const & b) const;

private:
  template <class Archive>
  void serialize(Archive & ar, unsigned int version);
};

// Hash function for Variant
struct VariantHash
{
  std::size_t operator()(Variant const & v) const;
};

std::vector<Variant> break_down_variant(Variant && variant,
                                        long const reach,
                                        bool const is_no_variant_overlapping,
                                        bool const is_all_biallelic);

std::vector<Variant> break_down_skyr(Variant && var, long const reach);
std::vector<Variant> extract_sequences_from_aligned_variant(Variant const && variant, std::size_t const THRESHOLD);
std::vector<Variant> simplify_complex_haplotype(Variant && variant, std::size_t const THRESHOLD);
std::vector<Variant> break_multi_snps(Variant const && var);
void find_variant_sequences(gyper::Variant & new_var, gyper::Variant const & old_var);

} // namespace gyper
