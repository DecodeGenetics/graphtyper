#include <algorithm> // std::swap
#include <limits>
#include <string> // std::string
#include <sstream> // std::stringstream
#include <vector> // std::vector

#include <boost/log/trivial.hpp>
#include <boost/functional/hash.hpp> // boost::hash_range

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

#include <paw/align.hpp>

#include <graphtyper/graph/absolute_position.hpp> // gyper::AbsolutePosition
#include <graphtyper/graph/graph.hpp> // gyper::Graph
#include <graphtyper/graph/sv.hpp> // gyper::SVTYPE
#include <graphtyper/typer/variant.hpp> // gyper::Variant
#include <graphtyper/typer/variant_candidate.hpp> // gyper::VariantCandidate
#include <graphtyper/utilities/graph_help_functions.hpp> // gyper::to_pair, gyper::to_index
#include <graphtyper/utilities/options.hpp> // gyper::Options
#include <graphtyper/utilities/sequence_operations.hpp> // gyper::remove_common_prefix, gyper::remove_common_suffix


namespace
{

void
update_strand_bias(std::size_t const num_seqs,
                   std::size_t new_num_seqs,
                   std::vector<uint16_t> const & old_phred_to_new_phred,
                   gyper::Variant const & var,
                   gyper::Variant & new_var)
{
  using gyper::get_strand_bias;
  using gyper::join_strand_bias;

  std::vector<uint64_t> seq_depths = var.get_seq_depth_of_all_alleles();

  std::vector<uint32_t> sbf = get_strand_bias(var.infos, "SBF");
  std::vector<uint32_t> sbr = get_strand_bias(var.infos, "SBR");

  std::vector<uint32_t> sbf1 = get_strand_bias(var.infos, "SBF1");
  std::vector<uint32_t> sbf2 = get_strand_bias(var.infos, "SBF2");
  std::vector<uint32_t> sbr1 = get_strand_bias(var.infos, "SBR1");
  std::vector<uint32_t> sbr2 = get_strand_bias(var.infos, "SBR2");

  assert(sbf.size() > 0);
  assert(sbf.size() == num_seqs);
  assert(sbr.size() == num_seqs);
  assert(sbf1.size() == num_seqs);
  assert(sbf2.size() == num_seqs);
  assert(sbr1.size() == num_seqs);
  assert(sbr2.size() == num_seqs);

//  std::vector<uint32_t> ra_count = get_strand_bias(var.infos, "RACount");
//  std::vector<uint32_t> ra_dist = get_strand_bias(var.infos, "RADist");

  if (sbf.size() > 0 && sbr.size() > 0)
  {
    std::vector<uint32_t> new_sbf(new_num_seqs, 0u);
    std::vector<uint32_t> new_sbr(new_num_seqs, 0u);

    std::vector<uint32_t> new_sbf1(new_num_seqs, 0u);
    std::vector<uint32_t> new_sbf2(new_num_seqs, 0u);
    std::vector<uint32_t> new_sbr1(new_num_seqs, 0u);
    std::vector<uint32_t> new_sbr2(new_num_seqs, 0u);

//    std::vector<uint32_t> new_ra_count(new_num_seqs, 0u);
//    std::vector<uint32_t> new_ra_dist(new_num_seqs, 0u);

    for (uint32_t y = 0; y < num_seqs; ++y)
    {
      uint32_t new_y = old_phred_to_new_phred[y];

      assert(new_y < new_num_seqs);
      assert(y < seq_depths.size());
      assert(y < sbf.size());
      assert(y < sbr.size());
      assert(y < sbf1.size());
      assert(y < sbf2.size());
      assert(y < sbr1.size());
      assert(y < sbr2.size());

      new_sbf[new_y] += sbf[y];
      new_sbr[new_y] += sbr[y];

      new_sbf1[new_y] += sbf1[y];
      new_sbf2[new_y] += sbf2[y];
      new_sbr1[new_y] += sbr1[y];
      new_sbr2[new_y] += sbr2[y];

//      new_ra_count[new_y] += ra_count[y];
//      new_ra_dist[new_y] += ra_dist[y];
    }

    new_var.infos["SBF"] = join_strand_bias(new_sbf);
    new_var.infos["SBR"] = join_strand_bias(new_sbr);

    new_var.infos["SBF1"] = join_strand_bias(new_sbf1);
    new_var.infos["SBF2"] = join_strand_bias(new_sbf2);
    new_var.infos["SBR1"] = join_strand_bias(new_sbr1);
    new_var.infos["SBR2"] = join_strand_bias(new_sbr2);

//    new_var.infos["RACount"] = join_strand_bias(new_ra_count);
//    new_var.infos["RADist"] = join_strand_bias(new_ra_dist);
  }
}


} // anon namesapce


namespace gyper
{

Variant::Variant() noexcept
  : abs_pos(0)
{}

Variant::Variant(Variant const & var) noexcept
  : abs_pos(var.abs_pos)
  , seqs(var.seqs)
  , calls(var.calls)
  , infos(var.infos)
  , suffix_id(var.suffix_id)
  , is_info_generated(var.is_info_generated)
{}

Variant::Variant(Variant && var) noexcept
  : abs_pos(std::forward<uint32_t>(var.abs_pos))
  , seqs(std::forward<std::vector<std::vector<char> > >(var.seqs))
  , calls(std::forward<std::vector<SampleCall> >(var.calls))
  , infos(std::forward<std::map<std::string, std::string> >(var.infos))
  , suffix_id(std::forward<std::string>(var.suffix_id))
  , is_info_generated(var.is_info_generated)
{}


Variant &
Variant::operator=(Variant const & o) noexcept
{
  abs_pos = o.abs_pos;
  seqs = o.seqs;
  calls = o.calls;
  infos = o.infos;
  suffix_id = o.suffix_id;
  is_info_generated = o.is_info_generated;
  return *this;
}


Variant &
Variant::operator=(Variant && o) noexcept
{
  abs_pos = o.abs_pos;
  seqs = std::move(o.seqs);
  calls = std::move(o.calls);
  infos = std::move(o.infos);
  suffix_id = std::move(o.suffix_id);
  is_info_generated = o.is_info_generated;
  return *this;
}


Variant::Variant(Genotype const & gt)
{
  seqs = graph.get_all_sequences_of_a_genotype(gt);
  abs_pos = gt.id - 1; // -1 cause we always fetch one position back as well
}


Variant::Variant(std::vector<Genotype> const & gts, std::vector<uint16_t> const & hap_calls)
{
  assert(gts.size() > 0);

  for (long i = 0; i < static_cast<long>(hap_calls.size()); ++i)
    seqs.push_back(graph.get_sequence_of_a_haplotype_call(gts, hap_calls[i]));

  abs_pos = gts[0].id - 1; // -1 cause we always fetch one position back as well
}


Variant::Variant(VariantCandidate const & var_candidate) noexcept
  : abs_pos(var_candidate.abs_pos)
  , seqs(var_candidate.seqs)
{}


/******************
 * CLASS MODIFERS *
 ******************/

void
Variant::update_camou_phred(long const ploidy)
{
  // Does not make sense if ploidy is 2 or less
  assert(ploidy > 2);

  for (auto & call : calls)
  {
    auto const & cov = call.coverage;
    assert(cov.size() >= 2);
    long const total_coverage = std::accumulate(cov.begin(), cov.end(), 0l);
    long const cnum = cov.size();
    std::vector<uint8_t> phred;

    if (total_coverage == 0)
    {
      // Not enough data
      phred.resize(cnum * (cnum + 1) / 2, 0);
    }
    else
    {
      phred.resize(cnum * (cnum + 1) / 2, 99);
      assert(phred.size() >= 3);
      phred[0] = 0;
      std::vector<long> normalized_cov;
      normalized_cov.reserve(cnum);

      for (long k{0}; k < cnum; ++k)
        normalized_cov.push_back(cov[k] * ploidy / 2l);

      for (long y{1}; y < cnum; ++y)
      {
        assert(y < static_cast<long>(normalized_cov.size()));

        long constexpr ERROR = 4;
        long const norm_cov = normalized_cov[y];
        long phred00 = norm_cov * ERROR;
        long phred01_or_11 = cov[0];

        {
          long const min_phred = std::min(phred00, phred01_or_11);
          phred00 = std::min(99l, (phred00 - min_phred) * 3l);
          phred01_or_11 = std::min(99l, (phred01_or_11 - min_phred) * 3l);
        }

        if (phred00 > phred[0])
          phred[0] = phred00;

        for (long x{0}; x < cnum; ++x)
        {
          long const index = x <= y ? to_index(x, y) : to_index(y, x);
          assert(index < static_cast<long>(phred.size()));

          if (phred01_or_11 < phred[index])
            phred[index] = phred01_or_11;
        }
      }
    }

    // We should not be changing the size of the phred vector
    assert(call.phred.size() == phred.size());
    call.phred = std::move(phred);
  }
}


void
Variant::generate_infos()
{
  // First check if INFO had already been generated
  if (is_info_generated)
    return;

  assert(seqs.size() >= 2);
  infos["RefLen"] = std::to_string(seqs[0].size());

  // Calculate AC, AN, SeqDepth, MaxVS, ABHom, ABHet, ABHomMulti, ABHetMulti
  std::vector<uint32_t> ac(seqs.size() - 1, 0u);
  std::vector<uint32_t> pass_ac(seqs.size() - 1, 0u);
  assert(ac.size() > 0);
  long an = 0;
  long pass_an = 0;
  uint64_t seqdepth = 0;
  std::pair<uint32_t, uint32_t> het_allele_depth = {0ul, 0ul};   // First is the first call, second is the second call
  std::pair<uint32_t, uint32_t> hom_allele_depth = {0ul, 0ul};   // First is the called allele, second is not the called one

  std::vector<std::pair<uint32_t, uint32_t> > het_multi_allele_depth(seqs.size() - 1, {0u, 0u});
  std::vector<std::pair<uint32_t, uint32_t> > hom_multi_allele_depth(seqs.size(), {0u, 0u});

  // Max variant support and variant support ratio
  std::vector<uint16_t> maximum_variant_support(seqs.size() - 1, 0);
  std::vector<double> maximum_alternative_support_ratio(seqs.size() - 1, 0.0);

  // Num REF/REF, REF/ALT, ALT/ALT
  uint64_t n_ref_ref = 0u;
  uint64_t n_ref_alt = 0u;
  uint64_t n_alt_alt = 0u;

  // Maximum number of alt proper pairs (MaxAltPP)
  uint8_t n_max_alt_proper_pairs = 0u;

  // For calculation how many calls passed
  long n_passed_calls = 0;

  for (auto const & sample_call : calls)
  {
    // Check if the there is any call (i.e. there is a non-zero phred score) and ignore the
    // sample call if that is the case
    if (std::find_if(sample_call.phred.cbegin(), sample_call.phred.cend(), [](uint8_t const val){
        return val != 0;
      }) == sample_call.phred.end())
    {
      continue;
    }

    // MaxAltPP
    n_max_alt_proper_pairs = std::max(n_max_alt_proper_pairs, sample_call.alt_proper_pair_depth);
    uint32_t const total_depth = std::accumulate(sample_call.coverage.cbegin(),
                                                 sample_call.coverage.cend(),
                                                 0u
                                                 );

    // Skip zero (the reference)
    for (uint16_t c = 1; c < sample_call.coverage.size(); ++c)
    {
      maximum_variant_support[c - 1] = std::max(maximum_variant_support[c - 1],
                                                sample_call.coverage[c]
                                                );

      if (total_depth > 0)
      {
        double const current_ratio = static_cast<double>(sample_call.coverage[c]) /
                                     static_cast<double>(total_depth);
        maximum_alternative_support_ratio[c - 1] =
          std::max(maximum_alternative_support_ratio[c - 1],
                   current_ratio
                   );
      }
    }

    // Get the GT
    std::pair<uint16_t, uint16_t> call = sample_call.get_gt_call();
    assert(seqs.size() == sample_call.coverage.size());
    assert(call.first < seqs.size());
    assert(call.second < seqs.size());
    long filter = sample_call.check_filter(sample_call.get_gq());

    if (filter == 0)
      ++n_passed_calls;

    // Restrict to unique calls, that is calls with only one PHRED=0 (@Hannes: removed for now..
    // I don't get why this would make sense)

    //if (std::count(sample_call.phred.begin(), sample_call.phred.end(), 0) == 1)
    {
      // ABHet and ABHom
      {
        // Check if heterozygous
        if (call.first != call.second)
        {
          het_allele_depth.first += sample_call.coverage[call.first];
          het_allele_depth.second += sample_call.coverage[call.second];
        }
        else
        {
          hom_allele_depth.first += sample_call.coverage[call.first];
          hom_allele_depth.second += std::accumulate(sample_call.coverage.cbegin(),
                                                     sample_call.coverage.cend(),
                                                     0ull
                                                     ) - sample_call.coverage[call.first];
        }
      }

      // ABHetMulti and ABHomMulti is calculated on multiallelic sites only
      if (seqs.size() > 2)
      {
        // Check if heterozygous
        if (call.first != call.second)
        {
          // Only consider heterozygous reference calls
          if (call.first == 0u)
          {
            auto const & c = call.second;
            assert((c - 1) < static_cast<long>(het_multi_allele_depth.size()));
            het_multi_allele_depth[c - 1].first += sample_call.coverage[0];
            het_multi_allele_depth[c - 1].second += sample_call.coverage[c];
          }
        }
        else
        {
          auto const & c = call.first;
          assert(c < hom_multi_allele_depth.size());
          hom_multi_allele_depth[c].first += sample_call.coverage[c];
          hom_multi_allele_depth[c].second += std::accumulate(sample_call.coverage.cbegin(),
                                                              sample_call.coverage.cend(),
                                                              0u
                                                              ) - sample_call.coverage[c];
        }
      }
    }

    // SeqDepth and SeqDepth excluding reference allele and reference calls
    if (sample_call.coverage.size() > 0)
      seqdepth += sample_call.get_depth();

    // AN
    {
      an += 2;

      if (filter == 0) // PASS
        pass_an += 2;
    }

    // AC (AF can be derived from AC and AN)
    {
      assert(call.first - 1 < static_cast<long>(ac.size()));
      assert(call.second - 1 < static_cast<long>(ac.size()));

      if (call.first != 0u)
        ++ac[call.first - 1];

      if (call.second != 0u)
        ++ac[call.second - 1];

      if (filter == 0) // PASS
      {
        if (call.first != 0u)
          ++pass_ac[call.first - 1];

        if (call.second != 0u)
          ++pass_ac[call.second - 1];
      }

      if (call.first == 0u && call.second == 0u)
        ++n_ref_ref;
      else if (call.first == 0u || call.second == 0u)
        ++n_ref_alt;
      else
        ++n_alt_alt;
    }
  }   // for sample_call

  // Write MaxAAS
  {
    assert(maximum_variant_support.size() > 0);
    assert(maximum_variant_support.size() == seqs.size() - 1);
    std::stringstream max_vs_ss;
    max_vs_ss << static_cast<uint16_t>(maximum_variant_support[0]);

    for (uint32_t e = 1; e < maximum_variant_support.size(); ++e)
      max_vs_ss << ',' << static_cast<uint16_t>(maximum_variant_support[e]);

    infos["MaxAAS"] = max_vs_ss.str();
  }

  // Write MaxAASR
  {
    assert(maximum_alternative_support_ratio.size() > 0);
    assert(maximum_alternative_support_ratio.size() == seqs.size() - 1);
    std::stringstream ss_max_asr;
    ss_max_asr.precision(4);
    ss_max_asr << maximum_alternative_support_ratio[0];

    for (uint32_t e = 1; e < maximum_alternative_support_ratio.size(); ++e)
      ss_max_asr << ',' << maximum_alternative_support_ratio[e];

    infos["MaxAASR"] = ss_max_asr.str();
  }

  // Write MaxAltPP
  if (is_sv())
  {
    infos["MaxAltPP"] = std::to_string(static_cast<uint16_t>(n_max_alt_proper_pairs));
  }

  // Write AC
  {
    assert(ac.size() > 0);
    assert(ac.size() == seqs.size() - 1);
    std::ostringstream ac_ss;
    ac_ss << ac[0];

    for (uint32_t e = 1; e < ac.size(); ++e)
      ac_ss << ',' << ac[e];

    infos["AC"] = ac_ss.str();
  }

  // Write AN
  infos["AN"] = std::to_string(an);

  // Write PASS_AC
  {
    assert(pass_ac.size() > 0);
    assert(pass_ac.size() == seqs.size() - 1);
    std::ostringstream ac_ss;
    ac_ss << pass_ac[0];

    for (uint32_t e = 1; e < pass_ac.size(); ++e)
      ac_ss << ',' << pass_ac[e];

    infos["PASS_AC"] = ac_ss.str();
  }

  // Write PASS_AN
  infos["PASS_AN"] = std::to_string(pass_an);

  // Write NRP
  {
    bool is_any_hq_alt = std::any_of(pass_ac.begin(), pass_ac.end(), [](long ac){
        return ac > 0;
      });

    if (!is_any_hq_alt)
      infos["NRP"] = "";
  }

  // Write PASS_ratio
  if (calls.size() > 0)
  {
    double pr = static_cast<double>(n_passed_calls) / static_cast<double>(calls.size());
    infos["PASS_ratio"] = std::to_string(pr);
  }

  // Write Num REF/REF, REF/ALT, ALT/ALT, homozygous and heterozygous calls
  {
    std::ostringstream n_gt;
    n_gt << n_ref_ref << "," << n_ref_alt << "," << n_alt_alt;
    infos["NHomRef"] = std::to_string(n_ref_ref);
    infos["NHet"] = std::to_string(n_ref_alt);
    infos["NHomAlt"] = std::to_string(n_alt_alt);
  }

  // Write SeqDepth
  {
    infos["SeqDepth"] = std::to_string(seqdepth);
  }

  // Write ABHet
  {
    std::stringstream ss_abhet;
    ss_abhet.precision(4);
    uint32_t const total_het_depth = het_allele_depth.first + het_allele_depth.second;

    if (total_het_depth > 0)
    {
      ss_abhet << (static_cast<double>(het_allele_depth.second) /
                   static_cast<double>(total_het_depth)
                   );
    }
    else
    {
      ss_abhet << "-1";
    }

    infos["ABHet"] = ss_abhet.str();
  }

  // Write ABHom
  {
    std::stringstream ss_abhom;
    ss_abhom.precision(4);
    uint32_t const total_hom_depth = hom_allele_depth.first + hom_allele_depth.second;

    if (total_hom_depth > 0)
    {
      ss_abhom << (static_cast<double>(hom_allele_depth.first) /
                   static_cast<double>(total_hom_depth));
    }
    else
    {
      ss_abhom << "-1";
    }

    infos["ABHom"] = ss_abhom.str();
  }

  // Write Strand Bias (SB)
  {
    uint32_t total_f = get_accumulated_strand_bias(infos, "SBF");
    uint32_t total_r = get_accumulated_strand_bias(infos, "SBR");

    std::stringstream ss_sb;
    ss_sb.precision(4);

    if (total_f + total_r == 0)
      ss_sb << "-1";
    else
      ss_sb << (static_cast<double>(total_f) / static_cast<double>(total_f + total_r));

    infos["SB"] = ss_sb.str();
  }

  // Write Alternative Strand Bias (SBAlt)
  {
    uint32_t total_f = get_accumulated_alt_strand_bias(infos, "SBF");
    uint32_t total_r = get_accumulated_alt_strand_bias(infos, "SBR");

    std::stringstream ss_sb;
    ss_sb.precision(4);

    if (total_f + total_r == 0)
      ss_sb << "-1";
    else
      ss_sb << (static_cast<double>(total_f) / static_cast<double>(total_f + total_r));

    infos["SBAlt"] = ss_sb.str();
  }

  if (seqs.size() > 2)
  {
    // Write ABHetMulti
    {
      assert(het_multi_allele_depth.size() > 0);
      assert(het_multi_allele_depth.size() == seqs.size() - 1);
      std::stringstream ss_abhet;
      ss_abhet.precision(4);

      for (unsigned i = 0; i < het_multi_allele_depth.size(); ++i)
      {
        if (i > 0)
          ss_abhet << ",";

        auto const & het_depth_pair = het_multi_allele_depth[i];
        uint32_t const total_het_depth = het_depth_pair.first + het_depth_pair.second;

        if (total_het_depth > 0)
        {
          ss_abhet << (static_cast<double>(het_depth_pair.second) /
                       static_cast<double>(total_het_depth)
                       );
        }
        else
        {
          ss_abhet << "-1";
        }
      }

      infos["ABHetMulti"] = ss_abhet.str();
    }

    // Write ABHomMulti
    {
      assert(hom_multi_allele_depth.size() > 0);
      assert(hom_multi_allele_depth.size() == seqs.size());
      std::stringstream ss_abhom;
      ss_abhom.precision(4);

      for (unsigned i = 0; i < hom_multi_allele_depth.size(); ++i)
      {
        if (i > 0)
          ss_abhom << ",";

        auto const & hom_depth_pair = hom_multi_allele_depth[i];
        uint32_t const total_hom_depth = hom_depth_pair.first + hom_depth_pair.second;

        if (total_hom_depth > 0)
        {
          ss_abhom << (static_cast<double>(hom_depth_pair.first) /
                       static_cast<double>(total_hom_depth));
        }
        else
        {
          ss_abhom << "-1";
        }
      }

      infos["ABHomMulti"] = ss_abhom.str();
    }
  }

  // Write VarType
  {
    infos["VarType"] = this->determine_variant_type();
  }

  // Calculate QD
  {
    std::stringstream ss;
    ss.precision(4);
    ss << get_qual_by_depth();
    infos["QD"] = ss.str();
  }

  // Calculate MQ
  {
    auto find_mapq_squared_it = infos.find("MQsquared");

    if (find_mapq_squared_it != infos.end())
    {
      double mapq_squared = std::stoull(find_mapq_squared_it->second);

      if (seqdepth > 0)
      {
        double mapq = std::sqrt(mapq_squared / static_cast<double>(seqdepth));
        infos["MQ"] = std::to_string(static_cast<uint64_t>(mapq + 0.499999));
      }
      else
      {
        infos["MQ"] = "255";
      }
    }
  }

  is_info_generated = true;
}


bool
Variant::is_sv() const
{
  // If the variant has an SV breakpoint, skip trimming
  for (long s = 1; s < static_cast<long>(seqs.size()); ++s)
  {
    auto const & seq = seqs[s];

    if (seq.size() < 5)
      continue;

    if (seq[0] == '<' || (seq.size() > 100 && std::find(seq.begin(), seq.end(), '<') != seq.end()))
      return true; // Found an SV
  }

  return false;
}


void
Variant::remove_common_prefix(bool const keep_one_match)
{
  gyper::remove_common_prefix(abs_pos, seqs, keep_one_match);
}


void
Variant::trim_sequences(bool const keep_one_match)
{
  add_base_in_front();

  if (!is_sv())
    remove_common_suffix(seqs);

  gyper::remove_common_prefix(abs_pos, seqs, keep_one_match);
}


bool
Variant::add_base_in_front(bool const add_N)
{
  uint32_t abs_pos_copy = abs_pos;
  uint32_t new_abs_pos = abs_pos - 1;
  std::vector<char> first_base = graph.get_generated_reference_genome(new_abs_pos, abs_pos_copy);

  if (first_base.size() != 1 || abs_pos_copy != abs_pos || new_abs_pos != abs_pos - 1)
    return false; // The base in front could not be extracted

  if (!add_N && first_base[0] == 'N')
    return false;

  // Insert the new first base in front of all the sequences
  for (auto & seq : seqs)
  {
    // Do not add base on missing allele
    if (seq.size() == 0 || seq.size() > 1 || seq[0] != '*')
    {
      seq.insert(seq.begin(), first_base[0]);
    }
  }

  // Update the abs pos
  --abs_pos;

  return true;
}


/*
bool
Variant::add_multiple_bases_in_front(long const num, bool const add_N)
{
  long abs_pos_copy = abs_pos;
  long new_abs_pos = std::max(0, abs_pos - num);
  std::vector<char> bases = graph.get_generated_reference_genome(new_abs_pos, abs_pos_copy);

  if (first_base.size() != 1 || abs_pos_copy != abs_pos || new_abs_pos != abs_pos - 1)
    return false; // The base in front could not be extracted

  if (!add_N && first_base[0] == 'N')
    return false;

  // Insert the new first base in front of all the sequences
  for (auto & seq : seqs)
  {
    // Do not add base on missing allele
    if (seq.size() == 0 || seq.size() > 1 || seq[0] != '*')
    {
      seq.insert(seq.begin(), first_base[0]);
    }
  }

  // Update the abs pos
  --abs_pos;

  return true;
}
*/


bool
Variant::add_base_in_back(bool const add_N)
{
  assert(seqs.size() >= 1);
  uint32_t abs_pos_copy = abs_pos + static_cast<uint32_t>(seqs[0].size());
  uint32_t abs_pos_end = abs_pos_copy + 1;
  std::vector<char> last_base = graph.get_generated_reference_genome(abs_pos_copy, abs_pos_end);

  if (last_base.size() != 1 || abs_pos_copy != (abs_pos + seqs[0].size()) ||
      abs_pos_end != (abs_pos + seqs[0].size() + 1))
    return false; // The base in back could not be extracted";

  if (!add_N && last_base[0] == 'N')
    return false;

  // Insert the new first base in back of all the sequences
  for (auto & seq : seqs)
  {
    // Do not add base on missing allele
    //if (seq.size() == 0 || seq.size() > 1 || seq[0] != '*')
    {
      seq.push_back(last_base[0]);
    }
  }

  return true;
}


void
Variant::expanded_normalized()
{
  // First normalize
  normalize();

  // Then expand on the right side if it is an indel
  if (!is_snp_or_snps())
  {
    long i = 0;
    bool is_done = false;

    while (!is_done && add_base_in_back(false))
    {
      ++i;
      assert(i < static_cast<long>(seqs[0].size()));
      char const ref_base = seqs[0][i];

      for (long s = 1; s < static_cast<long>(seqs.size()); ++s)
      {
        if (seqs[s][i] != ref_base)
        {
          is_done = true;
          break;
        }
      }
    }
  }
}


void
Variant::normalize()
{
  if (seqs.size() < 2)
    return;

  auto const & ref = seqs[0];

  // Check if some sequence is of length 0 or prefix is not matching
  for (long i = 0; i < static_cast<long>(seqs.size()); ++i)
  {
    auto const & seq = seqs[i];

    if (seq.size() == 0 || seq[0] != ref[0])
      return;

    // If a sequence is the same as the ref then we will loop forever
    if (i > 0 && seq == ref)
      return;
  }

  remove_common_suffix(seqs);

  // A lambda function which checks if all last bases matches
  auto all_last_bases_match =
    [&ref, this]()
    {
      for (long i = 1; i < static_cast<long>(this->seqs.size()); ++i)
      {
        if (this->seqs[i].back() != ref.back())
          return false;
      }

      return true;
    };

  while (all_last_bases_match())
  {
    bool const success_adding_base = add_base_in_front();

    if (not success_adding_base)
      break;

    // Remove the last base
    remove_common_suffix(seqs);
  }

  gyper::remove_common_prefix(abs_pos, seqs, false); // Don't keep one match
}


/***********************
 * VARIANT INFORMATION *
 ***********************/
bool
Variant::is_normalized() const
{
  Variant new_var;
  new_var.abs_pos = this->abs_pos;
  new_var.seqs = this->seqs;
  new_var.normalize();
  return new_var == *this;
}


bool
Variant::is_snp_or_snps() const
{
  // Check if all the alleles are of the same size
  assert(seqs.size() >= 2);
  return std::find_if(seqs.begin() + 1,
                      seqs.end(),
                      [&](std::vector<char> const & seq){
      return seq.size() != seqs[0].size();
    }
                      ) == seqs.end();
}


bool
Variant::is_with_matching_first_bases() const
{
  assert(seqs.size() > 0);
  assert(seqs[0].size() > 0);
  char const first_base = seqs[0][0];

  for (std::size_t i = 1; i < seqs.size(); ++i)
  {
    assert(seqs[i].size() > 0);

    if (seqs[i][0] != first_base)
      return false;
  }

  return true;
}


uint64_t
Variant::get_seq_depth() const
{
  uint64_t seqdepth = 0;

  for (auto const & sample_call : calls)
  {
    seqdepth += std::accumulate(sample_call.coverage.cbegin(),
                                sample_call.coverage.cend(),
                                static_cast<uint32_t>(sample_call.ambiguous_depth)
                                );
  }

  return seqdepth;
}


uint64_t
Variant::get_seq_depth_of_allele(uint16_t const allele_id) const
{
  uint64_t seqdepth = 0;

  for (auto const & sample_call : calls)
  {
    assert(allele_id < sample_call.coverage.size());
    seqdepth += sample_call.coverage[allele_id];
  }

  return seqdepth;
}


std::vector<uint64_t>
Variant::get_seq_depth_of_all_alleles() const
{
  std::vector<uint64_t> seq_depths;
  seq_depths.reserve(this->seqs.size());
  long const num_seqs = this->seqs.size();

  for (long i = 0; i < num_seqs; ++i)
  {
    seq_depths.push_back(this->get_seq_depth_of_allele(static_cast<uint16_t>(i)));
  }

  return seq_depths;
}


std::string
Variant::print() const
{
  std::stringstream os;
  auto contig_pos = absolute_pos.get_contig_position(this->abs_pos, gyper::graph.contigs);
  os << contig_pos.first << "\t" << contig_pos.second;

  if (this->seqs.size() > 0)
    os << "\t" << std::string(this->seqs[0].begin(), this->seqs[0].end());

  if (this->seqs.size() > 1)
    os << "\t" << std::string(this->seqs[1].begin(), this->seqs[1].end());
  else
    os << "\tN";

  for (unsigned i = 2; i < this->seqs.size(); ++i)
    os << "," << std::string(this->seqs[i].begin(), this->seqs[i].end());

  for (auto const & sample : this->calls)
  {
    os << "\t";

    if (sample.phred.size() > 0)
      os << static_cast<uint16_t>(sample.phred[0]);

    for (unsigned p = 1; p < sample.phred.size(); ++p)
      os << "," << static_cast<uint16_t>(sample.phred[p]);
  }

  return os.str();
}


std::string
Variant::determine_variant_type() const
{
  assert(seqs.size() >= 2);

  // Count the number of alleles which are not of size 1
  std::size_t num_non_ones = 0;

  // Check if the variant is an SV
  gyper::SVTYPE sv_type = NOT_SV;

  for (auto it = seqs.begin(); it != seqs.end(); ++it)
  {
    assert(it->size() != 0);

    if (it->size() > 1)
    {
      if (it->size() > 4 && (*it)[0] == '<')
      {
        std::string const type(it->begin() + 1, it->begin() + 4);

        if (type == "DEL" && (sv_type == NOT_SV || sv_type == DEL))
          sv_type = DEL;
        else if (type == "DUP" && (sv_type == NOT_SV || sv_type == DUP))
          sv_type = DUP;
        else if (type == "INS" && (sv_type == NOT_SV || sv_type == INS))
          sv_type = INS;
        else
          sv_type = OTHER;
      }
      else if (std::find_if(it->begin(), it->end(), [](char c){
          return c == '[' || c == ']';
        }) != it->end())
      {
        if (sv_type == NOT_SV || sv_type == BND)
          sv_type = BND;
        else
          sv_type = OTHER;
      }
      else
      {
        ++num_non_ones;
      }
    }
  }

  // Check if the variant is a SNP
  if (sv_type != NOT_SV)
  {
    switch (sv_type)
    {
    case DEL: return "DG";

    case DUP: return "UG";

    case INS: return "FG";

    case INV: return "NG";

    case BND: return "OG";

    default: return "TG";
    }
  }
  else if (num_non_ones == 0)
  {
    return "SG";
  }
  else if (seqs.size() - num_non_ones == 1) // Check if the variant is an indel
  {
    return "IG";
  }
  else if (seqs.size() - num_non_ones == 2 && seqs[seqs.size() - 1].size() == 1 && seqs[seqs.size() - 1][0] == '*')
  {
    // Indel with missing allele
    return "IG";
  }
  else // Otherwise the variant is complex
  {
    return "XG";
  }
}


uint64_t
Variant::get_qual() const
{
  uint64_t variant_qual = 0;

  for (auto const & sample_call : calls)
  {
    assert(sample_call.phred.size() > 0);
    variant_qual += sample_call.phred[0];
  }

  return variant_qual;
}


double
Variant::get_qual_by_depth() const
{
  long total_qual = 0;
  long total_depth = 0;

  for (auto const & sample_call : calls)
  {
    // Skip homozygous ref calls
    if (sample_call.phred[0] > 0)
    {
      auto const & cov = sample_call.coverage;
      total_qual += static_cast<long>(sample_call.phred[0]);
      long const depth = std::min(10l, std::accumulate(cov.begin() + 1, cov.end(), 0l));
      total_depth += depth;
    }
  }

  if (total_depth == 0)
    return 0.0;
  else
    return static_cast<double>(total_qual) / static_cast<double>(total_depth);
}


std::vector<Variant>
make_biallelic(Variant && var)
{
  std::vector<Variant> out;

  if (var.seqs.size() == 2)
  {
    // Already biallelic
    out.push_back(std::move(var));
    return out;
  }

  for (long a = 1; a < static_cast<long>(var.seqs.size()); ++a)
  {
    Variant new_var;
    new_var.seqs.push_back(var.seqs[0]);
    new_var.seqs.push_back(var.seqs[a]);
    new_var.abs_pos = var.abs_pos; // Add the pos
    new_var.infos = var.infos; // Copy the INFOs
    new_var.suffix_id = var.suffix_id; // Copy the suffix ID
    std::vector<uint16_t> old_phred_to_new_phred(var.seqs.size(), 0);
    old_phred_to_new_phred[a] = 1;

    for (long i = 0; i < static_cast<long>(var.calls.size()); ++i)
    {
      auto const & call = var.calls[i];

      SampleCall new_sample_call;
      new_sample_call.phred = std::vector<uint8_t>(3, 255u);
      new_sample_call.coverage = std::vector<uint16_t>(2, 0u);
      new_sample_call.ambiguous_depth = call.ambiguous_depth;
      new_sample_call.ref_total_depth = call.ref_total_depth;
      new_sample_call.alt_total_depth = call.alt_total_depth;
      new_sample_call.alt_proper_pair_depth = call.alt_proper_pair_depth;

      for (uint32_t y = 0; y < var.seqs.size(); ++y)
      {
        uint32_t const new_y = old_phred_to_new_phred[y];

        for (uint32_t x = 0; x <= y; ++x)
        {
          uint32_t const index = static_cast<uint32_t>(to_index(x, y));
          uint32_t const new_x = old_phred_to_new_phred[x];
          long new_index;

          if (new_x > new_y)
            new_index = to_index(new_y, new_x);
          else
            new_index = to_index(new_x, new_y);

          new_sample_call.phred[new_index] = std::min(new_sample_call.phred[new_index], call.phred[index]);
        }

        // Check overflow of coverage
        if (static_cast<uint32_t>(new_sample_call.coverage[new_y]) +
            static_cast<uint32_t>(call.coverage[y]) < 0xFFFFul)
        {
          new_sample_call.coverage[new_y] += call.coverage[y];
        }
        else
        {
          new_sample_call.coverage[new_y] = 0xFFFFul;
        }
      }

      new_var.calls.push_back(std::move(new_sample_call));
    }

    update_strand_bias(var.seqs.size(), new_var.seqs.size(), old_phred_to_new_phred, var, new_var);
    out.push_back(std::move(new_var));
  }

  return out;
}


// Variants functions
std::vector<Variant>
break_down_variant(Variant && var,
                   long const reach,
                   bool const is_no_variant_overlapping,
                   bool const is_all_biallelic)
{
  std::vector<Variant> broken_down_vars;

  if (Options::const_instance()->no_decompose ||
      (var.seqs.size() == 2 &&
       std::any_of(var.seqs[1].begin(), var.seqs[1].end(), [](char const c){
      return c == '<' || c == '[' || c == ']';
    })))
  {
    broken_down_vars.push_back(std::move(var)); // Don't break down SVs or if "no_decompose" was given
    return broken_down_vars;
  }

  bool const all_same_size =
    std::find_if(var.seqs.begin() + 1, var.seqs.end(), [&](std::vector<char> const & seq){
      return var.seqs[0].size() != seq.size();
    }) == var.seqs.end();

  if (all_same_size)
  {
    // We need to make sure there is a matching first base
    if (not var.is_with_matching_first_bases())
    {
      if (!var.add_base_in_front())
      {
        // Could not add a first base. Add N
        for (auto & seq : var.seqs)
          seq.insert(seq.begin(), 'N');
      }
    }

    std::vector<Variant> new_broken_down_snps = break_multi_snps(std::move(var));

    std::move(new_broken_down_snps.begin(),
              new_broken_down_snps.end(),
              std::back_inserter(broken_down_vars)
              );
  }
  else if (!is_no_variant_overlapping)
  {
    // Use the skyr
    BOOST_LOG_TRIVIAL(debug) << "Using the skyr";
    std::vector<Variant> new_broken_down_vars = break_down_skyr(std::move(var), reach);
    BOOST_LOG_TRIVIAL(debug) << "skyr finished.";

    std::move(new_broken_down_vars.begin(),
              new_broken_down_vars.end(),
              std::back_inserter(broken_down_vars));
  }
  else
  {
    broken_down_vars.push_back(std::move(var));
  }

  if (is_all_biallelic)
  {
    // Make everything biallelic
    std::vector<Variant> broken_down_vars2;

    for (auto && var : broken_down_vars)
    {
      auto new_vars = make_biallelic(std::move(var));
      std::move(new_vars.begin(), new_vars.end(), std::back_inserter(broken_down_vars2));
    }

    // replace origin broken down vars with biallelic variants
    broken_down_vars = std::move(broken_down_vars2);
  }

  return broken_down_vars;
}


std::vector<Variant>
extract_sequences_from_aligned_variant(Variant const && var, std::size_t const THRESHOLD)
{
  std::vector<Variant> new_vars;
  uint32_t const original_pos = var.abs_pos;
  assert(var.seqs.size() >= 2);
  assert(var.seqs[0].size() > 0);
  assert(var.seqs[0].size() == var.seqs[1].size());
  char first_base = var.seqs[0][0];
  uint32_t pos = var.abs_pos;

  // match_length == -1 means that we have not found a mismatch yet
  long match_length = -1;
  unsigned ref_gaps = 0;

  Variant new_var;
  new_var.abs_pos = pos;
  new_var.seqs =
    std::vector<std::vector<char> >(var.seqs.size(), std::vector<char>(1, first_base));

  new_var.infos = var.infos; // Copy INFOs
//  new_var.phase = var.phase; // Copy phase
  new_var.suffix_id = var.suffix_id; // Copy suffix ID

  // Lambda function which handles any matches found
  auto matches_handle =
    [](std::vector<Variant> & new_vars, Variant && new_var) -> void
    {
      // Remove sequences with Ns
      assert(new_var.seqs.size() > 1);
      new_var.seqs.erase(std::remove_if(new_var.seqs.begin() + 1,
                                        new_var.seqs.end(),
                                        [](std::vector<char> const & seq){
        return std::find(seq.begin(), seq.end(), 'N') != seq.end();
      }), new_var.seqs.end());

      if (new_var.seqs.size() <= 1 ||
          std::find(new_var.seqs[0].begin(), new_var.seqs[0].end(), 'N') != new_var.seqs[0].end())
      {
        // Do not do anything if there is only the reference or the reference has Ns
        return;
      }

      new_var.trim_sequences(false);  // Keep one match
      new_vars.push_back(std::move(new_var));
    };

  std::vector<char> const & reference = var.seqs[0];

  for (long i = 1; i < static_cast<long>(reference.size()); ++i)
  {
    assert(new_var.seqs.size() >= 2);

    if (reference[i] == '-')
      ++ref_gaps;
    else
      new_var.seqs[0].push_back(reference[i]);

    bool all_match = true;

    for (int a = 1; a < static_cast<long>(new_var.seqs.size()); ++a)
    {
      assert(a < static_cast<int>(var.seqs.size()));

      if (var.seqs[a][i] != '-')
        new_var.seqs[a].push_back(var.seqs[a][i]);

      if (var.seqs[a][i] != reference[i])
        all_match = false;
    }

    if (all_match)
    {
      if (match_length >= 0)
        ++match_length;
    }
    else
    {
      match_length = 0;
    }

    if (match_length >= static_cast<long>(THRESHOLD))
    {
      //first_base = var.seqs[0][i];
      find_variant_sequences(new_var, var);
      matches_handle(new_vars, std::move(new_var));   // Uses new_var, so don't move this anywhere else
      match_length = -1;
      //new_var = Variant();   // Create a new one
      new_var.abs_pos = original_pos + i - ref_gaps + 1;
      new_var.seqs = std::vector<std::vector<char> >(var.seqs.size());
      //new_var.abs_pos = original_pos + i - ref_gaps;
      //new_var.seqs =
      //  std::vector<std::vector<char> >(var.seqs.size(), std::vector<char>(1, first_base));

//      new_var.infos = var.infos;
//      new_var.phase = var.phase;
      new_var.suffix_id = var.suffix_id;
    }
  }

  assert(new_var.seqs.size() >= 2);

  // Process any leftovers which might be remaining
  if (new_var.seqs[0].size() > 0)
  {
    find_variant_sequences(new_var, var);

    if (new_var.seqs.size() >= 2)
      matches_handle(new_vars, std::move(new_var));
  }

  return new_vars;
}


void
find_variant_sequences(gyper::Variant & new_var, gyper::Variant const & old_var)
{
  using gyper::get_strand_bias;
  using gyper::join_strand_bias;
  using gyper::to_index;
  using gyper::SampleCall;

  assert(new_var.calls.size() == 0);
  auto & seqs = new_var.seqs;

  assert(seqs.size() >= 2);
  assert(seqs[0].size() > 0);
  assert(seqs[1].size() > 0);

  std::vector<uint16_t> old_phred_to_new_phred(1, 0);
  std::vector<std::vector<char> > new_seqs(1, seqs[0]);

  for (std::size_t a = 1; a < seqs.size(); ++a)
  {
    auto find_it = std::find(new_seqs.begin(), new_seqs.end(), seqs[a]);

    if (find_it == new_seqs.end() and seqs[a].size() != 0)
    {
      old_phred_to_new_phred.push_back(static_cast<uint16_t>(new_seqs.size()));
      new_seqs.push_back(seqs[a]);
    }
    else
    {
      old_phred_to_new_phred.push_back(
        static_cast<uint16_t>(std::distance(new_seqs.begin(), find_it))
        );
    }
  }

  assert(old_phred_to_new_phred.size() == seqs.size());

  if (new_seqs.size() <= 1)
  {
    // No variant found
    seqs.clear();
    return;
  }

  for (std::size_t i = 0; i < old_var.calls.size(); ++i)
  {
    auto const & old_call = old_var.calls[i];
    SampleCall new_call;

    // phred
    {
      assert(old_call.phred.size() == (old_var.seqs.size() * (old_var.seqs.size() + 1) / 2));
      new_call.phred = std::vector<uint8_t>(new_seqs.size() * (new_seqs.size() + 1) / 2, 255u);

      for (uint16_t y = 0; y < old_var.seqs.size(); ++y)
      {
        for (uint16_t x = 0; x <= y; ++x)
        {
          uint16_t const index = static_cast<uint16_t>(to_index(x, y));
          uint16_t new_x = old_phred_to_new_phred[x];
          uint16_t new_y = old_phred_to_new_phred[y];

          //if (old_var.phase.size() > 0 && old_call.phred[index] == 0)
          //{
          //  assert(i < old_var.phase.size());
          //
          //  if (new_x > new_y)
          //    new_var.phase[i] = !old_var.phase[i];
          //}

          if (new_x > new_y)
            std::swap<uint16_t>(new_x, new_y);

          uint16_t const new_index = static_cast<uint16_t>(to_index(new_x, new_y));
          assert(new_index < new_call.phred.size());
          assert(index < old_call.phred.size());
          new_call.phred[new_index] = std::min(new_call.phred[new_index], old_call.phred[index]);
        }
      }
    }

    // coverage
    {
      assert(old_call.coverage.size() == old_var.seqs.size());
      new_call.coverage = std::vector<uint16_t>(new_seqs.size(), 0u);

      for (uint32_t y = 0; y < old_var.seqs.size(); ++y)
      {
        uint32_t const new_y = old_phred_to_new_phred[y];

        // Check overflow of coverage
        if (static_cast<uint32_t>(new_call.coverage[new_y]) +
            static_cast<uint32_t>(old_call.coverage[y]) < 0xFFFFul)
        {
          new_call.coverage[new_y] += old_call.coverage[y];
        }
        else
        {
          new_call.coverage[new_y] = 0xFFFFul;
        }
      }
    }

    // other
    {
      new_call.ambiguous_depth = old_call.ambiguous_depth;
      new_call.ref_total_depth = old_call.ref_total_depth;
      new_call.alt_total_depth = old_call.alt_total_depth;
      new_call.alt_proper_pair_depth = old_call.alt_proper_pair_depth;
    }

    // Update strand bias
    update_strand_bias(old_var.seqs.size(),
                       new_seqs.size(),
                       old_phred_to_new_phred,
                       old_var,
                       new_var);

    new_var.calls.push_back(std::move(new_call));
  }

  new_var.seqs = std::move(new_seqs);
}


std::vector<Variant>
break_multi_snps(Variant const && var)
{
  uint32_t const pos = var.abs_pos;
  std::vector<std::vector<char> > const & seqs = var.seqs;
  std::vector<Variant> new_vars;

  for (long j = 0; j < static_cast<long>(seqs[0].size()); ++j)
  {
    std::vector<char> new_seqs(1, seqs[0][j]);
    std::vector<uint16_t> old_phred_to_new_phred(1, 0);

    for (long k = 1; k < static_cast<long>(seqs.size()); ++k)
    {
      auto find_it = std::find(new_seqs.begin(), new_seqs.end(), seqs[k][j]);

      if (find_it == new_seqs.end())
      {
        assert(new_seqs.size() <= std::numeric_limits<uint16_t>::max());
        old_phred_to_new_phred.push_back(static_cast<uint16_t>(new_seqs.size()));
        new_seqs.push_back(seqs[k][j]);
      }
      else
      {
        long const dist = std::distance(new_seqs.begin(), find_it);
        assert(dist >= 0);
        assert(dist <= 0xFFFF);
        old_phred_to_new_phred.push_back(static_cast<uint16_t>(dist));
      }
    }

    if (new_seqs.size() == 1)
      continue; // All sequences are the same => No SNP at this location

    // Create a variant record for the SNP
    Variant new_var;

    // Add the sequences of the SNP
    for (auto seq_it = new_seqs.begin(); seq_it != new_seqs.end(); ++seq_it)
      new_var.seqs.push_back(std::vector<char>(1, *seq_it));

    new_var.abs_pos = pos + j; // Add the pos of the SNP
    new_var.infos = var.infos; // Copy the INFOs
    //new_var.phase = var.phase; // Copy the old phase
    new_var.suffix_id = var.suffix_id; // Copy the suffix ID

    for (long i = 0; i < static_cast<long>(var.calls.size()); ++i)
    {
      auto const & call = var.calls[i];
      SampleCall new_sample_call;
      new_sample_call.phred = std::vector<uint8_t>(new_var.seqs.size() * (new_var.seqs.size() + 1) / 2, 255u);
      new_sample_call.coverage = std::vector<uint16_t>(new_var.seqs.size(), 0u);
      new_sample_call.ambiguous_depth = call.ambiguous_depth;
      new_sample_call.ref_total_depth = call.ref_total_depth;
      new_sample_call.alt_total_depth = call.alt_total_depth;
      new_sample_call.alt_proper_pair_depth = call.alt_proper_pair_depth;

      for (uint32_t y = 0; y < seqs.size(); ++y)
      {
        for (uint32_t x = 0; x <= y; ++x)
        {
          uint32_t const index = static_cast<uint32_t>(to_index(x, y));
          uint32_t new_y = old_phred_to_new_phred[y];
          uint32_t new_x = old_phred_to_new_phred[x];

          if (new_x > new_y)
            std::swap<uint32_t>(new_x, new_y);

          long const new_index = to_index(new_x, new_y);
          new_sample_call.phred[new_index] = std::min(new_sample_call.phred[new_index], call.phred[index]);
        }

        uint32_t new_y = old_phred_to_new_phred[y];

        // Check overflow of coverage
        if (static_cast<uint32_t>(new_sample_call.coverage[new_y]) +
            static_cast<uint32_t>(call.coverage[y]) < 0xFFFFul)
        {
          new_sample_call.coverage[new_y] += call.coverage[y];
        }
        else
        {
          new_sample_call.coverage[new_y] = 0xFFFFul;
        }
      }

      new_var.calls.push_back(std::move(new_sample_call));
    }

    if (var.infos.count("SBF1") == 1)
      update_strand_bias(seqs.size(), new_var.seqs.size(), old_phred_to_new_phred, var, new_var);

    new_vars.push_back(std::move(new_var));
  }

  return new_vars;
}


std::vector<Variant>
break_down_skyr(Variant && var, long const reach)
{
  std::vector<Variant> new_vars;

  long constexpr OPTIMAL_EXTRA = 50l;
  bool const use_asterisks = !Options::instance()->no_asterisks;

  long const extra_bases_before = use_asterisks ?
                                  std::min(OPTIMAL_EXTRA, static_cast<long>(var.abs_pos) - reach - 1l) :
                                  OPTIMAL_EXTRA;

  long constexpr extra_bases_after = OPTIMAL_EXTRA;

  for (long i = 0; i < extra_bases_before && var.add_base_in_front(false); ++i)
  {}

  for (long i = 0; i < extra_bases_after && var.add_base_in_back(false); ++i)
  {}

  auto const & seqs = var.seqs;
  paw::Skyr skyr(var.seqs);

  bool constexpr is_normalize = true; // The events must the normalized
  skyr.find_all_edits(is_normalize);
  skyr.find_variants_from_edits();
  skyr.populate_variants_with_calls(use_asterisks);

  for (auto & new_edit : skyr.vars)
  {
    Variant new_var;

    // Add the sequences of the SNP
    for (auto const & edit_seq : new_edit.seqs)
      new_var.seqs.push_back(std::vector<char>(edit_seq.begin(), edit_seq.end()));

    new_var.abs_pos = var.abs_pos + new_edit.pos; // Add the pos of the SNP

    if (!new_var.is_snp_or_snps())
      new_var.add_base_in_front(true); // Add N is true

    new_var.infos = var.infos; // Copy the INFOs
    new_var.suffix_id = var.suffix_id; // Copy the suffix ID

    auto const & old_phred_to_new_phred = new_edit.calls;

    for (long i = 0; i < static_cast<long>(var.calls.size()); ++i)
    {
      auto const & call = var.calls[i];
      SampleCall new_sample_call;
      new_sample_call.phred = std::vector<uint8_t>(new_var.seqs.size() * (new_var.seqs.size() + 1) / 2, 255u);
      new_sample_call.coverage = std::vector<uint16_t>(new_var.seqs.size(), 0u);
      new_sample_call.ambiguous_depth = call.ambiguous_depth;
      new_sample_call.ref_total_depth = call.ref_total_depth;
      new_sample_call.alt_total_depth = call.alt_total_depth;
      new_sample_call.alt_proper_pair_depth = call.alt_proper_pair_depth;

      for (uint32_t y = 0; y < seqs.size(); ++y)
      {
        for (uint32_t x = 0; x <= y; ++x)
        {
          uint32_t const index = static_cast<uint32_t>(to_index(x, y));
          uint32_t new_y = old_phred_to_new_phred[y];
          uint32_t new_x = old_phred_to_new_phred[x];

          if (new_x > new_y)
            std::swap<uint32_t>(new_x, new_y);

          long const new_index = to_index(new_x, new_y);
          new_sample_call.phred[new_index] = std::min(new_sample_call.phred[new_index], call.phred[index]);
        }

        uint32_t new_y = old_phred_to_new_phred[y];

        // Check overflow of coverage
        if (static_cast<uint32_t>(new_sample_call.coverage[new_y]) +
            static_cast<uint32_t>(call.coverage[y]) < 0xFFFFul)
        {
          new_sample_call.coverage[new_y] += call.coverage[y];
        }
        else
        {
          new_sample_call.coverage[new_y] = 0xFFFFul;
        }
      }

      new_var.calls.push_back(std::move(new_sample_call));
    }

    update_strand_bias(seqs.size(), new_var.seqs.size(), old_phred_to_new_phred, var, new_var);
    new_vars.push_back(std::move(new_var));
  }

  return new_vars;
}


/*********************
 * OPERATOR OVERLOAD *
 *********************/
bool
Variant::operator==(Variant const & v) const
{
  return (abs_pos == v.abs_pos) &&
         (seqs.size() == v.seqs.size()) &&
         (seqs == v.seqs);
}


bool
Variant::operator!=(Variant const & b) const
{
  return !(*this == b);
}


bool
Variant::operator<(Variant const & b) const
{
  return abs_pos < b.abs_pos || (abs_pos == b.abs_pos && seqs < b.seqs);
}


//Variant &
//Variant::operator=(Variant const & var)
//{
//  abs_pos = var.abs_pos;
//  seqs = var.seqs;
//  calls = var.calls;
//  infos = var.infos;
//  return *this;
//}
//
//
//Variant &
//Variant::operator=(Variant && var) noexcept
//{
//  abs_pos = var.abs_pos;
//  seqs = std::move(var.seqs);
//  calls = std::move(var.calls);
//  infos = std::move(var.infos);
//  return *this;
//}


std::size_t
VariantHash::operator()(Variant const & v) const
{
  assert(v.seqs.size() == 2);
  std::size_t h1 = std::hash<uint32_t>()(v.abs_pos);
  std::size_t h2 = boost::hash_range(v.seqs[0].begin(), v.seqs[0].end());
  std::size_t h3 = 42 + boost::hash_range(v.seqs[1].begin(), v.seqs[1].end());
  return h1 ^ (h2 << 1) ^ (h3 + 0x9e3779b9);
}


} // namespace gyper
