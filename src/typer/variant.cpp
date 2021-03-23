#include <algorithm> // std::swap
#include <limits>
#include <string> // std::string
#include <sstream> // std::stringstream
#include <vector> // std::vector

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/log/trivial.hpp>
#include <boost/functional/hash.hpp> // boost::hash_range
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

#include <paw/align.hpp>

#include <graphtyper/graph/absolute_position.hpp> // gyper::AbsolutePosition
#include <graphtyper/graph/graph.hpp> // gyper::Graph
#include <graphtyper/graph/sv.hpp> // gyper::SVTYPE
#include <graphtyper/typer/variant.hpp> // gyper::Variant
#include <graphtyper/typer/variant_candidate.hpp> // gyper::VariantCandidate
#include <graphtyper/typer/logistic_constants.hpp> // gyper::LOGF_ constants
#include <graphtyper/utilities/graph_help_functions.hpp> // gyper::to_pair, gyper::to_index
#include <graphtyper/utilities/options.hpp> // gyper::Options
#include <graphtyper/utilities/sequence_operations.hpp> // gyper::remove_common_prefix, gyper::remove_common_suffix


namespace
{

void
update_per_allele_stats(std::size_t const num_seqs,
                        std::size_t new_num_seqs,
                        std::vector<uint16_t> const & old_phred_to_new_phred,
                        gyper::Variant const & var,
                        gyper::Variant & new_var)
{
  //assert(var.infos.count("SBF1") > 0);

  // Check if any infos to update
  if (var.stats.per_allele.size() == 0)
    return;

  gyper::VarStats const & old_stats = var.stats;//(num_seqs);
  //old_stats.read_stats(var.infos);

  new_var.stats = gyper::VarStats(new_num_seqs);
  gyper::VarStats & new_stats = new_var.stats;
  new_stats.clipped_reads = old_stats.clipped_reads;
  new_stats.mapq_squared = old_stats.mapq_squared;
  //new_stats.are_calls_scanned = false;

  for (uint32_t y = 0; y < num_seqs; ++y)
  {
    uint32_t new_y = old_phred_to_new_phred[y];

    assert(new_y < new_num_seqs);
    assert(new_y < new_stats.per_allele.size());
    assert(new_y < new_stats.read_strand.size());
    assert(y < old_stats.per_allele.size());
    assert(y < old_stats.read_strand.size());

    auto const & old_a = old_stats.per_allele[y];
    auto & new_a = new_stats.per_allele[new_y];
    auto const & old_rs = old_stats.read_strand[y];
    auto & new_rs = new_stats.read_strand[new_y];

    new_a.clipped_bp += old_a.clipped_bp;
    new_a.mapq_squared += old_a.mapq_squared;
    new_a.score_diff += old_a.score_diff;
    new_a.mismatches += old_a.mismatches;

    new_rs.r1_forward += old_rs.r1_forward;
    new_rs.r2_forward += old_rs.r2_forward;
    new_rs.r1_reverse += old_rs.r1_reverse;
    new_rs.r2_reverse += old_rs.r2_reverse;
  }

  //new_stats.add_stats(new_var.infos);
}


} // anon namesapce


namespace gyper
{

Variant::Variant(Variant const & var) noexcept
  : abs_pos(var.abs_pos)
  , seqs(var.seqs)
  , calls(var.calls)
  , stats(var.stats)
  , infos(var.infos)
  , suffix_id(var.suffix_id)
  , hap_id(var.hap_id)
  , type(var.type)
{}

Variant::Variant(Variant && var) noexcept
  : abs_pos(std::forward<uint32_t>(var.abs_pos))
  , seqs(std::forward<std::vector<std::vector<char> > >(var.seqs))
  , calls(std::forward<std::vector<SampleCall> >(var.calls))
  , stats(std::forward<VarStats>(var.stats))
  , infos(std::forward<std::map<std::string, std::string> >(var.infos))
  , suffix_id(std::forward<std::string>(var.suffix_id))
  , hap_id(var.hap_id)
  , type(var.type)
{}


Variant &
Variant::operator=(Variant const & o) noexcept
{
  abs_pos = o.abs_pos;
  seqs = o.seqs;
  calls = o.calls;
  stats = o.stats;
  infos = o.infos;
  suffix_id = o.suffix_id;
  hap_id = o.hap_id;
  type = o.type;
  return *this;
}


Variant &
Variant::operator=(Variant && o) noexcept
{
  abs_pos = o.abs_pos;
  seqs = std::move(o.seqs);
  calls = std::move(o.calls);
  stats = std::move(o.stats);
  infos = std::move(o.infos);
  suffix_id = std::move(o.suffix_id);
  hap_id = o.hap_id;
  type = o.type;
  return *this;
}


Variant::Variant(Genotype const & gt)
  : stats(0)
{
  seqs = graph.get_all_sequences_of_a_genotype(gt);
  stats = VarStats(seqs.size());
  abs_pos = graph.genomic_region.get_absolute_position(gt.id); // -1 cause we always fetch one position back as well
}


/*
Variant::Variant(std::vector<Genotype> const & gts, std::vector<uint16_t> const & hap_calls)
  : stats(0)
{
  assert(gts.size() > 0);

  for (long i = 0; i < static_cast<long>(hap_calls.size()); ++i)
    seqs.push_back(graph.get_sequence_of_a_haplotype_call(gts, hap_calls[i]));

  stats = VarStats(seqs.size());
  abs_pos = graph.genomic_region.get_absolute_position(gts[0].id); // -1 cause we always fetch one position back as well
}
*/


Variant::Variant(VariantCandidate const & var_candidate) noexcept
  : abs_pos(var_candidate.abs_pos)
  , seqs(var_candidate.seqs)
  , stats(0)
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

        long constexpr ERROR{4};
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
Variant::scan_calls()
{
  if (stats.seqdepth > 0 || stats.n_calls > 0)
    return;

  if (stats.per_allele.size() == 0)
  {
    stats.per_allele.resize(seqs.size());
    stats.read_strand.resize(seqs.size());
  }

  long const num_alts = seqs.size() - 1l;

  //BOOST_LOG_TRIVIAL(info) << __HERE__ << " " << calls.size();

  stats.n_calls += calls.size();

  for (auto const & sample_call : calls)
  {
    // Skip homozygous ref calls
    if (sample_call.phred.size() > 0 && sample_call.phred[0] > 0)
    {
      std::pair<uint16_t, uint16_t> gt_call = sample_call.get_gt_call();
      auto const & cov = sample_call.coverage;
      assert(gt_call.first < cov.size());
      assert(gt_call.second < cov.size());

      if (gt_call.first > 0)
      {
        assert(gt_call.first < stats.per_allele.size());
        auto & per_al = stats.per_allele[gt_call.first];
        long const depth = std::min(10l, static_cast<long>(cov[gt_call.first] + sample_call.ambiguous_depth));

        if (depth > 0)
        {
          per_al.qd_qual += std::min(25l * depth,
                                     static_cast<long>(sample_call.get_lowest_phred_not_with(gt_call.first)));

          per_al.qd_depth += depth;
        }
      }

      if (gt_call.first != gt_call.second)
      {
        assert(gt_call.second < stats.per_allele.size());
        auto & per_al = stats.per_allele[gt_call.second];
        assert(gt_call.second > 0);
        long const depth = std::min(10l, static_cast<long>(cov[gt_call.second] + sample_call.ambiguous_depth));

        if (depth > 0)
        {
          per_al.qd_qual +=
            std::min(25l * depth, static_cast<long>(sample_call.get_lowest_phred_not_with(gt_call.second)));

          per_al.qd_depth += depth;
        }
      }
    }

    // do stuff for all calls
    // MaxAltPP
    stats.n_max_alt_proper_pairs = std::max(stats.n_max_alt_proper_pairs, sample_call.alt_proper_pair_depth);

    uint32_t const total_depth = std::accumulate(sample_call.coverage.cbegin(),
                                                 sample_call.coverage.cend(),
                                                 0u);

    assert(static_cast<long>(sample_call.coverage.size()) == num_alts + 1);

    for (long c{0}; c < num_alts; ++c)
    {
      assert((c + 1) < static_cast<long>(stats.per_allele.size()));
      // Skip zero (the reference)
      auto & per_al = stats.per_allele[c + 1];
      per_al.maximum_alt_support = std::max(per_al.maximum_alt_support, sample_call.coverage[c + 1]);

      if (total_depth > 0)
      {
        double const current_ratio = static_cast<double>(sample_call.coverage[c + 1]) /
                                     static_cast<double>(total_depth);
        per_al.maximum_alt_support_ratio = std::max(per_al.maximum_alt_support_ratio, current_ratio);
      }
    }

    // Get the GT
    std::pair<uint16_t, uint16_t> call = sample_call.get_gt_call();
    assert(seqs.size() == sample_call.coverage.size());
    assert(call.first < seqs.size());
    assert(call.second < seqs.size());
    long const filter = sample_call.check_filter(sample_call.get_gq());

    {
      auto pl_non_zero = [](uint8_t const pl) -> bool
                         {
                           return pl != 0;
                         };

      if (std::find_if(sample_call.phred.begin(), sample_call.phred.end(), pl_non_zero) != sample_call.phred.end())
        ++stats.n_genotyped;
    }

    if (filter == 0)
      ++stats.n_passed_calls;

    // ABHet and ABHom
    {
      // Check if heterozygous
      if (call.first != call.second)
      {
        stats.het_allele_depth.first += sample_call.coverage[call.first];
        stats.het_allele_depth.second += sample_call.coverage[call.second];
      }
      else
      {
        stats.hom_allele_depth.first += sample_call.coverage[call.first];
        stats.hom_allele_depth.second += std::accumulate(sample_call.coverage.cbegin(),
                                                         sample_call.coverage.cend(),
                                                         0ull
                                                         ) - sample_call.coverage[call.first];
      }
    }

    {
      uint32_t const call_depth = sample_call.get_unique_depth();

      // Check if heterozygous
      if (call.first != call.second)
      {
        auto const & c1 = call.first;
        auto const & c2 = call.second;
        assert(c1 < static_cast<long>(stats.per_allele.size()));
        assert(c2 < static_cast<long>(stats.per_allele.size()));

        stats.per_allele[c1].het_multi_allele_depth.first += sample_call.coverage[c1];
        stats.per_allele[c1].het_multi_allele_depth.second += call_depth - sample_call.coverage[c1];

        stats.per_allele[c2].het_multi_allele_depth.first += sample_call.coverage[c2];
        stats.per_allele[c2].het_multi_allele_depth.second += call_depth - sample_call.coverage[c2];
      }
      else
      {
        auto const & c = call.first;
        assert(static_cast<long>(c) < static_cast<long>(stats.per_allele.size()));
        stats.per_allele[c].hom_multi_allele_depth.first += sample_call.coverage[c];
        stats.per_allele[c].hom_multi_allele_depth.second += call_depth - sample_call.coverage[c];
      }
    }

    // SeqDepth
    if (sample_call.coverage.size() > 0)
    {
      stats.seqdepth += sample_call.get_depth();
      assert(sample_call.coverage.size() == stats.per_allele.size());

      for (long c{1}; c < static_cast<long>(sample_call.coverage.size()); ++c)
        stats.per_allele[c].total_depth += sample_call.coverage[c];
    }

    // AC (AF can be derived from AC and AN)
    {
      assert(call.first < static_cast<long>(stats.per_allele.size()));
      assert(call.second < static_cast<long>(stats.per_allele.size()));

      ++stats.per_allele[call.first].ac;
      ++stats.per_allele[call.second].ac;

      if (filter == 0) // PASS
      {
        ++stats.per_allele[call.first].pass_ac;
        ++stats.per_allele[call.second].pass_ac;
      }

      if (call.first == 0u && call.second == 0u)
        ++stats.n_ref_ref;
      else if (call.first == 0u || call.second == 0u)
        ++stats.n_ref_alt;
      else
        ++stats.n_alt_alt;
    }
  }
}


std::vector<int8_t>
Variant::generate_infos()
{
  assert(seqs.size() >= 2);

  long const num_seqs = seqs.size();
  long const num_alts = seqs.size() - 1l;

  // Write stats from VarStat object
  bool const is_stats = stats.per_allele.size() != 0;

  if (is_stats && stats.per_allele.size() != seqs.size())
  {
    BOOST_LOG_TRIVIAL(error) << __HERE__ << " stats.per_allele.size() != seqs.size(), "
                             << stats.per_allele.size() << " != " << seqs.size();
    std::exit(1);
  }

  if (is_stats)
  {
    scan_calls();
    stats.write_stats(infos);
  }
  else
  {
    stats.per_allele.resize(num_seqs);
    stats.read_strand.resize(num_seqs);
    scan_calls();
  }

  std::vector<int8_t> is_good_alt(seqs.size() - 1, 1);

  infos["RefLen"] = std::to_string(seqs[0].size());

  // Make sure END is greater or equal to POS
  {
    auto find_it = infos.find("END");

    if (find_it != infos.end())
    {
      auto contig_pos = absolute_pos.get_contig_position(abs_pos, gyper::graph.contigs);
      long end = std::strtol(find_it->second.c_str(), nullptr, 10);

      if (end < contig_pos.second)
        end = contig_pos.second;

      infos["END"] = std::to_string(end);
    }
  }


  // Calculate AC, AN, SeqDepth, MaxVS, ABHom, ABHet, ABHomMulti, ABHetMulti
  //std::vector<std::pair<uint32_t, uint32_t> > het_multi_allele_depth(seqs.size(), {0u, 0u});
  //std::vector<std::pair<uint32_t, uint32_t> > hom_multi_allele_depth(seqs.size(), {0u, 0u});

  // Write MaxAAS
  {
    //auto const & per_al =
    assert(stats.per_allele.size() > 0);
    assert(static_cast<long>(stats.per_allele.size()) == num_seqs);
    std::stringstream max_vs_ss;
    // skip 0 (reference)
    max_vs_ss << static_cast<uint16_t>(stats.per_allele[1].maximum_alt_support);

    for (long e{2}; e < num_seqs; ++e)
      max_vs_ss << ',' << static_cast<uint16_t>(stats.per_allele[e].maximum_alt_support);

    infos["MaxAAS"] = max_vs_ss.str();
  }

  // Write MaxAASR
  {
    assert(stats.per_allele.size() > 0);
    assert(static_cast<long>(stats.per_allele.size()) == num_seqs);
    std::stringstream ss_max_asr;
    ss_max_asr.precision(4);
    ss_max_asr << stats.per_allele[1].maximum_alt_support_ratio;

    for (long e{2}; e < num_seqs; ++e)
      ss_max_asr << ',' << stats.per_allele[e].maximum_alt_support_ratio;

    infos["MaxAASR"] = ss_max_asr.str();
  }

  // Write MaxAltPP
  if (is_sv())
  {
    infos["MaxAltPP"] = std::to_string(static_cast<uint16_t>(stats.n_max_alt_proper_pairs));
  }

  // Write AC
  {
    assert(stats.per_allele.size() > 0);
    assert(stats.per_allele.size() == seqs.size());
    std::ostringstream ac_ss;
    ac_ss << stats.per_allele[1].ac;

    for (uint32_t e{2}; e < stats.per_allele.size(); ++e)
      ac_ss << ',' << stats.per_allele[e].ac;

    infos["AC"] = ac_ss.str();
  }

  // Write AN
  infos["AN"] = std::to_string(2 * stats.n_genotyped);

  // Write AF
  {
    std::ostringstream af_ss;
    af_ss.precision(4);

    if (stats.n_genotyped > 0)
      af_ss <<
        static_cast<double>(static_cast<double>(stats.per_allele[1].ac) / static_cast<double>(2 * stats.n_genotyped));
    else
      af_ss << "0.0";

    for (long e{2}; e < static_cast<long>(stats.per_allele.size()); ++e)
    {
      af_ss << ',';

      if (stats.n_genotyped > 0)
        af_ss <<
          static_cast<double>(static_cast<double>(stats.per_allele[e].ac) / static_cast<double>(2 * stats.n_genotyped));
      else
        af_ss << "0.0";
    }

    infos["AF"] = af_ss.str();
  }

  // Write PASS_AC
  {
    assert(stats.per_allele.size() == seqs.size());
    std::ostringstream ac_ss;
    ac_ss << stats.per_allele[1].pass_ac;

    for (long e{2}; e < static_cast<long>(stats.per_allele.size()); ++e)
      ac_ss << ',' << stats.per_allele[e].pass_ac;

    infos["PASS_AC"] = ac_ss.str();
  }

  // Write PASS_AN
  infos["PASS_AN"] = std::to_string(2 * stats.n_passed_calls);

  double info_pass_ratio{0.0};

  // Write PASS_ratio
  if (stats.n_genotyped > 0)
  {
    info_pass_ratio = static_cast<double>(stats.n_passed_calls) / static_cast<double>(stats.n_genotyped);
    std::stringstream ss_pr;
    ss_pr.precision(4);
    ss_pr << info_pass_ratio;
    infos["PASS_ratio"] = ss_pr.str();
  }

  // Write Num REF/REF, REF/ALT, ALT/ALT, homozygous and heterozygous calls
  {
    std::ostringstream n_gt;
    n_gt << stats.n_ref_ref << "," << stats.n_ref_alt << "," << stats.n_alt_alt;
    infos["NHomRef"] = std::to_string(stats.n_ref_ref);
    infos["NHet"] = std::to_string(stats.n_ref_alt);
    infos["NHomAlt"] = std::to_string(stats.n_alt_alt);
  }

  // Write SeqDepth
  {
    infos["SeqDepth"] = std::to_string(stats.seqdepth);
  }

  double info_ab_het{0.5};

  // Write ABHet
  {
    std::stringstream ss_abhet;
    ss_abhet.precision(4);
    uint32_t const total_het_depth = stats.het_allele_depth.first + stats.het_allele_depth.second;

    if (total_het_depth > 0)
    {
      info_ab_het = static_cast<double>(stats.het_allele_depth.second) / static_cast<double>(total_het_depth);
      ss_abhet << info_ab_het;
    }
    else
    {
      ss_abhet << "-1";
    }

    infos["ABHet"] = ss_abhet.str();
  }

  double info_abhom{0.985};

  // Write ABHom
  {
    std::stringstream ss_abhom;
    ss_abhom.precision(4);
    uint32_t const total_hom_depth = stats.hom_allele_depth.first + stats.hom_allele_depth.second;

    if (total_hom_depth > 0)
    {
      info_abhom = static_cast<double>(stats.hom_allele_depth.first) / static_cast<double>(total_hom_depth);
      ss_abhom << info_abhom;
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

  double info_sbalt{0.0};

  // Write Alternative Strand Bias (SBAlt)
  {
    uint32_t total_f = get_accumulated_alt_strand_bias(infos, "SBF");
    uint32_t total_r = get_accumulated_alt_strand_bias(infos, "SBR");

    std::stringstream ss_sb;
    ss_sb.precision(4);

    if (total_f + total_r == 0)
    {
      ss_sb << "-1";
    }
    else
    {
      info_sbalt = static_cast<double>(total_f) / static_cast<double>(total_f + total_r);
      ss_sb << info_sbalt;
    }

    infos["SBAlt"] = ss_sb.str();
  }

  // Write ABHetMulti
  {
    assert(stats.per_allele.size() == seqs.size());
    std::stringstream ss_abhet;
    ss_abhet.precision(4);

    for (unsigned i{0}; i < stats.per_allele.size(); ++i)
    {
      if (i > 0)
        ss_abhet << ",";

      auto const & het_depth_pair = stats.per_allele[i].het_multi_allele_depth;
      uint32_t const total_het_depth = het_depth_pair.first + het_depth_pair.second;

      if (total_het_depth > 0)
      {
        ss_abhet << (static_cast<double>(het_depth_pair.second) /
                     static_cast<double>(total_het_depth));
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
    assert(stats.per_allele.size() == seqs.size());
    std::stringstream ss_abhom;
    ss_abhom.precision(4);

    for (unsigned i{0}; i < stats.per_allele.size(); ++i)
    {
      if (i > 0)
        ss_abhom << ",";

      auto const & hom_depth_pair = stats.per_allele[i].hom_multi_allele_depth;
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

  // Write VarType
  {
    infos["VarType"] = this->determine_variant_type();
  }

  double info_qd{0.0};

  // Calculate QD
  {
    std::stringstream ss;
    ss.precision(4);
    info_qd = get_qual_by_depth();
    ss << info_qd;
    infos["QD"] = ss.str();
  }

  std::vector<double> qd_alt = get_qual_by_depth_per_alt_allele();

  // Calculate QDalt
  {
    std::stringstream ss;
    ss.precision(4);
    ss << qd_alt[0];

    for (long q{1}; q < static_cast<long>(qd_alt.size()); ++q)
    {
      ss.precision(4);
      ss << ',' << qd_alt[q];
    }

    infos["QDalt"] = ss.str();
  }

  long info_mq{60};

  // Calculate MQ
  {
    //auto find_mapq_squared_it = infos.find("MQsquared");
    //
    //if (find_mapq_squared_it != infos.end())
    {
      //double mapq_squared = std::stoull(find_mapq_squared_it->second);
      double mapq_squared = stats.mapq_squared;

      if (stats.seqdepth > 0)
      {
        double mapq = std::sqrt(mapq_squared / static_cast<double>(stats.seqdepth));
        info_mq = static_cast<long>(mapq + 0.5);
        infos["MQ"] = std::to_string(info_mq);
      }
      else
      {
        infos["MQ"] = "0";
      }
    }
  }

  if (Options::const_instance()->is_segment_calling)
  {
    assert(stats.per_allele.size() == seqs.size());

    for (long a{1}; a < static_cast<long>(stats.per_allele.size()); ++a)
    {
      is_good_alt[a - 1] = stats.per_allele[a].ac > 0;
    }

    // erase stats that don't make sense in segment calling
    infos.erase("ABHet");
    infos.erase("ABHom");
    infos.erase("ABHetMulti");
    infos.erase("ABHomMulti");
    infos.erase("MaxAAS");
    infos.erase("MaxAASR");
    infos.erase("QD");
    infos.erase("QDalt");
    infos.erase("SB");
    infos.erase("SBAlt");
    infos.erase("SeqDepth");

    return is_good_alt;
  }


  assert(stats.per_allele.size() == seqs.size());

  // Calculate SDalt, MMalt, CRalt and MQalt
  if (is_stats)
  {
    if (stats.per_allele.size() == seqs.size())
    {
      assert(stats.read_strand.size() == seqs.size());
      assert(stats.per_allele.size() == stats.read_strand.size());

      std::ostringstream sd_ss;
      std::ostringstream mm_ss;
      std::ostringstream cr_ss;
      std::ostringstream mq_ss;

      for (long s{1}; s < static_cast<long>(stats.per_allele.size()); ++s)
      {
        auto const & per_al = stats.per_allele[s];
        uint64_t const alt_depth{per_al.total_depth};

        if (s > 1)
        {
          sd_ss << ',';
          mm_ss << ',';
          cr_ss << ',';
          mq_ss << ',';
        }

        if (alt_depth > 0)
        {
          sd_ss << (static_cast<double>(per_al.score_diff) / static_cast<double>(alt_depth));
          mm_ss << (static_cast<double>(per_al.mismatches) / static_cast<double>(alt_depth) / 10.0);
          cr_ss << (static_cast<double>(per_al.clipped_bp) / static_cast<double>(alt_depth) / 10.0);
          mq_ss << static_cast<long>(std::sqrt(static_cast<double>(per_al.mapq_squared) /
                                               static_cast<double>(alt_depth)) + 0.5);
        }
        else
        {
          sd_ss << "0.0";
          mm_ss << "0.0";
          cr_ss << "0.0";
          mq_ss << "0";
        }
      }

      infos["SDalt"] = sd_ss.str();
      infos["MMalt"] = mm_ss.str();
      infos["CRalt"] = cr_ss.str();
      infos["MQalt"] = mq_ss.str();
    }

    std::vector<uint64_t> sb_alt(seqs.size() - 1, 0);
    assert(stats.read_strand.size() == seqs.size());

    for (long s{0}; s < static_cast<long>(sb_alt.size()); ++s)
      sb_alt[s] = stats.read_strand[s + 1].get_reverse_count();

    std::vector<double> aa_score;

    // Alternative allele score (AAScore)
    if (graph.is_sv_graph)
    {
      aa_score = std::vector<double>(num_alts, 1.0);
    }
    else
    {
      aa_score.resize(num_alts);

      for (long s{0}; s < num_alts; ++s)
      {
        assert((s + 1) < static_cast<long>(stats.per_allele.size()));
        auto const & per_al = stats.per_allele[s + 1];

        assert(s < static_cast<long>(qd_alt.size()));
        //assert(s < static_cast<long>(maximum_variant_support.size()));
        double const qd = qd_alt[s];

        if (per_al.total_depth > 0 &&
            qd > 0.1 &&
            per_al.maximum_alt_support >= 3 &&
            per_al.maximum_alt_support_ratio >= 0.175)
        {
          double const alt_seq_depth = per_al.total_depth;
          double const _sb = 2.0 * ((static_cast<double>(sb_alt[s]) / alt_seq_depth) - 0.5);
          double const sb = _sb >= 0.0 ? _sb : -_sb;
          assert(sb >= 0.0);
          double const mm = static_cast<double>(per_al.mismatches) / alt_seq_depth / 10.0;
          long const sd = static_cast<double>(per_al.score_diff) / alt_seq_depth + 0.5;
          double const cr = static_cast<double>(per_al.clipped_bp) / alt_seq_depth / 10.0;
          long const mq = std::sqrt(static_cast<double>(per_al.mapq_squared) / alt_seq_depth) + 0.5;

          //BOOST_LOG_TRIVIAL(info) << __HERE__ << " " << sb << "," << alt_seq_depth << "," << mm
          //                        << "," << sd << "," << cr << "," << mq;

          double score = get_aa_score(info_abhom, sb, mm, sd, qd, cr, mq);

          if (mm > 1.5)
          {
            double m = 1.0 - ((mm - 1.5) / 20.0);
            m = m <= 0.5 ? 0.5 : m;
            score *= m;
          }

          if ((cr + mm) > 2.5)
          {
            double m = 1.0 - ((cr + mm - 2.5) / 40.0);
            m = m <= 0.5 ? 0.5 : m;
            score *= m;
          }

          aa_score[s] = score;
        }
        else
        {
          aa_score[s] = 0.0;
        }
      }

      std::ostringstream ss;
      ss.precision(4);
      ss << aa_score[0];

      for (long s{1}; s < static_cast<long>(aa_score.size()); ++s)
      {
        ss.precision(4);
        ss << "," << aa_score[s];
      }

      infos["AAScore"] = ss.str();

      // Logistic regression quality filter (LOGF).
      {
        long const info_cr = infos.count("CR") ? std::stol(infos.at("CR")) : 0;
        long const ab_het_bin = static_cast<long>(info_ab_het * 10.0 + 0.00001);
        long const sbalt_bin = static_cast<long>(info_sbalt * 10.0 + 0.00001);
        double const cr_by_seqdepth = static_cast<double>(info_cr) / static_cast<double>(stats.seqdepth);
        double const gt_yield = static_cast<double>(stats.n_genotyped) / static_cast<double>(stats.n_calls);

        // variables:
        //  crBySeqdepth=info_cr/seqdepth
        //  info_qd
        //  info_abhom
        //  info_mq
        //  info_pass_ratio
        //  gtYield=genotyped/an
        assert(ab_het_bin >= 0l);
        assert(ab_het_bin <= 10l);
        assert(sbalt_bin >= 0l);
        assert(sbalt_bin <= 10l);

        double const logf = get_logf(info_abhom, cr_by_seqdepth, info_mq,
                                     info_pass_ratio, gt_yield, info_qd, ab_het_bin, sbalt_bin);

        std::ostringstream ss;
        ss.precision(4);
        ss << logf;
        infos["LOGF"] = ss.str();
      }
    }

    // check which alts are good
    for (long a{0}; a < static_cast<long>(seqs.size()) - 1l; ++a)
    {
      //assert(maximum_variant_support.size() == is_good_alt.size());
      //assert(maximum_alternative_support_ratio.size() == is_good_alt.size());
      assert(aa_score.size() == is_good_alt.size());
      assert(qd_alt.size() == is_good_alt.size());
      assert(stats.per_allele.size() == is_good_alt.size() + 1);
      auto const & per_al = stats.per_allele[a + 1];

      if (per_al.total_depth == 0)
      {
        is_good_alt[a] = 0;
        continue;
      }

      assert(a < static_cast<long>(qd_alt.size()));
      //assert(s < static_cast<long>(maximum_variant_support.size()));
      double const qd = qd_alt[a];

      is_good_alt[a] = aa_score[a] >= 0.05 &&
                       qd > 1.0 &&
                       per_al.maximum_alt_support >= 3 &&
                       per_al.maximum_alt_support_ratio >= 0.175 &&
                       (seqs.size() < 71 || (qd > 1.25 && per_al.maximum_alt_support_ratio >= 0.2)) &&
                       (seqs.size() < 131 || (qd > 1.5 && per_al.maximum_alt_support_ratio >= 0.225));

#ifndef NDEBUG
      if (is_good_alt[a] == 0)
      {
        BOOST_LOG_TRIVIAL(debug) << __HERE__ << " In variant=" << to_string(true) // skip calls
                                 << " bad alt="
                                 << std::string(seqs[a + 1].begin(), seqs[a + 1].end())
                                 << " MaxAAS=" << per_al.maximum_alt_support
                                 << " MaxAASR=" << per_al.maximum_alt_support_ratio
                                 << " AAScore=" << aa_score[a]
                                 << " ABHom=" << info_abhom
                                 << " QDAlt=" << qd_alt[a];
      }
#endif // NDEBUG
    }
  }

  return is_good_alt;
}


bool
Variant::is_sv() const
{
  assert(seqs.size() > 0);

  if (Options::const_instance()->is_segment_calling)
    return false;

  // If the variant has an SV breakpoint, skip trimming
  for (long s{1}; s < static_cast<long>(seqs.size()); ++s)
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
  //uint32_t abs_pos_copy = abs_pos - graph.genomic_region.get_absolute_begin_position() + 1;
  //uint32_t abs_pos_copy = abs_pos - graph.genomic_region.get_absolute_begin_position() + 1; //graph.genomic_region.abs_pos;
  uint32_t contig_pos = absolute_pos.get_contig_position(abs_pos, graph.contigs).second;
  uint32_t contig_pos_cp = contig_pos;
  uint32_t new_contig_pos = contig_pos - 1;
  std::vector<char> first_base = graph.get_generated_reference_genome(new_contig_pos, contig_pos);

  if (first_base.size() != 1 || contig_pos_cp != contig_pos || new_contig_pos != contig_pos - 1)
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


long
Variant::normalize()
{
  if (seqs.size() < 2)
  {
    BOOST_LOG_TRIVIAL(warning) << __HERE__ << " Tried to normalize variant with " << seqs.size() << " alleles";
    return 0;
  }

  auto const & ref = seqs[0];

  // Check if some sequence is of length 0 or prefix is not matching
  for (long i = 0; i < static_cast<long>(seqs.size()); ++i)
  {
    auto const & seq = seqs[i];

    if (seq.size() == 0)
    {
      BOOST_LOG_TRIVIAL(warning) << __HERE__ << " Tried to normalize variant with size 0: " << to_string();
      return 0;
    }

    if (seq[0] != ref[0])
      return 0;

    // If a sequence is the same as the ref then we will loop forever
    if (i > 0 && seq == ref)
      return 0;
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

  long distance{0};

  while (all_last_bases_match())
  {
    bool const success_adding_base = add_base_in_front();

    if (not success_adding_base)
      break;

    ++distance;
    // Remove the last base
    remove_common_suffix(seqs);
  }

  gyper::remove_common_prefix(abs_pos, seqs, false); // Don't keep one match
  return distance;
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
Variant::to_string(bool is_skipping_calls) const
{
  std::stringstream os;
  auto contig_pos = absolute_pos.get_contig_position(this->abs_pos, gyper::graph.contigs);
  os << contig_pos.first << ":" << contig_pos.second;

  if (this->seqs.size() > 0)
    os << " " << std::string(this->seqs[0].begin(), this->seqs[0].end());

  if (this->seqs.size() > 1)
    os << " " << std::string(this->seqs[1].begin(), this->seqs[1].end());
  else
    os << " N";

  for (long i{2}; i < static_cast<long>(this->seqs.size()); ++i)
    os << "," << std::string(this->seqs[i].begin(), this->seqs[i].end());

  if (!is_skipping_calls)
  {
    for (auto const & sample : this->calls)
    {
      os << " ";

      if (sample.phred.size() > 0)
        os << static_cast<uint16_t>(sample.phred[0]);

      for (unsigned p = 1; p < sample.phred.size(); ++p)
        os << "," << static_cast<uint16_t>(sample.phred[p]);
    }
  }

  return os.str();
}


std::string
Variant::determine_variant_type() const
{
  assert(seqs.size() >= 2);

  if (Options::const_instance()->is_segment_calling && seqs[0].size() > 0 && seqs[0][0] == '<')
  {
    return "H";
  }

  // Count the number of alleles which are not of size 1
  std::size_t num_non_ones{0};

  // Check if the variant is an SV
  gyper::SVTYPE sv_type{NOT_SV};

  for (auto it = seqs.begin(); it != seqs.end(); ++it)
  {
#ifndef NDEBUG
    if (it->size() == 0)
    {
      BOOST_LOG_TRIVIAL(error) << __HERE__ << " " << this->to_string();
    }

    assert(it->size() != 0);
#endif
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
    if (sample_call.phred.size() > 0)
      variant_qual += sample_call.phred[0];
  }

  return variant_qual;
}


double
Variant::get_qual_by_depth() const
{
  long total_qual{0};
  long total_depth{0};

  for (auto const & sample_call : calls)
  {
    // Skip homozygous ref calls
    if (sample_call.phred.size() > 0 && sample_call.phred[0] > 0)
    {
      long const depth = std::min(10l, static_cast<long>(sample_call.get_alt_depth()));

      if (depth > 0)
      {
        total_qual += std::min(25l * depth, static_cast<long>(sample_call.phred[0]));
        total_depth += depth;
      }
    }
  }

  if (total_depth == 0)
    return 0.0;
  else
    return static_cast<double>(total_qual) / static_cast<double>(total_depth);
}


std::vector<double>
Variant::get_qual_by_depth_per_alt_allele() const
{
  assert(seqs.size() > 0);
  long const num_alts = seqs.size() - 1;
  std::vector<double> qd(num_alts, 0.0);

  for (long s{0}; s < num_alts; ++s)
  {
    auto const & per_al = stats.per_allele[s + 1];

    if (per_al.qd_depth > 0)
      qd[s] = static_cast<double>(per_al.qd_qual) / static_cast<double>(per_al.qd_depth);
  }

  return qd;
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

    update_per_allele_stats(var.seqs.size(), new_var.seqs.size(), old_phred_to_new_phred, var, new_var);
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
  //assert(var.stats.per_allele.size() == var.seqs.size());

  if (Options::const_instance()->no_decompose ||
      (var.seqs.size() == 2 &&
       std::any_of(var.seqs[1].begin(), var.seqs[1].end(), [](char const c){
      return c == '<' || c == '[' || c == ']';
    })))
  {
    //var.stats.are_calls_scanned = true;
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
              std::back_inserter(broken_down_vars));
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

  for (long i{1}; i < static_cast<long>(reference.size()); ++i)
  {
    assert(new_var.seqs.size() >= 2);

    if (reference[i] == '-')
      ++ref_gaps;
    else
      new_var.seqs[0].push_back(reference[i]);

    bool all_match{true};

    for (int a{1}; a < static_cast<long>(new_var.seqs.size()); ++a)
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


SampleCall
bin_phred(Variant const & new_var,
          Variant const & old_var,
          SampleCall const & old_call,
          std::vector<long> const & old2new)
{
  //BOOST_LOG_TRIVIAL(info) << __HERE__ << " old_var.seqs.size()=" << old_var.seqs.size()
  //                        << " new_var.seqs.size()=" << new_var.seqs.size()
  //                        << " old2new.size()=" << old2new.size();
  SampleCall new_call;

  long const new_cnum = new_var.seqs.size();
  long const old_cnum = old_var.seqs.size();

  assert(new_cnum > 1);
  assert(old_cnum > 1);

  new_call.coverage.resize(new_cnum);
  new_call.phred.resize((new_cnum * (new_cnum + 1)) / 2, 255u);

  for (long y{0}; y < old_cnum; ++y)
  {
    assert(y < static_cast<long>(old2new.size()));
    long const new_y = old2new[y];

    for (long x{0}; x <= y; ++x)
    {
      assert(x < static_cast<long>(old2new.size()));
      long const old_i = to_index(x, y);

      assert(old_i < static_cast<long>(old_call.phred.size()));
      uint8_t const old_phred = old_call.phred[old_i];

      long const new_x = old2new[x];
      long const new_i = to_index_safe(new_x, new_y);

      assert(new_i < static_cast<long>(new_call.phred.size()));
      uint8_t const new_phred = new_call.phred[new_i];

      if (old_phred < new_phred)
        new_call.phred[new_i] = old_phred;
    }
  }

  //std::cerr << __HERE__ << " ";
  //
  //for (auto p : new_call.phred)
  //  std::cerr << static_cast<long>(p) << ",";
  //
  //std::cerr << std::endl;

  return new_call;
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
    update_per_allele_stats(old_var.seqs.size(),
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
  //assert(var.infos.count("SBF1") == 1);

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

    //if (var.infos.count("SBF1") == 1)
    //assert(var.infos.count("SBF1") == 1);
    update_per_allele_stats(seqs.size(), new_var.seqs.size(), old_phred_to_new_phred, var, new_var);

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

  long constexpr extra_bases_after{OPTIMAL_EXTRA};

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

#ifndef NDEBUG
    for (auto const & seq : new_var.seqs)
    {
      if (seq.size() == 0)
      {
        BOOST_LOG_TRIVIAL(warning) << __HERE__ << " Variant with size 0: " << new_var.to_string();
      }
    }
#endif // NDEBUG

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

    update_per_allele_stats(seqs.size(), new_var.seqs.size(), old_phred_to_new_phred, var, new_var);
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
  return abs_pos < b.abs_pos ||
         (abs_pos == b.abs_pos && type < b.type) ||
         ((abs_pos == b.abs_pos && type == b.type && seqs < b.seqs));
}


template <typename Archive>
void
Variant::serialize(Archive & ar, unsigned const int /*version*/)
{
  ar & abs_pos;
  ar & seqs;
  ar & calls;
  ar & stats;
  ar & infos;
  ar & suffix_id;
  ar & hap_id;
  ar & type;
}


template void Variant::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive &,
                                                                  const unsigned int);
template void Variant::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive &,
                                                                  const unsigned int);


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
