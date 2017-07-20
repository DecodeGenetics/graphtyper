#include <string> // std::string
#include <sstream> // std::stringstream
#include <vector> // std::vector

#include <boost/log/trivial.hpp>
#include <boost/functional/hash.hpp> // boost::hash_range

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/graph_msa.h>
#include <seqan/score.h>
#include <seqan/stream.h>

#include <graphtyper/graph/absolute_position.hpp> // gyper::AbsolutePosition
#include <graphtyper/graph/graph.hpp> // gyper::Graph
#include <graphtyper/typer/variant.hpp> // gyper::Variant
#include <graphtyper/typer/variant_candidate.hpp> // gyper::VariantCandidate
#include <graphtyper/utilities/graph_help_functions.hpp> // gyper::to_pair, gyper::to_index
#include <graphtyper/utilities/sequence_operations.hpp> // gyper::remove_common_prefix, gyper::remove_common_suffix


namespace
{

void
update_strand_bias(std::size_t const num_seqs, std::size_t new_num_seqs, std::vector<uint16_t> const & old_phred_to_new_phred, gyper::Variant const & var, gyper::Variant & new_var)
{
  using gyper::Variant;
  using gyper::get_strand_bias;
  using gyper::join_strand_bias;

  std::vector<uint32_t> originally_cropped = get_strand_bias(var.infos, "CRAligner");

  std::vector<uint32_t> mq_per_allele = get_strand_bias(var.infos, "MQperAllele");
  std::vector<uint64_t> seq_depths = var.get_seq_depth_of_all_alleles();

  std::vector<uint32_t> sbf = get_strand_bias(var.infos, "SBF");
  std::vector<uint32_t> sbr = get_strand_bias(var.infos, "SBR");

  std::vector<uint32_t> sbf1 = get_strand_bias(var.infos, "SBF1");
  std::vector<uint32_t> sbf2 = get_strand_bias(var.infos, "SBF2");
  std::vector<uint32_t> sbr1 = get_strand_bias(var.infos, "SBR1");
  std::vector<uint32_t> sbr2 = get_strand_bias(var.infos, "SBR2");

  std::vector<uint32_t> ra_count = get_strand_bias(var.infos, "RACount");
  std::vector<uint32_t> ra_dist = get_strand_bias(var.infos, "RADist");

  if (sbf.size() > 0 && sbr.size() > 0)
  {
    std::vector<uint32_t> new_originally_cropped(new_num_seqs, 0u);

    std::vector<uint64_t> new_mq_rooted(new_num_seqs, 0u);
    std::vector<uint64_t> new_seq_depths(new_num_seqs, 0u);

    std::vector<uint32_t> new_sbf(new_num_seqs, 0u);
    std::vector<uint32_t> new_sbr(new_num_seqs, 0u);

    std::vector<uint32_t> new_sbf1(new_num_seqs, 0u);
    std::vector<uint32_t> new_sbf2(new_num_seqs, 0u);
    std::vector<uint32_t> new_sbr1(new_num_seqs, 0u);
    std::vector<uint32_t> new_sbr2(new_num_seqs, 0u);

    std::vector<uint32_t> new_ra_count(new_num_seqs, 0u);
    std::vector<uint32_t> new_ra_dist(new_num_seqs, 0u);

    for (uint32_t y = 0; y < num_seqs; ++y)
    {
      uint32_t new_y = old_phred_to_new_phred[y];

      assert(new_y < new_num_seqs);
      assert(y < mq_per_allele.size());
      assert(y < seq_depths.size());
      assert(y < sbf.size());
      assert(y < sbr.size());
      assert(y < sbf1.size());
      assert(y < sbf2.size());
      assert(y < sbr1.size());
      assert(y < sbr2.size());

      new_originally_cropped[new_y] += originally_cropped[y];

      new_mq_rooted[new_y] += mq_per_allele[y] * mq_per_allele[y] * seq_depths[y];
      new_seq_depths[new_y] += seq_depths[y];

      new_sbf[new_y] += sbf[y];
      new_sbr[new_y] += sbr[y];

      new_sbf1[new_y] += sbf1[y];
      new_sbf2[new_y] += sbf2[y];
      new_sbr1[new_y] += sbr1[y];
      new_sbr2[new_y] += sbr2[y];

      new_ra_count[new_y] += ra_count[y];
      new_ra_dist[new_y] += ra_dist[y];
    }

    {
      std::ostringstream ss;
      std::size_t m = 0;

      if (new_seq_depths[m] > 0)
        ss << static_cast<uint16_t>(sqrt(static_cast<double>(new_mq_rooted[m]) / static_cast<double>(new_seq_depths[m])));
      else
        ss << 255u;

      for (m = 1; m < new_num_seqs; ++m)
      {
        if (new_seq_depths[m] > 0)
          ss << ',' << static_cast<uint16_t>(sqrt(static_cast<double>(new_mq_rooted[m]) / static_cast<double>(new_seq_depths[m])));
        else
          ss << ',' << 255u;
      }

      new_var.infos["MQperAllele"] = ss.str();
    }

    new_var.infos["CRAligner"] = join_strand_bias(new_originally_cropped);

    new_var.infos["SBF"] = join_strand_bias(new_sbf);
    new_var.infos["SBR"] = join_strand_bias(new_sbr);

    new_var.infos["SBF1"] = join_strand_bias(new_sbf1);
    new_var.infos["SBF2"] = join_strand_bias(new_sbf2);
    new_var.infos["SBR1"] = join_strand_bias(new_sbr1);
    new_var.infos["SBR2"] = join_strand_bias(new_sbr2);

    new_var.infos["RACount"] = join_strand_bias(new_ra_count);
    new_var.infos["RADist"] = join_strand_bias(new_ra_dist);
  }
}


void
find_variant_sequences(gyper::Variant & new_var, gyper::Variant const & old_var)
{
  using gyper::get_strand_bias;
  using gyper::join_strand_bias;
  using gyper::to_index;

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
      old_phred_to_new_phred.push_back(new_seqs.size());
      new_seqs.push_back(seqs[a]);
    }
    else
    {
      old_phred_to_new_phred.push_back(std::distance(new_seqs.begin(), find_it));
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
    gyper::SampleCall new_call;

    // phred
    {
      assert(old_call.phred.size() == (old_var.seqs.size() * (old_var.seqs.size() + 1) / 2));
      new_call.phred = std::vector<uint8_t>(new_seqs.size() * (new_seqs.size() + 1) / 2, 255u);

      for (uint16_t y = 0; y < old_var.seqs.size(); ++y)
      {
        for (uint16_t x = 0; x <= y; ++x)
        {
          uint16_t const index = to_index(x, y);
          uint16_t new_x = old_phred_to_new_phred[x];
          uint16_t new_y = old_phred_to_new_phred[y];

          if (old_var.phase.size() > 0 && old_call.phred[index] == 0)
          {
            assert (i < old_var.phase.size());

            if (new_x > new_y)
              new_var.phase[i] = !old_var.phase[i];
          }

          if (new_x > new_y)
            std::swap(new_x, new_y);

          uint16_t const new_index = to_index(new_x, new_y);
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
        if (static_cast<uint32_t>(new_call.coverage[new_y]) + static_cast<uint32_t>(old_call.coverage[y]) < 0xFFFFul)
          new_call.coverage[new_y] += old_call.coverage[y];
        else
          new_call.coverage[new_y] = 0xFFFFul;
      }
    }

    // multi-depth
    {
      new_call.ambiguous_depth = old_call.ambiguous_depth;
    }

    // Update strand bias
    update_strand_bias(old_var.seqs.size(), new_seqs.size(), old_phred_to_new_phred, old_var, new_var);
    new_var.calls.push_back(std::move(new_call));
  }

  new_var.seqs = std::move(new_seqs);
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
  , phase(var.phase)
{}

Variant::Variant(Variant && var) noexcept
  : abs_pos(std::forward<uint32_t>(var.abs_pos))
  , seqs(std::forward<std::vector<std::vector<char> > >(var.seqs))
  , calls(std::forward<std::vector<SampleCall> >(var.calls))
  , infos(std::forward<std::map<std::string, std::string> >(var.infos))
  , phase(std::forward<std::vector<uint8_t> >(var.phase))
{}


Variant::Variant(Genotype const & gt)
{
  seqs = graph.get_all_sequences_of_a_genotype(gt);
  abs_pos = gt.id - 1; // -1 cause we always fetch one position back as well
}


Variant::Variant(std::vector<Genotype> const & gts, std::vector<uint32_t> const & hap_calls)
{
  assert(gts.size() > 0);
  for (std::size_t i = 0; i < hap_calls.size(); ++i)
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
Variant::generate_infos()
{
  assert(seqs.size() >= 2);
  infos["RefLen"] = std::to_string(seqs[0].size());

  // Calculate AC, AN, SeqDepth, MaxVS, ABHom, ABHet, ABHomMulti, ABHetMulti
  std::vector<uint32_t> ac(seqs.size() - 1, 0u);
  assert(ac.size() > 0);
  uint32_t an = 0;
  uint64_t seqdepth = 0;
  uint64_t seqdepth_nonref = 0;
  std::pair<uint32_t, uint32_t> het_allele_depth = {0ul, 0ul};   // First is the first call, second is the second call
  std::pair<uint32_t, uint32_t> hom_allele_depth = {0ul, 0ul};   // First is the called allele, second is not the called one

  std::vector<std::pair<uint32_t, uint32_t> > het_multi_allele_depth(seqs.size() - 1, {0u, 0u});
  std::vector<std::pair<uint32_t, uint32_t> > hom_multi_allele_depth(seqs.size(), {0u, 0u});

  // MaxVS
  std::vector<uint16_t> maximum_variant_support(seqs.size() - 1, 0);
  std::vector<double> maximum_alternative_support_ratio(seqs.size() - 1, 0.0);

  for (auto const & sample_call : calls)
  {
    uint32_t total_depth = std::accumulate(sample_call.coverage.cbegin(),
                                           sample_call.coverage.cend(),
                                           static_cast<uint32_t>(sample_call.ambiguous_depth)
      );

    // Skip zero (the reference)
    for (uint16_t c = 1; c < sample_call.coverage.size(); ++c)
    {
      maximum_variant_support[c - 1] = std::max(maximum_variant_support[c - 1], sample_call.coverage[c]);

      if (sample_call.coverage[c] > 0)
      {
        assert (total_depth > 0);
        double const current_ratio = static_cast<double>(sample_call.coverage[c]) / static_cast<double>(total_depth);
        maximum_alternative_support_ratio[c - 1] = std::max(maximum_alternative_support_ratio[c - 1], current_ratio);
      }
    }

    // Check if the there is any call (i.e. there is a non-zero phred score) and ignore the sample call if that is the case
    if (std::find_if(sample_call.phred.cbegin(), sample_call.phred.cend(), [](uint8_t const val){
        return val != 0;
      }) == sample_call.phred.end())
    {
      continue;
    }

    // Get the GT
    std::pair<uint16_t, uint16_t> call = sample_call.get_gt_call();
    assert(seqs.size() == sample_call.coverage.size());
    assert(call.first < seqs.size());
    assert(call.second < seqs.size());

    // Restrict to unique calls, that is calls with only one PHRED=0
    if (std::count(sample_call.phred.begin(), sample_call.phred.end(), 0) == 1)
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
          hom_allele_depth.second += std::accumulate(sample_call.coverage.cbegin(), sample_call.coverage.cend(), 0ull) - sample_call.coverage[call.first];
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
          hom_multi_allele_depth[c].second += std::accumulate(sample_call.coverage.cbegin(), sample_call.coverage.cend(), 0u) - sample_call.coverage[c];
        }
      }
    }

    // SeqDepth and SeqDepth excluding reference allele
    if (sample_call.coverage.size() > 0)
    {
      //assert (sample_call.coverage.size() > 0);
      const uint64_t non_ref_seqdepth = std::accumulate(sample_call.coverage.cbegin() + 1, sample_call.coverage.cend(), 0ull);
      seqdepth_nonref += std::min(static_cast<uint64_t>(10), non_ref_seqdepth); // Dont add too much because phred cannot go very high
      seqdepth += non_ref_seqdepth + sample_call.coverage[0] + sample_call.ambiguous_depth;
    }

    // AN
    {
      an += 2;
    }

    // AC (AF can be derived from AC and AN)
    {
      assert(call.first - 1 < static_cast<long>(ac.size()));
      assert(call.second - 1 < static_cast<long>(ac.size()));

      if (call.first != 0u)
        ++ac[call.first - 1];

      if (call.second != 0u)
        ++ac[call.second - 1];
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


  // Write AC
  {
    assert(ac.size() > 0);
    assert(ac.size() == seqs.size() - 1);
    std::stringstream ac_ss;
    ac_ss << ac[0];

    for (uint32_t e = 1; e < ac.size(); ++e)
      ac_ss << ',' << ac[e];

    infos["AC"] = ac_ss.str();
  }

  // Write AN
  {
    infos["AN"] = std::to_string(an);
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
      ss_abhet << (static_cast<double>(het_allele_depth.second) / static_cast<double>(total_het_depth));
    else
      ss_abhet << "-1";

    infos["ABHet"] = ss_abhet.str();
  }

  // Write ABHom
  {
    std::stringstream ss_abhom;
    ss_abhom.precision(4);
    uint32_t const total_hom_depth = hom_allele_depth.first + hom_allele_depth.second;

    if (total_hom_depth > 0)
      ss_abhom << (static_cast<double>(hom_allele_depth.first) / static_cast<double>(total_hom_depth));
    else
      ss_abhom << "-1";

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

  // Write RAMeanDist
  if (infos.find("RACount") != infos.end() && infos.find("RADist") != infos.end())
  {
    std::vector<uint32_t> count = split_bias_to_numbers(infos["RACount"]);
    std::vector<uint32_t> dist = split_bias_to_numbers(infos["RADist"]);

    assert (count.size() == dist.size());
    assert (count.size() > 0);
    assert (count.size() == seqs.size());
    std::stringstream ss;

    for (unsigned i = 0; i < count.size(); ++i)
    {
      if (i > 0)
        ss << ",";

      if (count[i] > 0)
        ss << static_cast<int>(0.5 + static_cast<double>(dist[i])/static_cast<double>(count[i]));
      else
        ss << "-1";
    }

    infos["RAMeanDist"] = ss.str();
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
          ss_abhet << (static_cast<double>(het_depth_pair.second) / static_cast<double>(total_het_depth));
        else
          ss_abhet << "-1";
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
          ss_abhom << (static_cast<double>(hom_depth_pair.first) / static_cast<double>(total_hom_depth));
        else
          ss_abhom << "-1";
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
    const uint64_t variant_qual = get_qual();

    if (seqdepth_nonref > 0)
      infos["QD"] = std::to_string(static_cast<double>(variant_qual) / static_cast<double>(seqdepth_nonref));
    else
      infos["QD"] = std::string("0.0");
  }
}


void
Variant::trim_sequences(bool const keep_one_match)
{
  add_base_in_front();
  remove_common_suffix(seqs);
  remove_common_prefix(abs_pos, seqs, keep_one_match);
}


void
Variant::remove_uncalled_alleles()
{
  // Check if any AC found in the INFO field
  if (infos.count("AC") == 0)
    generate_infos();

  std::vector<uint16_t> uncalled_alleles = get_list_of_uncalled_alleles(infos["AC"]);

  if (uncalled_alleles.size() == seqs.size() - 1)
    seqs.clear(); // If all alleles are uncalled, just remove the variant
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
    seq.insert(seq.begin(), first_base[0]);

  // Update the abs pos
  --abs_pos;

  return true;
}


bool
Variant::add_base_in_back(bool const add_N)
{
  assert(seqs.size() >= 1);
  uint32_t abs_pos_copy = abs_pos + seqs[0].size();
  uint32_t abs_pos_end = abs_pos_copy + 1;
  std::vector<char> last_base = graph.get_generated_reference_genome(abs_pos_copy, abs_pos_end);

  if (last_base.size() != 1 || abs_pos_copy != (abs_pos + seqs[0].size()) || abs_pos_end != (abs_pos + seqs[0].size() + 1))
    return false; // The base in back could not be extracted";

  if (!add_N && last_base[0] == 'N')
    return false;

  // Insert the new first base in back of all the sequences
  for (auto & seq : seqs)
    seq.push_back(last_base[0]);

  return true;
}


void
Variant::normalize()
{
  remove_common_suffix(seqs);

  // A lambda function which checks if all last bases matches
  auto all_last_bases_match = [&](){
                                for (std::size_t i = 1; i < seqs.size(); ++i)
                                {
                                  if (seqs[i].back() != seqs[0].back())
                                  {
                                    return false;
                                  }
                                  else
                                  {
                                    assert(seqs[i] != seqs[0]); // Make sure the sequences are not identical
                                  }
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

  remove_common_prefix(abs_pos, seqs, false); // Don't keep one match
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
  assert (seqs.size() >= 2);
  return std::find_if(seqs.begin() + 1, seqs.end(), [&](std::vector<char> const & seq){return seq.size() != seqs[0].size();}) == seqs.end();
}


bool
Variant::is_snp() const
{
  auto is_of_size_1 = [](std::vector<char> const & seq){
                        return seq.size() == 1;
                      };

  return seqs.size() >= 2 && seqs.size() <= 4 && std::all_of(seqs.begin(), seqs.end(), is_of_size_1);
}


bool
Variant::is_transition_snp() const
{
  return is_snp() && seqs.size() == 2 && ((seqs[0][0] == 'A' && seqs[1][0] == 'G') ||
                                          (seqs[0][0] == 'G' && seqs[1][0] == 'A') ||
                                          (seqs[0][0] == 'C' && seqs[1][0] == 'T') ||
                                          (seqs[0][0] == 'T' && seqs[1][0] == 'C')
                                          );
}


bool
Variant::is_transversion_snp() const
{
  return is_snp() && seqs.size() == 2 && !is_transition_snp();
}


bool
Variant::is_indel() const
{
  auto is_of_size_1 = [](std::vector<char> const & seq){
                        return seq.size() == 1;
                      };

  auto is_larger_than_1 = [](std::vector<char> const & seq){
                            return seq.size() > 1;
                          };

  return seqs.size() >= 2 && ((seqs[0].size() == 1 && std::all_of(seqs.begin() + 1, seqs.end(), is_larger_than_1)) || // Insertion
                              (seqs[0].size() > 1 && std::all_of(seqs.begin() + 1, seqs.end(), is_of_size_1)) // Deletion
                              );
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

  for (std::size_t i = 0; i < this->seqs.size(); ++i)
  {
    seq_depths.push_back(this->get_seq_depth_of_allele(i));
  }

  return seq_depths;
}


uint64_t
Variant::get_rooted_mapq() const
{
  uint64_t const seq_depth = this->get_seq_depth();
  auto find_it = infos.find("MQ");

  if (find_it != infos.end())
  {
    uint64_t const mapq = std::strtoull(find_it->second.c_str(), NULL, 10);
    return mapq * mapq * seq_depth;
  }

  BOOST_LOG_TRIVIAL(warning) << "[graphtyper::variant] Could not extract MQ.";
  return 60ul * 60ul * seq_depth; // This is nasty, but I think this is the nicest way to go without making things so much more complicated.
}

std::vector<uint64_t>
Variant::get_rooted_mapq_per_allele() const
{
  auto find_it = infos.find("MQperAllele");

  if (find_it == infos.end())
  {
    uint64_t const seq_depth = this->get_seq_depth();
    return std::vector<uint64_t>(this->seqs.size(), 60ul * 60ul * seq_depth);
  }

  std::vector<uint64_t> rooted_mapq_per_allele;
  std::vector<uint32_t> mapq_per_allele = get_strand_bias(this->infos, "MQperAllele");
  assert(mapq_per_allele.size() == seqs.size());
  rooted_mapq_per_allele.reserve(seqs.size());

  for (std::size_t i = 0; i < seqs.size(); ++i)
  {
    uint64_t const seq_depth = this->get_seq_depth_of_allele(i);
    rooted_mapq_per_allele.push_back(static_cast<uint64_t>(mapq_per_allele[i] * mapq_per_allele[i]) * seq_depth);
  }

  return rooted_mapq_per_allele;
}


std::string
Variant::print() const
{
  std::stringstream os;
  auto contig_pos = absolute_pos.get_contig_position(this->abs_pos);
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
  std::size_t non_ones = 0;

  for (auto it = seqs.begin(); it != seqs.end(); ++it)
  {
    assert (it->size() != 0);

    if (it->size() > 1)
      ++non_ones;
  }

  // Check if the variant is a SNP
  if (non_ones == 0)
    return "SG";
  else if (seqs.size() - non_ones == 1) // Check if the variant is an indel
    return "IG";
  else // Otherwise the variant is complex
    return "XG";
}


uint64_t
Variant::get_qual() const
{
  uint64_t variant_qual = 0;

  for (auto const & sample_call : calls)
  {
    if (std::count(sample_call.phred.begin(), sample_call.phred.end(), 0) == 1)
    {
      std::pair<uint16_t, uint16_t> call = sample_call.get_gt_call();

      // Update QUAL if the call is non-reference
      if (call.first != 0 || call.second != 0)
      {
        assert (sample_call.phred.size() > 0);
        variant_qual += sample_call.phred[0];
      }
    }
  }

  return variant_qual;
}


// Variants functions
std::vector<Variant>
break_down_variant(Variant && var, std::size_t const THRESHOLD)
{
  // We need to make sure there is a matching first base
  if (not var.is_with_matching_first_bases())
  {
    if (!var.add_base_in_front())
    {
      // Could not add a first base. Add Ns
      for (auto & seq : var.seqs)
        seq.insert(seq.begin(), 'N');
    }
  }

  std::vector<Variant> broken_down_vars;
  bool const all_same_size = std::find_if(var.seqs.begin() + 1, var.seqs.end(), [&](std::vector<char> const & seq){
      return var.seqs[0].size() != seq.size();
  }) == var.seqs.end();

  //BOOST_LOG_TRIVIAL(debug) << "[graphtyper::variant] all_same_size = " << all_same_size;
  //BOOST_LOG_TRIVIAL(debug) << "[graphtyper::variant] var = " << var.print();

  if (all_same_size)
  {
    std::vector<Variant> new_broken_down_snps = break_multi_snps(std::move(var));
    std::move(new_broken_down_snps.begin(), new_broken_down_snps.end(), std::back_inserter(broken_down_vars));
  }
  else if (var.seqs.size() > 50)
  {
    broken_down_vars.push_back(std::move(var));
  }
  else
  {
    bool const noone_very_large = std::find_if(var.seqs.begin(),
                                               var.seqs.end(),
                                               [](std::vector<char> const & seq){
        return seq.size() > 152;
      }) == var.seqs.end();

    if (noone_very_large)
    {
      assert(var.seqs.size() >= 2);
      assert(var.seqs[0].size() >= 2);
      std::vector<Variant> simpler_variants = simplify_complex_haplotype(std::move(var), THRESHOLD);
      std::move(simpler_variants.begin(), simpler_variants.end(), std::back_inserter(broken_down_vars));
    }
    else
    {
      broken_down_vars.push_back(std::move(var));
    }
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
  assert (var.seqs[0].size() == var.seqs[1].size());
  char first_base = var.seqs[0][0];
  uint32_t pos = var.abs_pos;

  // match_length == -1 means that we have not found a mismatch yet
  long match_length = -1;
  unsigned ref_gaps = 0;

  Variant new_var;
  new_var.abs_pos = pos;
  new_var.seqs = std::vector<std::vector<char> >(var.seqs.size(), std::vector<char>(1, first_base));
  new_var.infos = var.infos;
  new_var.phase = var.phase;

  // Lambda function which handles any matches found
  auto matches_handle =
    [&]()
    {
      // Remove sequences with Ns
      assert(new_var.seqs.size() > 1);
      new_var.seqs.erase(std::remove_if(new_var.seqs.begin() + 1,
                                        new_var.seqs.end(),
                                        [](std::vector<char> const & seq){
        return std::find(seq.begin(), seq.end(), 'N') != seq.end();
      }), new_var.seqs.end());

      if (new_var.seqs.size() <= 1 ||
         std::find(new_var.seqs[0].begin(), new_var.seqs[0].end(), 'N') != new_var.seqs[0].end() /*also check if reference has Ns*/
         )
      {
        return;
      }

      new_var.trim_sequences(true);   // Keep one match

      // Only break down further if THRESHOLD is 1
      if (THRESHOLD == 1)
      {
        bool const all_same_size = std::find_if(new_var.seqs.begin() + 1, new_var.seqs.end(), [&](std::vector<char> const & seq){
          return new_var.seqs[0].size() != seq.size();
        }) == new_var.seqs.end();

        if (all_same_size)
        {
          std::vector<Variant> new_broken_down_snps = break_multi_snps(std::move(new_var));
          std::move(new_broken_down_snps.begin(), new_broken_down_snps.end(), std::back_inserter(new_vars));
        }
        else if (new_var.seqs.size() == 2)
        {
          // This should check for a simple case of indel + snps, e.g. ref=AA, var=ACG. Then I want to output ref1=A, var1=AC, ref2=A, var2=G
          Variant snp_var;
          snp_var.calls = new_var.calls;
          snp_var.infos = new_var.infos;
          snp_var.phase = new_var.phase;

          if (new_var.seqs[0].size() > new_var.seqs[1].size())
          {
            // deletion
            std::size_t const deletion_size = new_var.seqs[0].size() - new_var.seqs[1].size() + 1;
            snp_var.abs_pos = new_var.abs_pos + deletion_size;
            snp_var.seqs.push_back(std::vector<char>(new_var.seqs[0].begin() + deletion_size, new_var.seqs[0].end()));
            snp_var.seqs.push_back(std::vector<char>(new_var.seqs[1].begin() + 1, new_var.seqs[1].end()));
            assert(snp_var.seqs[0].size() == snp_var.seqs[1].size());
            new_var.seqs[0].resize(deletion_size);
            new_var.seqs[1].resize(1);
          }
          else
          {
            // insertion
            assert(new_var.seqs[1].size() > new_var.seqs[0].size());
            std::size_t const insertion_size = new_var.seqs[1].size() - new_var.seqs[0].size() + 1;
            snp_var.abs_pos = new_var.abs_pos + 1;
            snp_var.seqs.push_back(std::vector<char>(new_var.seqs[0].begin() + 1, new_var.seqs[0].end()));
            snp_var.seqs.push_back(std::vector<char>(new_var.seqs[1].begin() + insertion_size, new_var.seqs[1].end()));
            assert(snp_var.seqs[0].size() == snp_var.seqs[1].size());
            new_var.seqs[0].resize(1);
            new_var.seqs[1].resize(insertion_size);
          }

          // Add the indel
          new_vars.push_back(std::move(new_var));

          // Add SNPs, if any
          if (snp_var.seqs[0].size() > 0)
          {
            std::vector<Variant> new_broken_down_snps = break_multi_snps(std::move(snp_var));
            std::move(new_broken_down_snps.begin(), new_broken_down_snps.end(), std::back_inserter(new_vars));
          }
        }
        else if (new_var.seqs.size() > 2)
        {
          // Only apply this if there is exactly one sequence which is of the same size as the reference
          std::size_t const num_same_size_sequences = std::count_if(
            new_var.seqs.begin() + 1,
            new_var.seqs.end(),
            [&](std::vector<char> const & seq){
            return seq.size() == new_var.seqs[0].size();
          }
            );

          if (num_same_size_sequences != 1)
          {
            new_vars.push_back(std::move(new_var));
            return;
          }

          // Find all alleles which are SNPs and call them seperately
          std::vector<Variant> snp_vars;

          for (std::size_t s = 1; s < new_var.seqs.size(); ++s)
          {
            // Only consider variants of the same size
            if (new_var.seqs[0].size() != new_var.seqs[s].size())
              continue;

            uint16_t ref_alt_mismatches = 0;

            for (std::size_t i = 0; i < new_var.seqs[0].size(); ++i)
            {
              assert(i < new_var.seqs[s].size());
              ref_alt_mismatches += new_var.seqs[0][i] != new_var.seqs[s][i];
            }

            assert(ref_alt_mismatches > 0);

            if (ref_alt_mismatches == 1)
            {
              Variant snp_var;
              snp_var.abs_pos = new_var.abs_pos;
              snp_var.infos = new_var.infos;
              snp_var.phase = new_var.phase;

              // Copy reference everywhere except this SNP
              snp_var.seqs = std::vector<std::vector<char> >(new_var.seqs.size(), new_var.seqs[0]);
              snp_var.seqs[s] = new_var.seqs[s];
              find_variant_sequences(snp_var, new_var);
              assert(snp_var.calls.size() == new_var.calls.size());

              // Replace this SNP variants with a reference in the old variant
              Variant old_variant(new_var);
              new_var.seqs[s] = old_variant.seqs[0];
              new_var.calls.clear();   // find_variant_sequences assumes there are no calls in the variant yet
              find_variant_sequences(new_var, old_variant);
              assert(new_var.calls.size() == old_variant.calls.size());
              assert(new_var.seqs.size() == old_variant.seqs.size() - 1);
              --s;   // We just removed a sequence, we need to lower s to keep it correct
              snp_vars.push_back(std::move(snp_var));
            }
          }

          // Trim SNP variants and sort them
          for (auto & a_snp_var : snp_vars)
          {
            a_snp_var.trim_sequences(false);   // don't keep a base in front
            assert(a_snp_var.seqs.size() == 2);
            assert(a_snp_var.seqs[0].size() == 1);
            assert(a_snp_var.seqs[1].size() == 1);
          }

          std::sort(snp_vars.begin(), snp_vars.end(), [](Variant const & v1, Variant const & v2){
            return v1.abs_pos < v2.abs_pos;
          });

          // Add the previous variant first
          new_vars.push_back(std::move(new_var));

          // Add all variants
          if (snp_vars.size() > 0)
            std::move(snp_vars.begin(), snp_vars.end(), std::back_inserter(new_vars));
        }
        else
        {
          new_vars.push_back(std::move(new_var));
        }
      }
      else
      {
        new_vars.push_back(std::move(new_var));
      }
    };

  std::vector<char> const & reference = var.seqs[0];

  for (unsigned i = 1; i < reference.size(); ++i)
  {
    assert(new_var.seqs.size() >= 2);

    if (reference[i] == '-')
      ++ref_gaps;
    else
      new_var.seqs[0].push_back(reference[i]);

    bool all_match = true;

    for (unsigned a = 1; a < new_var.seqs.size(); ++a)
    {
      assert(a < var.seqs.size());

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
      first_base = var.seqs[0][i];
      find_variant_sequences(new_var, var);
      matches_handle();   // Uses new_var, so don't move this anywhere else
      match_length = -1;
      new_var = Variant();   // Create a new one
      new_var.abs_pos = original_pos + i - ref_gaps;
      new_var.seqs = std::vector<std::vector<char> >(var.seqs.size(), std::vector<char>(1, first_base));
      new_var.infos = var.infos;
      new_var.phase = var.phase;
    }
  }

  assert(new_var.seqs.size() >= 2);

  // Process any leftovers which might be remaining
  if (new_var.seqs[0].size() > 0)
  {
    find_variant_sequences(new_var, var);

    if (new_var.seqs.size() >= 2)
      matches_handle();
  }

  return new_vars;
}


std::vector<Variant>
simplify_complex_haplotype(Variant && var, std::size_t const THRESHOLD)
{
  // {
  //   bool const all_same_size = std::find_if(var.seqs.begin() + 1, var.seqs.end(), [&](std::vector<char> const & seq){
  //       return var.seqs[0].size() != seq.size();
  //     }) == var.seqs.end();
  //
  //   // Assume if they are all of the same size they are already aligned
  //   if (all_same_size)
  //     return extract_sequences_from_aligned_variant(std::move(var), THRESHOLD);
  // }

  using namespace seqan;

  auto const & old_seqs = var.seqs;
  assert(old_seqs.size() > 1);
  assert(old_seqs[0].size() > 1);
  char first_base = old_seqs[0][0];
  Align<Dna5String> align;
  resize(rows(align), old_seqs.size());

  for (std::size_t i = 0; i < old_seqs.size(); ++i)
  {
    // Skip first base
    seqan::Dna5String dna_string = std::string(old_seqs[i].begin() + 1, old_seqs[i].end());

    if (seqan::length(dna_string) == 0)
      return {std::move(var)}; // If some string is empty we cannot simplify

    assignSource(row(align, i), dna_string);
  }

  Score<int> score_scheme(2 /*match*/, -2 /*mismatch*/, -1 /*gap extend*/, -7 /*gap open*/);

  globalMsaAlignment(align, score_scheme);
  auto aligned_rows = seqan::rows(align);

  Variant aligned_variant;
  aligned_variant.abs_pos = var.abs_pos;
  aligned_variant.seqs = std::vector<std::vector<char> >(old_seqs.size(), std::vector<char>(1, first_base));
  aligned_variant.infos = var.infos;
  aligned_variant.calls = std::move(var.calls);
  aligned_variant.phase = std::move(var.phase);

  for (unsigned i = 0; i < seqan::length(aligned_rows); ++i)
  {
    std::stringstream ss;
    ss << aligned_rows[i];
    std::string ss_str(ss.str());
    std::move(ss_str.begin(), ss_str.end(), std::back_inserter(aligned_variant.seqs[i]));
  }

  return extract_sequences_from_aligned_variant(std::move(aligned_variant), THRESHOLD);
}


std::vector<Variant>
break_multi_snps(Variant const && var)
{
  uint32_t const pos = var.abs_pos;
  std::vector<std::vector<char> > const & seqs = var.seqs;
  std::vector<Variant> new_vars;

  for (unsigned j = 0; j < seqs[0].size(); ++j)
  {
    std::vector<char> new_seqs(1, seqs[0][j]);
    std::vector<uint16_t> old_phred_to_new_phred(1, 0);

    for (unsigned k = 1; k < seqs.size(); ++k)
    {
      auto find_it = std::find(new_seqs.begin(), new_seqs.end(), seqs[k][j]);

      if (find_it == new_seqs.end())
      {
        old_phred_to_new_phred.push_back(new_seqs.size());
        new_seqs.push_back(seqs[k][j]);
      }
      else
      {
        old_phred_to_new_phred.push_back(std::distance(new_seqs.begin(), find_it));
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
    new_var.phase = var.phase; // Copy the old phase

    for (std::size_t i = 0; i < var.calls.size(); ++i)
    {
      auto const & call = var.calls[i];
      SampleCall new_sample_call;
      new_sample_call.phred = std::vector<uint8_t>(new_var.seqs.size() * (new_var.seqs.size() + 1) / 2, 255u);
      new_sample_call.coverage = std::vector<uint16_t>(new_var.seqs.size(), 0u);
      new_sample_call.ambiguous_depth = call.ambiguous_depth;

      for (uint32_t y = 0; y < seqs.size(); ++y)
      {
        for (uint32_t x = 0; x <= y; ++x)
        {
          uint32_t const index = to_index(x, y);
          uint32_t new_y = old_phred_to_new_phred[y];
          uint32_t new_x = old_phred_to_new_phred[x];

          if (var.phase.size() > 0 && call.phred[index] == 0)
          {
            assert (i < var.phase.size());

            if (new_x > new_y)
              new_var.phase[i] = !var.phase[i];
          }

          if (new_x > new_y)
            std::swap(new_x, new_y);

          uint32_t const new_index = to_index(new_x, new_y);
          new_sample_call.phred[new_index] = std::min(new_sample_call.phred[new_index], call.phred[index]);
        }

        uint32_t new_y = old_phred_to_new_phred[y];

        // Check overflow of coverage
        if (static_cast<uint32_t>(new_sample_call.coverage[new_y]) + static_cast<uint32_t>(call.coverage[y]) < 0xFFFFul)
          new_sample_call.coverage[new_y] += call.coverage[y];
        else
          new_sample_call.coverage[new_y] = 0xFFFFul;
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
  return abs_pos == v.abs_pos && seqs == v.seqs;
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


Variant &
Variant::operator=(Variant const & var)
{
  abs_pos = var.abs_pos;
  seqs = var.seqs;
  calls = var.calls;
  return *this;
}


std::size_t
VariantHash::operator()(Variant const & v) const
{
  assert(v.seqs.size() == 2);
  std::size_t h1 = std::hash<uint32_t>()(v.abs_pos);
  std::size_t h2 = boost::hash_range(v.seqs[0].begin(), v.seqs[0].end());
  std::size_t h3 = boost::hash_range(v.seqs[1].begin(), v.seqs[1].end());
  return h1 ^ (h2 << 1) ^ (h3 + 0x9e3779b9);
}


} // namespace gyper
