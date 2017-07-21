#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <set>

#include <boost/log/trivial.hpp>

#include <graphtyper/graph/graph.hpp> // gyper::Graph
#include <graphtyper/utilities/graph_help_functions.hpp>
#include <graphtyper/utilities/options.hpp> // *gyper::Options::instance()


namespace gyper
{

Haplotype::Haplotype() noexcept
  : gts(0)
  , hap_samples(0)
#ifndef NDEBUG
  , unique_gts(0)
#endif // NDEBUG
{}


void
Haplotype::add_genotype(Genotype && gt)
{
  gts.push_back(std::move(gt));
  explains.push_back(std::bitset<MAX_NUMBER_OF_HAPLOTYPES>(0));
  coverage.push_back(0xFFFFu);
  var_stats.push_back(VarStats(gts.back().num));
}


#ifndef NDEBUG // Only check for duplicate haplotypes in debug mode
void
Haplotype::check_for_duplicate_haplotypes()
{
  // No need to check if there is only one genotype
  if (gts.size() == 1)
    return;

  std::vector<uint32_t> duplicate_haplotypes;
  uint32_t cnum = this->get_genotype_num();

  // Create a set of all sequences in this haplotype, and find out which are duplicates
  std::set<std::vector<char> > haplotype_sequences;

  // Insert reference
  haplotype_sequences.insert(graph.get_sequence_of_a_haplotype_call(gts, 0));

  // Iterate backwards, since we cant to keep the highest haplotype number. Also do not check the reference.
  for (uint32_t i = cnum - 1; i > 0; --i)
  {
    std::vector<char> hap_sequence = graph.get_sequence_of_a_haplotype_call(gts, i);
    auto find_it = haplotype_sequences.find(hap_sequence);

    if (find_it == haplotype_sequences.end())
      haplotype_sequences.insert(std::move(hap_sequence));
    else
      duplicate_haplotypes.push_back(i);
  }

  // Get unique genotypes
  if (duplicate_haplotypes.size() > 0)
  {
    std::cerr << "[graphtyper::haplotype] WARNING: Found a duplicated haplotype, i.e. two haplotypes with the same sequence." << std::endl;
    unique_gts.resize(gts.size());

    for (auto const & c : duplicate_haplotypes)
    {
      uint32_t rem = c;
      uint32_t q = cnum;

      for (uint32_t i = 0; i < gts.size(); ++i)
      {
        assert(gts[i].num > 0);
        q /= gts[i].num;

        assert(q > 0);
        std::bitset<MAX_NUMBER_OF_HAPLOTYPES> new_bits;
        new_bits.set(rem / q); // Set the called allele
        new_bits.flip(); // Flip all the bits, so all are set except the called allele
        unique_gts[i] |= new_bits;
        rem %= q;
      }
    }

    // Mom told me I would just be flippin' burgers all my life. MAMA LOOK AT ME NOW FLIPPIN' 'EM BITSETS
    for (auto uniq_gt : unique_gts)
      uniq_gt.flip();
  }
}
#endif // NDEBUG


void
Haplotype::clear_and_resize_samples(std::size_t const new_size)
{
  hap_samples.clear();
  calls.clear();

  for (std::size_t i = 0; i < new_size; ++i)
  {
    HapSample new_sample;
    std::size_t cnum = get_genotype_num();
    new_sample.log_score.resize(cnum * (cnum + 1) / 2, 0u);

    for (auto const & gt : gts)
      new_sample.gt_coverage.push_back(std::vector<uint16_t>(gt.num, 0u));

    if (Options::instance()->stats.size() > 0)
    {
      new_sample.stats = std::unique_ptr<HapStats>(new HapStats);
      new_sample.stats->hap_coverage.resize(cnum, 0u);
      new_sample.stats->hap_unique_coverage.resize(cnum, 0u);
      new_sample.stats->predecessor.resize(cnum);
      new_sample.stats->successor.resize(cnum);
      new_sample.stats->pair_info.resize(cnum);
    }

    hap_samples.push_back(std::move(new_sample));
  }
}


void
Haplotype::clear()
{
  gts.clear();
  explains.clear();
  coverage.clear();
  hap_samples.clear();
  hap_samples.shrink_to_fit();
  std::vector<HapSample>().swap(hap_samples); // This will always deallocate any memory used by the hap_samples.
  var_stats.clear();
}


uint32_t
Haplotype::get_genotype_num() const
{
  uint32_t cnum = 1;

  for (auto gt_it = gts.cbegin(); gt_it != gts.cend(); ++gt_it)
    cnum *= gt_it->num;

  return cnum;
}


bool
Haplotype::has_too_many_genotypes() const
{
  uint32_t cnum = 1;

  for (uint32_t i = 0; i < gts.size(); ++i)
  {
    cnum *= gts[i].num;

    if (Options::instance()->phased_output)
    {
      if (cnum >= (MAX_NUMBER_OF_HAPLOTYPES - 1))
        return true;
    }
    else if (Options::instance()->is_segment_calling)
    {
      if (cnum >= 337)
        return true;
    }
    else if (cnum > 97u) // 97 is a prime
    {
      return true;
    }
  }

  return false;
}


void
Haplotype::add_coverage(uint32_t const i, uint16_t const c)
{
  // i is local genotype id
  // c is coverage for this local genotype id
  assert(i < gts.size());
  assert(i < coverage.size());

  if (c == 0xFFFEu)
  {
    coverage[i] = c;
  }
  else
  {
    assert(c < gts[i].num);
    assert(c < 0xFFFEu);

    if (coverage[i] == 0xFFFFu)
      coverage[i] = c;   // Coverage has not been set, so we set it to c
    else if (coverage[i] != c)
      coverage[i] = 0xFFFEu;   // Multiple matches
  }
}


void
Haplotype::clipped_reads_to_stats(bool const fully_aligned)
{
  if (!fully_aligned)
  {
    assert (var_stats.size() == gts.size());
    assert (var_stats.size() == coverage.size());

    for (std::size_t i = 0; i < var_stats.size(); ++i)
    {
      // coverage[i] == 0xFFFFu means there is no coverage
      if (coverage[i] != 0xFFFFu) // && coverage[i] != 0xFFFEu)
        ++var_stats[i].clipped_reads;
    }
  }
}


void
Haplotype::graph_complexity_to_stats()
{
  assert(var_stats.size() == gts.size());

  for (std::size_t i = 0; i < gts.size(); ++i)
    var_stats[i].graph_complexity = graph.get_10log10_num_paths(gts[i].first_variant_node);
}


void
Haplotype::mapq_to_stats(uint8_t const new_mapq)
{
  assert(var_stats.size() == gts.size());
  assert(var_stats.size() == coverage.size());

  for (std::size_t i = 0; i < var_stats.size(); ++i)
  {
    // coverage[i] == 0xFFFFu means there is no coverage
    if (coverage[i] != 0xFFFFu && coverage[i] != 0xFFFEu)
      var_stats[i].add_mapq(coverage[i], new_mapq);
  }
}


void
Haplotype::realignment_to_stats(bool const is_unaligned_read,
                                bool const is_originally_clipped,
                                uint32_t const original_pos,
                                uint32_t const new_pos
  )
{
  assert(var_stats.size() == gts.size());
  assert(var_stats.size() == coverage.size());

  if (is_unaligned_read)
  {
    for (std::size_t i = 0; i < var_stats.size(); ++i)
    {
      ++var_stats[i].unaligned_reads;
    }
  }
  else
  {
    for (std::size_t i = 0; i < var_stats.size(); ++i)
    {
      if (coverage[i] != 0xFFFFu && coverage[i] != 0xFFFEu)
      {
        assert(coverage[i] < var_stats[i].originally_clipped.size());

        if (is_originally_clipped)
          ++var_stats[i].originally_clipped[coverage[i]];

        var_stats[i].add_realignment_distance(coverage[i] /*allele_id*/, original_pos, new_pos);
      }
    }
  }
}


void
Haplotype::strand_to_stats(bool const forward_strand, bool const is_first_in_pair)
{
  assert(var_stats.size() == gts.size());
  assert(var_stats.size() == coverage.size());

  for (std::size_t i = 0; i < var_stats.size(); ++i)
  {
    // coverage[i] == 0xFFFFu means there is no coverage
    if (coverage[i] != 0xFFFFu && coverage[i] != 0xFFFEu)
    {
      assert(coverage[i] < var_stats[i].r1_strand_forward.size());
      assert(coverage[i] < var_stats[i].r2_strand_forward.size());

      assert(var_stats[i].r1_strand_forward.size() == var_stats[i].r1_strand_reverse.size());
      assert(var_stats[i].r2_strand_forward.size() == var_stats[i].r2_strand_reverse.size());

      if (forward_strand)
      {
        if (is_first_in_pair)
          ++var_stats[i].r1_strand_forward[coverage[i]];
        else
          ++var_stats[i].r2_strand_forward[coverage[i]];
      }
      else
      {
        if (is_first_in_pair)
          ++var_stats[i].r1_strand_reverse[coverage[i]];
        else
          ++var_stats[i].r2_strand_reverse[coverage[i]];
      }
    }
  }
}


void
Haplotype::coverage_to_gts(std::size_t const pn_index)
{
  assert(coverage.size() == gts.size());
  assert(pn_index < hap_samples.size());

  // Loop over the coverage values and add the coverage to its respected genotype
  for (unsigned i = 0; i < coverage.size(); ++i)
  {
    // coverage[i] == 0xFFFFu means there is no coverage
    // coverage[i] == 0xFFFEu means there are more than one coverage match
    if (coverage[i] < 0xFFFEu) // same as coverage[i] != 0xFFFFu && coverage[i] != 0xFFFEu
    {
      assert(i < hap_samples[pn_index].gt_coverage.size());
      assert(coverage[i] < gts[i].num);
      assert(coverage[i] < hap_samples[pn_index].gt_coverage[i].size());

      // Check for overflow
      if (hap_samples[pn_index].gt_coverage[i][coverage[i]] < 0xFFFFul)
        ++hap_samples[pn_index].gt_coverage[i][coverage[i]];
    }
    else if (coverage[i] == 0xFFFEu)
    {
      // More than one coverage match means we need to update ambiguous depth
      if (hap_samples[pn_index].ambiguous_depth < 255)
        ++hap_samples[pn_index].ambiguous_depth;
    }
  }

  // Reset coverage
  coverage = std::vector<uint16_t>(gts.size(), 0xFFFFu);
}


std::bitset<MAX_NUMBER_OF_HAPLOTYPES>
Haplotype::explain_to_path_explain()
{
  uint32_t const cnum = get_genotype_num();

  // Find gts with no explanation
  for (auto & e : explains)
  {
    if (e.none())
      e.set(); // Flip 'em all, cause they can all explain this read
  }

  std::bitset<MAX_NUMBER_OF_HAPLOTYPES> path_explains = find_which_haplotypes_explain_the_read(cnum);

  // Clear all bitsets
  for (std::size_t i = 0; i < explains.size(); ++i)
    explains[i] = std::bitset<MAX_NUMBER_OF_HAPLOTYPES>(0); // This is better than using reset()!

  return path_explains;
}


void
Haplotype::add_explanation(uint32_t const i, std::bitset<MAX_NUMBER_OF_HAPLOTYPES> const & e)
{
  // i is local genotype id
  // e is explain bitset for this local genotype id
  assert(i < explains.size());
  assert(gts.size() == explains.size());
  explains[i] |= e;
}


std::bitset<MAX_NUMBER_OF_HAPLOTYPES>
Haplotype::find_which_haplotypes_explain_the_read(uint32_t const num) const
{
  std::bitset<MAX_NUMBER_OF_HAPLOTYPES> haplotype_explains(0);

  for (uint32_t c = 0; c < num; ++c)
  {
    uint32_t rem = c;
    uint32_t q = num;
    bool all_gts_explain = true;

    for (uint32_t i = 0; i < explains.size(); ++i)
    {
      assert (gts[i].num != 0);
      q /= gts[i].num;
      assert(q > 0);

      if (not explains[i].test(rem / q))
      {
        all_gts_explain = false;
        break;
      }

      rem %= q;
    }

    if (all_gts_explain)
      haplotype_explains.set(c);
  }

  return haplotype_explains;
}


std::vector<uint16_t>
Haplotype::find_with_how_many_errors_haplotypes_explain_the_read(uint32_t const cnum) const
{
  std::vector<uint16_t> haplotype_errors;
  haplotype_errors.reserve(cnum);
  assert(gts.size() == explains.size());

  for (uint32_t c = 0; c < cnum; ++c)
  {
    uint32_t rem = c;
    uint32_t q = cnum;
    uint16_t errors = 0;

    for (uint32_t i = 0; i < explains.size(); ++i)
    {
      assert(gts[i].num > 0);
      assert(explains[i].any());
      q /= gts[i].num;

      if (!explains[i].test(rem / q))
        ++errors;

      rem %= q;
    }

    haplotype_errors.push_back(errors);
  }

  assert(haplotype_errors.size() == cnum);
  return haplotype_errors;
}


void
Haplotype::explain_to_score(std::size_t const pn_index,
                            bool const has_low_quality_snp,
                            bool const non_unique_paths,
                            uint8_t const mapq,
                            bool const fully_aligned,
                            std::size_t const mismatches)
{
  // Update log_score
  uint16_t epsilon_exponent;

  if (Options::instance()->is_segment_calling)
  {
    uint16_t constexpr MISMATCH_PENALTY = 2;
    uint16_t constexpr MAX_MISMATCHES = 2;
    uint16_t constexpr NOT_FULLY_ALIGNED_READ_PENALTY = 3;
    uint16_t constexpr BAD_MAPQ_PENALTY = 1;
    uint16_t constexpr MAPQ_ZERO_PENALTY = 1;
    uint16_t constexpr MAX_PENALTY = MISMATCH_PENALTY * MAX_MISMATCHES + NOT_FULLY_ALIGNED_READ_PENALTY + BAD_MAPQ_PENALTY + MAPQ_ZERO_PENALTY;

    epsilon_exponent = std::max(MAX_PENALTY, static_cast<uint16_t>(12));

    // Penalize for mismatches with regards to the graph
    epsilon_exponent -= MISMATCH_PENALTY * std::min(mismatches, static_cast<std::size_t>(MAX_MISMATCHES));

    // Penalize reads which do not fully align to the graph
    epsilon_exponent -= NOT_FULLY_ALIGNED_READ_PENALTY;

    if (non_unique_paths /*local graph mapping qual*/ || mapq <= 25 /*global mapping qual*/)
      epsilon_exponent -= BAD_MAPQ_PENALTY;

    // Further increase epsilon if the mapping quality is zero
    if (mapq == 0)
      epsilon_exponent -= MAPQ_ZERO_PENALTY;
  }
  else
  {
    uint16_t constexpr LOW_QUALITY_SNP_PENALTY = 1;
    uint16_t constexpr BAD_MAPQ_PENALTY = 2;
    uint16_t constexpr MAPQ_ZERO_PENALTY = 1;
    uint16_t constexpr NOT_FULLY_ALIGNED_READ_PENALTY = 3;
    uint16_t constexpr MAX_PENALTY = LOW_QUALITY_SNP_PENALTY + BAD_MAPQ_PENALTY + MAPQ_ZERO_PENALTY + NOT_FULLY_ALIGNED_READ_PENALTY;

    // Make sure we do not underflow
    epsilon_exponent = std::max(MAX_PENALTY, Options::instance()->epsilon_0_exponent); // epsilon_0 = (1/2)^{epsilon_0_exponent} is the lowest epsilon

    if (has_low_quality_snp)
      epsilon_exponent -= LOW_QUALITY_SNP_PENALTY;

    if (non_unique_paths /*local graph mapping qual*/ || mapq <= 25 /*global mapping qual*/)
      epsilon_exponent -= BAD_MAPQ_PENALTY;

    // Further increase epsilon if the mapping quality is zero
    if (mapq == 0)
      epsilon_exponent -= MAPQ_ZERO_PENALTY;

    // Decrease epsilon if the path was not fully aligned, i.e. it has been clipped
    if (!fully_aligned)
      epsilon_exponent -= NOT_FULLY_ALIGNED_READ_PENALTY;
  }

  // Make sure epsilon is not too low
  epsilon_exponent = std::max(epsilon_exponent, static_cast<uint16_t>(9));

  // Get how many number of paths are in this haplotype
  uint32_t const cnum = get_genotype_num();

  // Find gts with no explanation
  assert(gts.size() == explains.size());

  for (uint32_t i = 0; i < explains.size(); ++i)
  {
    if (explains[i].none())
      explains[i].set(); // Flip 'em all, cause they can all explain this read
  }

  assert(gts.size() == explains.size());
  std::vector<uint16_t> haplotype_errors = find_with_how_many_errors_haplotypes_explain_the_read(cnum);

  // Update log scores
  assert (pn_index < hap_samples.size());
  auto & hap_sample = hap_samples[pn_index];

  // If we are generating statistics, we store the number of reads that support haplotype
  if (Options::instance()->stats.size() > 0)
  {
    assert (hap_sample.stats);
    // assert (haplotype_errors.size() == hap_sample.stats->hap_anti_coverage.size());

    if (std::count(haplotype_errors.begin(), haplotype_errors.end(), 0) == 1)
    {
      // Update unique coverage
      auto it = std::find(haplotype_errors.begin(), haplotype_errors.end(), 0);
      assert (it != haplotype_errors.end());
      size_t const i = std::distance(haplotype_errors.begin(), it);

      if (hap_sample.stats->hap_unique_coverage[i] < 0xFFul)
        ++hap_sample.stats->hap_unique_coverage[i];

      assert (std::find(it + 1, haplotype_errors.end(), 0) == haplotype_errors.end());

      // Also update total coverage
      if (hap_sample.stats->hap_coverage[i] < 0xFFul)
        ++hap_sample.stats->hap_coverage[i];
    }
    else
    {
      // Update total coverage
      for (size_t i = 0; i < haplotype_errors.size(); ++i)
      {
        if (haplotype_errors[i] == 0 && hap_sample.stats->hap_coverage[i] < 0xFFul)
          ++hap_sample.stats->hap_coverage[i];
      }
    }
  }

  // Stop adding reads if the maximum log score is maxed out (this can happend if read depth is more than about 6000x)
  if (hap_sample.max_log_score < 0xFFFFul - epsilon_exponent)
  {
    hap_sample.max_log_score += epsilon_exponent;
    uint32_t i = 0;

    for (uint32_t y = 0; y < cnum; ++y)
    {
      for (uint32_t x = 0; x <= y; ++x, ++i)
      {
        assert(to_index(x, y) < hap_sample.log_score.size());
        assert(i == to_index(x, y));

        if (haplotype_errors[x] == 0 && haplotype_errors[y] == 0)
          hap_sample.log_score[i] += epsilon_exponent;
        else if (haplotype_errors[x] == 0 || haplotype_errors[y] == 0)
          hap_sample.log_score[i] += epsilon_exponent - 1;
        else if (haplotype_errors[x] == 1 && haplotype_errors[y] == 1)
          hap_sample.log_score[i] += 4;
        else if (haplotype_errors[x] == 1 || haplotype_errors[y] == 1)
          hap_sample.log_score[i] += 3;
        else if (haplotype_errors[x] == 2 && haplotype_errors[y] == 2)
          hap_sample.log_score[i] += 2;
        else if (haplotype_errors[x] == 2 || haplotype_errors[y] == 2)
          hap_sample.log_score[i] += 1;
        // else the log_score does not change
      }
    }
  }

  // Clear all bitsets
  assert(gts.size() == explains.size());
  for (std::size_t i = 0; i < gts.size(); ++i)
    explains[i] = std::bitset<MAX_NUMBER_OF_HAPLOTYPES>(0);
}


void
Haplotype::update_max_log_score()
{
  for (auto & hap_sample : hap_samples)
  {
    assert (hap_sample.log_score.size() > 0);
    auto max_it = std::max_element(hap_sample.log_score.begin(), hap_sample.log_score.end());
    assert (max_it != hap_sample.log_score.end());
    hap_sample.max_log_score = *max_it;
    calls.push_back(to_pair(std::distance(hap_sample.log_score.begin(), max_it)));
  }

  assert (calls.size() == hap_samples.size());
}


std::vector<uint32_t>
Haplotype::get_genotype_ids() const
{
  std::vector<uint32_t> gt_ids;

  for (unsigned i = 0; i < gts.size(); ++i)
    gt_ids.push_back(gts[i].id);

  return gt_ids;
}


std::vector<uint16_t>
Haplotype::get_haplotype_calls() const
{
  std::vector<uint16_t> all_calls = {0, 0};
  std::size_t const cnum = get_genotype_num();
  std::size_t const MINIMUM_VARIANT_SUPPORT = Options::instance()->minimum_variant_support;

  for (auto const & hap_sample : hap_samples)
  {
    // Check if coverage is above the minimum_variant_support.
    // Note: It is possible to get false negatives from this calculations, but no false positives.
    auto coverage_above_cutoff = [&](uint16_t call)
    {
      uint32_t q = cnum;
      assert (gts.size() == hap_sample.gt_coverage.size());

      for (uint32_t i = 0; i < hap_sample.gt_coverage.size(); ++i)
      {
        q /= gts[i].num;
        assert (q != 0);
        assert (call / q < hap_sample.gt_coverage[i].size());

        // Require 'MINIMUM_VARIANT_SUPPORT'-many reads to overlap at least one of the alleles of the genotype call
        if (static_cast<uint32_t>(hap_sample.gt_coverage[i][call / q]) >= MINIMUM_VARIANT_SUPPORT)
          return true;

        call %= q;
      }

      return false;
    };

    if (Options::instance()->output_all_variants)
    {
      for (uint16_t call = 1; call < cnum; ++call)
      {
        if (!hap_sample.stats || hap_sample.stats->hap_unique_coverage[call] >= MINIMUM_VARIANT_SUPPORT)
          all_calls.push_back(call);
      }

      continue;
    }

    // Favor homozygous calls
    uint16_t call1 = 0xFFFFul;
    uint16_t call2 = 0xFFFFul;
    uint16_t local_max_log_score = 0;

    // Homozugous first, because we want to have the least amount af haplotypes extracted
    for (uint16_t y = 0; y < cnum; ++y)
    {
      if (hap_sample.log_score[to_index(y, y)] > local_max_log_score)
      {
        local_max_log_score = hap_sample.log_score[to_index(y, y)];
        call1 = y;
        call2 = y;
      }
    }

    for (uint16_t y = 1; y < cnum; ++y)
    {
      for (uint16_t x = 0; x < y; ++x)
      {
        if (hap_sample.log_score[to_index(x, y)] > local_max_log_score)
        {
          local_max_log_score = hap_sample.log_score[to_index(x, y)];
          call1 = x;
          call2 = y;
        }
      }
    }

    // call1 == 0xFFFFul means no read overlapped this haplotype
    if (call1 != 0xFFFFul)
    {
      if (call1 != 0 && coverage_above_cutoff(call1))
        all_calls.push_back(call1);

      if (call2 != 0 && coverage_above_cutoff(call2))
        all_calls.push_back(call2);
    }
  }

  return all_calls;
}


uint32_t
Haplotype::best_score_of_a_path(std::size_t const pn_index,
                                std::bitset<MAX_NUMBER_OF_HAPLOTYPES> const & e1,
                                std::bitset<MAX_NUMBER_OF_HAPLOTYPES> const & e2
                                /*std::vector<uint16_t> const & max_scores_per_y*/
                                ) const
{
  assert (pn_index < hap_samples.size());
  assert (hap_samples[pn_index].log_score.size() > 0);
  // std::size_t const e1_count = e1.count();
  // std::size_t const e2_count = e2.count();

  // assert (pn_index < calls.size());
  // assert (calls[pn_index].first != 0xFFFFu);
  // assert (calls[pn_index].second != 0xFFFFu);

  if (e1.none() or e2.none())
    return 0;

  auto const & hap_sample = hap_samples[pn_index];
  auto const & call = calls[pn_index];

  if ((e1[call.first] and e2[call.second]) or (e2[call.first] and e1[call.second]))
    return hap_sample.max_log_score;

  // if (e1_count * e2_count <= 4000)
  // {
  uint16_t max_path_score = 0;
  std::vector<uint16_t> e2_bit_on;
  uint32_t const cnum = this->get_genotype_num();

  for (uint16_t idx2 = 0; idx2 < cnum; ++idx2)
  {
    if (e2.test(idx2))
      e2_bit_on.push_back(idx2);
  }

  for (uint16_t idx1 = 0; idx1 < cnum; ++idx1)
  {
    if (!e1.test(idx1))
      continue;

    for (auto const idx2 : e2_bit_on)
    {
      uint16_t score;

      if (idx1 <= idx2)
        score = hap_sample.log_score[to_index(idx1, idx2)];
      else
        score = hap_sample.log_score[to_index(idx2, idx1)];

      if (score > max_path_score)
      {
        if (score == hap_sample.max_log_score)
          return hap_sample.max_log_score;

        max_path_score = score;
      }

    }
  }

  // std::cout << "max_path_score = " << max_path_score << std::endl;
  return max_path_score;

  // if (idx1 <= idx2)
  //   return hap_sample.log_score[to_index(idx1, idx2)];
  // else
  //   return hap_sample.log_score[to_index(idx2, idx1)];
  // }
  // else
  // {
  //   // std::cout << "High amount of possibilities: " << (e1.count() * e2.count()) << std::endl;
  //   // Approximate score for this case, otherwise the process takes too much time
  //   if (e1[call.first] or e2[call.first] or e1[call.second] or e2[call.second])
  //     return std::max(static_cast<int32_t>(0), static_cast<int32_t>(hap_sample.max_log_score) - 2 * 15);
  //   else
  //     return std::max(static_cast<int32_t>(0), static_cast<int32_t>(hap_sample.max_log_score) - 100);
  // }

  // // uint16_t max_path_score = std::max(static_cast<int32_t>(0), static_cast<int32_t>(hap_sample.max_log_score - 100));
  // uint16_t max_path_score = 0;
  // uint32_t const cnum = this->get_genotype_num();
  // // assert(max_scores_per_y.size() == cnum);
  //
  // for (uint32_t y = 0; y < cnum; ++y)
  // {
  //   assert (to_index(0, y) < hap_sample.log_score.size());
  //   assert (to_index(y, y) < hap_sample.log_score.size());
  //   // if (max_path_score > max_scores_per_y[y])
  //   //   continue;
  //
  //   // We split this check into two loops
  //   // if (!e1.test(y) && !e2.test(y))
  //   //   continue;
  //
  //   if (e1.test(y))
  //   {
  //     for (uint32_t x = 0; x <= y; ++x)
  //     {
  //       // Only use solutions which are explained by the bitsets
  //       if (e2.test(x))
  //       {
  //         uint16_t const score = hap_sample.log_score[to_index(x, y)];
  //
  //         if (score > max_path_score)
  //         {
  //           if (score == hap_sample.max_log_score)
  //             return score;
  //
  //           max_path_score = score;
  //         }
  //       }
  //     }
  //   }
  //
  //   if (e2.test(y))
  //   {
  //     for (uint32_t x = 0; x <= y; ++x)
  //     {
  //       // Only use solutions which are explained by the bitsets
  //       if (e1.test(x))
  //       {
  //         uint16_t const score = hap_sample.log_score[to_index(x, y)];
  //
  //         if (score > max_path_score)
  //         {
  //           if (score == hap_sample.max_log_score)
  //             return score;
  //
  //           max_path_score = score;
  //         }
  //       }
  //     }
  //   }
  // }
  //
  // return max_path_score;
}


} // namespace gyper
