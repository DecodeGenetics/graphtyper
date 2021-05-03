#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <set>
#include <vector>

#include <boost/log/trivial.hpp>

#include <graphtyper/graph/graph.hpp> // gyper::Graph
#include <graphtyper/utilities/graph_help_functions.hpp>
#include <graphtyper/utilities/options.hpp> // *gyper::Options::instance()


namespace gyper
{

uint16_t constexpr Haplotype::NO_COVERAGE;
uint16_t constexpr Haplotype::MULTI_ALT_COVERAGE;
uint16_t constexpr Haplotype::MULTI_REF_COVERAGE;


void
HapSample::increment_ambiguous_depth()
{
  // Check for overflow
  if (ambiguous_depth < 0xFFu)
    ++ambiguous_depth;
}


void
HapSample::increment_ambiguous_depth_alt()
{
  // Check for overflow
  if (ambiguous_depth_alt < 0xFFu)
    ++ambiguous_depth_alt;
}


void
HapSample::increment_allele_depth(std::size_t const allele_index)
{
  // Check for overflow
  if (gt_coverage[allele_index] < 0xFFFFu)
    ++gt_coverage[allele_index];
}


void
HapSample::increment_alt_proper_pair_depth()
{
  if (alt_proper_pair_depth < 0xFFu)
    ++alt_proper_pair_depth;
}


Haplotype::Haplotype() noexcept
  : hap_samples(0)
{}


void
Haplotype::add_genotype(Genotype && gt)
{
  var_stats = VarStats(gt.num);
  this->gt = std::move(gt);
  //explains = std::bitset<MAX_NUMBER_OF_HAPLOTYPES>(0);
  //coverage = Haplotype::NO_COVERAGE;
}


/*
void
Haplotype::check_for_duplicate_haplotypes()
{
  std::vector<std::bitset<MAX_NUMBER_OF_HAPLOTYPES> > unique_gts; // One per gt

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
    BOOST_LOG_TRIVIAL(warning) << __HERE__ << " Found a duplicated haplotype, "
                               << "i.e. two haplotypes with the same sequence.";
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
*/


void
Haplotype::clear_and_resize_samples(std::size_t const new_size)
{
  hap_samples.clear();
  //calls.clear();

  for (std::size_t i = 0; i < new_size; ++i)
  {
    HapSample new_sample;
    std::size_t const cnum = get_genotype_num();
    new_sample.log_score.resize(cnum * (cnum + 1) / 2, 0u);

#ifdef GT_DEV
    new_sample.connections.resize(cnum);
#endif // GT_DEV

    new_sample.gt_coverage = std::vector<uint16_t>(gt.num, 0u);

#ifndef NDEBUG
    if (Options::instance()->stats.size() > 0)
    {
      new_sample.stats = std::unique_ptr<HapStats>(new HapStats);
      new_sample.stats->hap_coverage.resize(cnum, 0u);
      new_sample.stats->hap_unique_coverage.resize(cnum, 0u);
      new_sample.stats->pair_info.resize(cnum);
    }
#endif // NDEBUG

    hap_samples.push_back(std::move(new_sample));
  }
}


void
Haplotype::clear()
{
  gt = Genotype();
  explains.clear();
  coverage = NO_COVERAGE;
  hap_samples.clear();
  hap_samples.shrink_to_fit();
  std::vector<HapSample>().swap(hap_samples); // This will always deallocate any memory used by the hap_samples.
  var_stats = VarStats();
}


uint32_t
Haplotype::get_genotype_num() const
{
  return gt.num;
}


bool
Haplotype::has_too_many_genotypes() const
{
  //long cnum = 1;
  //
  //for (unsigned i = 0; i < gts.size(); ++i)
  //{
  //  cnum *= static_cast<long>(gts[i].num);
  //
  //  if (cnum > 337u)
  //    return true;
  //}

  return false;
}


void
Haplotype::add_coverage(uint32_t const /*i*/, uint16_t const c)
{

  switch (coverage)
  {
  case NO_COVERAGE:
  {
    // Coverage has not been set, so we set it (uniquely) to c
    coverage = c;
    break;
  }

  case MULTI_ALT_COVERAGE:
  {
    if (c == 0)
    {
      // Coverage for the reference allele
      coverage = MULTI_REF_COVERAGE;
    }

    break;
  }

  case MULTI_REF_COVERAGE:
  {
    break;
  }

  default:
  {
    if (coverage != c)
    {
      // Multiple alleles are covered
      if (coverage == 0 || c == 0)
      {
        // One of the allele is the reference allele
        coverage = MULTI_REF_COVERAGE;
      }
      else
      {
        // None of the alleles are the reference allele
        coverage = MULTI_ALT_COVERAGE;
      }
    }

    break;
  }
  }
}


void
Haplotype::clipped_reads_to_stats(int const clipped_bp, int const read_length)
{
  if (clipped_bp == 0)
    return;

  long const scaled_clipped_bp = (clipped_bp * 1000l) / read_length;

  if (coverage != NO_COVERAGE)
    ++var_stats.clipped_reads;

  if (coverage < MULTI_REF_COVERAGE)   // true when the coverage is unique to a particular allele
  {
    assert(coverage < var_stats.per_allele.size());
    var_stats.per_allele[coverage].clipped_bp += scaled_clipped_bp;
  }
}


void
Haplotype::mapq_to_stats(uint8_t const new_mapq)
{
  if (new_mapq == 255)
    return; // means that MQ is unavailable

  uint64_t const mapq_squared = static_cast<uint64_t>(new_mapq) * static_cast<uint64_t>(new_mapq);

  if (coverage != NO_COVERAGE)
    var_stats.mapq_squared += mapq_squared;

  if (coverage < MULTI_REF_COVERAGE)   // true when the coverage is unique to a particular allele
  {
    assert(coverage < var_stats.per_allele.size());
    var_stats.per_allele[coverage].mapq_squared += mapq_squared;
  }
}


void
Haplotype::strand_to_stats(uint16_t const flags)
{
  bool const forward_strand = (flags & IS_SEQ_REVERSED) == 0;
  bool const is_first_in_pair = (flags & IS_FIRST_IN_PAIR) != 0;

  if (coverage < MULTI_REF_COVERAGE)   // true when the coverage is unique to a particular allele
  {
    assert(coverage < var_stats.read_strand.size());

    if (forward_strand)
    {
      if (is_first_in_pair)
        ++var_stats.read_strand[coverage].r1_forward;
      else
        ++var_stats.read_strand[coverage].r2_forward;
    }
    else
    {
      if (is_first_in_pair)
        ++var_stats.read_strand[coverage].r1_reverse;
      else
        ++var_stats.read_strand[coverage].r2_reverse;
    }
  }
}


void
Haplotype::mismatches_to_stats(uint8_t const mismatches, int const read_length)
{
  if (mismatches == 0)
    return;

  long const scaled_mismatches = (mismatches * 1000l) / read_length;

  if (coverage < MULTI_REF_COVERAGE)   // true when the coverage is unique to a particular allele
  {
    assert(coverage < var_stats.per_allele.size());
    var_stats.per_allele[coverage].mismatches += scaled_mismatches;
  }
}


void
Haplotype::score_diff_to_stats(uint8_t const score_diff)
{
  if (score_diff == 0)
    return;

  if (coverage < MULTI_REF_COVERAGE)   // true when the coverage is unique to a particular allele
  {
    assert(coverage < var_stats.per_allele.size());
    var_stats.per_allele[coverage].score_diff += score_diff;
  }
}


void
Haplotype::coverage_to_gts(std::size_t const pn_index, bool const is_proper_pair)
{
  //assert(coverage.size() == gts.size());
  assert(pn_index < hap_samples.size());
  auto & hap_sample = hap_samples[pn_index];

  // Loop over the coverage values and add the coverage to its respected genotype
  switch (coverage)
  {
  case NO_COVERAGE:
  {
    // Do nothing
    break;
  }

  case MULTI_REF_COVERAGE:
  {
    // More than one allele has coverage means we need to update ambiguous depth, including the reference allele
    hap_sample.increment_ambiguous_depth();
    break;
  }

  case MULTI_ALT_COVERAGE:
  {
    // The reference alleles does not have coverage but 2 or more alternative alleles have coverage
    hap_sample.increment_ambiguous_depth();
    hap_sample.increment_ambiguous_depth_alt();

    if (is_proper_pair)
      hap_sample.increment_alt_proper_pair_depth();

    break;
  }


  default:
  {
    // when the coverage is unique to a particular allele
    assert(coverage < gt.num);
    hap_sample.increment_allele_depth(coverage /*allele_index*/);

    if (coverage > 0 && is_proper_pair)
      hap_sample.increment_alt_proper_pair_depth();

    break;
  }

  }
}


/*
std::bitset<MAX_NUMBER_OF_HAPLOTYPES>
Haplotype::explain_to_path_explain()
{
  uint32_t const cnum = get_genotype_num();

  // Find gts with no explanation
  if (explains.none())
    explains.set(); // Flip 'em all, cause they can all explain this read

  std::bitset<MAX_NUMBER_OF_HAPLOTYPES> path_explains = find_which_haplotypes_explain_the_read(cnum);

  // Clear all bitsets
  explains = std::bitset<MAX_NUMBER_OF_HAPLOTYPES>(0); // This is better than using reset()!

  return path_explains;
}
*/

/*
void
Haplotype::add_explanation(uint32_t const i, std::bitset<MAX_NUMBER_OF_HAPLOTYPES> const & e)
{
  // i is local genotype id
  // e is explain bitset for this local genotype id
  //assert(i < explains.size());
  //assert(gts.size() == explains.size());
  explains |= e;
}
*/

/*
std::bitset<MAX_NUMBER_OF_HAPLOTYPES>
Haplotype::find_which_haplotypes_explain_the_read(uint32_t const num) const
{
  return explains;
  //std::bitset<MAX_NUMBER_OF_HAPLOTYPES> haplotype_explains(0);
  //
  //for (uint32_t c = 0; c < num; ++c)
  //{
  //  uint32_t rem = c;
  //  uint32_t q = num;
  //  bool all_gts_explain = true;
  //
  //  for (uint32_t i = 0; i < explains.size(); ++i)
  //  {
  //    assert(gts[i].num != 0);
  //    q /= gts[i].num;
  //    assert(q > 0);
  //
  //    if (not explains[i].test(rem / q))
  //    {
  //      all_gts_explain = false;
  //      break;
  //    }
  //
  //    rem %= q;
  //  }
  //
  //  if (all_gts_explain)
  //    haplotype_explains.set(c);
  //}
  //
  //return haplotype_explains;
}


std::vector<uint16_t>
Haplotype::find_with_how_many_errors_haplotypes_explain_the_read(uint32_t const cnum) const
{
  std::vector<uint16_t> haplotype_errors;
  haplotype_errors.reserve(cnum);

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
*/

void
Haplotype::explain_to_score(std::size_t const pn_index,
                            bool const non_unique_paths,
                            uint16_t const flags,
                            bool const fully_aligned,
                            bool const is_read_overlapping,
                            bool const is_low_qual,
                            std::size_t const mismatches)
{
  long constexpr MISMATCH_PENALTY = 1; // penalty for every mismatch
  long constexpr NON_UNIQUE_PATHS_PENALTY = 3; // penalty for read matching to multiple places
  long constexpr BAD_MAPQ_PENALTY = 2; // penalty for reads with low MQ
  long constexpr NOT_FULLY_ALIGNED_READ_PENALTY = 3; // penalty for reads which were clipped
  long constexpr IS_READ_OVERLAPPING_PENALTY = 1; // penalty for reads that didn't fully overlap the variant breakpoint
  long constexpr IS_LOW_QUAL = 2; // penalty for reads that have low quality bases overlapping the variant

  // Make sure we do not underflow
  long epsilon_exponent_long = EPSILON_0_EXPONENT;

  // Mismatch penalty
  epsilon_exponent_long -= MISMATCH_PENALTY * mismatches;

  if (non_unique_paths)
    epsilon_exponent_long -= NON_UNIQUE_PATHS_PENALTY;

  if ((flags & IS_MAPQ_BAD) != 0u)
    epsilon_exponent_long -= BAD_MAPQ_PENALTY;

  // Decrease epsilon exponent if the path was not fully aligned, i.e. it has been clipped
  if (!fully_aligned)
    epsilon_exponent_long -= NOT_FULLY_ALIGNED_READ_PENALTY;

  if (!is_read_overlapping)
    epsilon_exponent_long -= IS_READ_OVERLAPPING_PENALTY;

  if (is_low_qual)
    epsilon_exponent_long -= IS_LOW_QUAL;

  // Make sure epsilon is not too low
  uint16_t epsilon_exponent = std::max(epsilon_exponent_long, static_cast<long>(8)) - 4;
  // the -4 is for historical reasons...

  // Get how many number of paths are in this haplotype
  uint32_t const cnum = get_genotype_num();

#ifndef NDEBUG
  // TODO(h2) THE FOLLOWING DOESNT MAKE SENSE
  /*assert(!explains.none());

  if (explains.none())
    explains.set(); // Flip 'em all, cause they can all explain this read
    */
#endif // NDEBUG

  //std::vector<uint16_t> haplotype_errors = find_with_how_many_errors_haplotypes_explain_the_read(cnum);

  // Update log scores
  assert(pn_index < hap_samples.size());
  auto & hap_sample = hap_samples[pn_index];

  // If we are generating statistics, we store the number of reads that support haplotype
#ifndef NDEBUG
  /*
  if (Options::instance()->stats.size() > 0)
  {
    assert(hap_sample.stats);

    if (explains.count() == 1)
    {
      // Update unique coverage
      auto it = std::find(haplotype_errors.begin(), haplotype_errors.end(), 0);
      assert(it != haplotype_errors.end());
      long const i = std::distance(haplotype_errors.begin(), it);
      assert(i >= 0);
      assert(i < static_cast<long>(hap_sample.stats->hap_unique_coverage.size()));

      if (hap_sample.stats->hap_unique_coverage[i] < 0xFFul)
        ++hap_sample.stats->hap_unique_coverage[i];

      assert(std::find(it + 1, haplotype_errors.end(), 0) == haplotype_errors.end());

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
  */
#endif // NDEBUG

  // Stop adding reads if the maximum log score is maxed out (this can happend if read depth is more than about 6000x)
  if (hap_sample.max_log_score < (0xFFFFul - epsilon_exponent))
  {
    hap_sample.max_log_score += epsilon_exponent;
    int i{0};

    for (std::size_t y = 0; y < cnum; ++y)
    {
      for (std::size_t x = 0; x <= y; ++x, ++i)
      {
        assert(to_index(x, y) < static_cast<long>(hap_sample.log_score.size()));
        assert(i == to_index(x, y));

        if (explains.contains(x) && explains.contains(y))
          hap_sample.log_score[i] += epsilon_exponent;
        else if (explains.contains(x) || explains.contains(y))
          hap_sample.log_score[i] += epsilon_exponent - 1;
        // else the log_score does not change
      }
    }
  }
}


void
Haplotype::update_max_log_score()
{
  for (auto & hap_sample : hap_samples)
  {
    assert(hap_sample.log_score.size() > 0);
    auto max_it = std::max_element(hap_sample.log_score.begin(), hap_sample.log_score.end());
    assert(max_it != hap_sample.log_score.end());
    hap_sample.max_log_score = *max_it;
    //calls.push_back(to_pair(std::distance(hap_sample.log_score.begin(), max_it)));
  }

  //assert(calls.size() == hap_samples.size());
}


/*
std::vector<uint32_t>
Haplotype::get_genotype_ids() const
{
  std::vector<uint32_t> gt_ids;

  for (unsigned i = 0; i < gts.size(); ++i)
    gt_ids.push_back(gts[i].id);

  return gt_ids;
}
*/

/*
std::vector<uint16_t>
Haplotype::get_haplotype_calls() const
{
  std::vector<uint16_t> all_calls{0u}; // Add the reference
  int const cnum = get_genotype_num();

  /// Filter options for haplotype extraction
  // the minimum support of any allele in a haplotype
  int const MINIMUM_VARIANT_SUPPORT = Options::const_instance()->minimum_extract_variant_support;

  // score*3 is approximately genotype quality compared to reference
  int const REQUIRED_SCORE_OVER_HOMREF = Options::const_instance()->minimum_extract_score_over_homref;
  ///

  for (auto const & hap_sample : hap_samples)
  {
    // Check if coverage is above the minimum variant support.
    // Note: The function may return false negatives.
    auto coverage_above_cutoff =
      [&](uint16_t call)
      {
        assert(call > 0);
        int q = cnum;
        assert(gts.size() == hap_sample.gt_coverage.size());

        for (int i = 0; i < static_cast<int>(hap_sample.gt_coverage.size()); ++i)
        {
          auto const & cov = hap_sample.gt_coverage[i];
          q /= gts[i].num;
          assert(q != 0);
          auto const a = call / q;
          assert(a < static_cast<int>(cov.size()));

          // Require 'MINIMUM_VARIANT_SUPPORT'-many reads to overlap at least one of the alleles of
          // the genotype call
          if (a > 0 && static_cast<int>(cov[a]) >= MINIMUM_VARIANT_SUPPORT)
            return true;

          call %= q;
        }

        return false;
      };


    // Favor homozygous calls
    int call1 = 0;
    int call2 = 0;

    int local_max_log_score = hap_sample.log_score[to_index(0, 0)] + REQUIRED_SCORE_OVER_HOMREF;

    // Homozugous first, because we want to have the least amount af haplotypes extracted
    for (int y = 1; y < cnum; ++y)
    {
      // Skip 0/0
      int const score = hap_sample.log_score[to_index(y, y)];

      if (score >= local_max_log_score)
      {
        // For homozygous calls, also check if there are enough read supporting the alternative allele
        local_max_log_score = score;
        call1 = y;
        call2 = y;
      }
    }

    //  heterzygous calls
    for (int y = 1; y < cnum; ++y)
    {
      for (int x = 0; x < y; ++x)
      {
        int const score = hap_sample.log_score[to_index(x, y)];

        if (score > local_max_log_score)
        {
          local_max_log_score = score;
          call1 = x;
          call2 = y;
        }
      }
    }

    // Add non-reference calls above coverage cut-off
    if (call1 != 0 && coverage_above_cutoff(call1))
      all_calls.push_back(call1);

    if (call2 == call1)
      all_calls.push_back(call2);
    else if (call2 != 0 && coverage_above_cutoff(call2))
      all_calls.push_back(call2);
  }

  return all_calls;
}
*/

} // namespace gyper
