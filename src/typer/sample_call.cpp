#include <algorithm>
#include <iostream>
#include <numeric>

#include <graphtyper/graph/haplotype.hpp> // AlleleCoverage
#include <graphtyper/utilities/graph_help_functions.hpp> // to_index
#include <graphtyper/typer/sample_call.hpp>


namespace gyper
{

SampleCall::SampleCall() noexcept
   : ambiguous_depth(0)
{}

SampleCall::SampleCall(std::vector<uint8_t> const & _phred,
                       std::vector<uint16_t> const & _coverage,
                       uint8_t const _ambiguous_depth,
                       uint8_t const _ambiguous_depth_alt,
                       uint8_t const _alt_proper_pair_depth
  ) noexcept
  : phred(_phred)
  , coverage(_coverage)
  , ambiguous_depth(_ambiguous_depth)
  , alt_proper_pair_depth(_alt_proper_pair_depth)
{
  assert(coverage.size() > 1);
  assert(ambiguous_depth >= _ambiguous_depth_alt);

  // Only (ambiguous_depth - _ambiguous_depth_alt) covers the reference
  uint32_t const ref_depth = coverage[0] + ambiguous_depth - _ambiguous_depth_alt;

  assert(ref_depth >= coverage[0]);

  ref_total_depth = std::min(static_cast<uint32_t>(0xFFFFu),
                             static_cast<uint32_t>(ref_depth)
                             );

  // All ambiguous depth supports alt
  uint32_t const alt_depth = std::accumulate(coverage.begin() + 1, coverage.end(), 0u) + ambiguous_depth;

  alt_total_depth = std::min(static_cast<uint32_t>(0xFFFFu),
                             static_cast<uint32_t>(alt_depth)
                             );

  assert(alt_total_depth >= ambiguous_depth);
  assert(ref_total_depth >= coverage[0]);
}


std::pair<uint16_t, uint16_t>
SampleCall::get_gt_call() const
{
  std::size_t i = 0;

  for (std::size_t y = 0; y < coverage.size(); ++y)
  {
    for (std::size_t x = 0; x <= y; ++x, ++i)
    {
      if (phred[i] == 0)
      {
        return std::make_pair<uint16_t, uint16_t>(x, y);
      }
    }
  }

  assert(false);
  return std::make_pair<uint16_t, uint16_t>(0xFFFFu, 0xFFFFu);
}


uint8_t
SampleCall::get_gq() const
{
  bool seen_zero = false;
  uint8_t next_lowest_phred = 255;

  for (auto const p : phred)
  {
    if (p == 0)
    {
      if (!seen_zero)
        seen_zero = true;
      else
        return 0;
    }
    else
    {
      next_lowest_phred = std::min(next_lowest_phred, p);
    }
  }

  return next_lowest_phred;
}


void
SampleCall::change_to_ref_vs_all()
{
  assert(coverage.size() >= 2);
  assert(phred.size() >= 3);

  coverage.resize(2); // Change to 2 alleles
  phred.resize(3); // Change to 3 phred scores
  assert(alt_total_depth >= ambiguous_depth);
  assert(ref_total_depth >= coverage[0]);

  // Re-calculate ambiguous_depth_alt
  int32_t ambiguous_depth_alt = coverage[0] + ambiguous_depth - ref_total_depth;

  // Can happen if the depth overflew..
  ambiguous_depth_alt = std::min(static_cast<int32_t>(ambiguous_depth), ambiguous_depth_alt);

  // What was before ambiguous_depth_alt is no longer ambiguous
  assert(ambiguous_depth >= ambiguous_depth_alt);
  ambiguous_depth -= ambiguous_depth_alt;
  coverage[1] = alt_total_depth - ambiguous_depth;
  int32_t const alt_not_proper = coverage[1] > alt_proper_pair_depth ? coverage[1] - alt_proper_pair_depth : 0;
  int32_t const alt_proper = coverage[1] - alt_not_proper;

  // Only output 3 phred scores (0/0, 0/1 and 1/1)
  uint64_t constexpr ERROR_PHRED_PROPER = 25;
  uint64_t constexpr ERROR_PHRED_NOT_PROPER = 10;
  uint64_t const gt_00 = alt_proper * ERROR_PHRED_PROPER + alt_not_proper * ERROR_PHRED_NOT_PROPER;
  uint64_t const gt_01 = 3 * (coverage[0] + coverage[1]);
  uint64_t const gt_11 = coverage[0] * ERROR_PHRED_PROPER;

  uint64_t min_gt = std::min(gt_00, std::min(gt_01, gt_11));
  phred[0] = std::min(static_cast<uint64_t>(0xFFu), gt_00 - min_gt);
  phred[1] = std::min(static_cast<uint64_t>(0xFFu), gt_01 - min_gt);
  phred[2] = std::min(static_cast<uint64_t>(0xFFu), gt_11 - min_gt);
}


SampleCall
make_bi_allelic_call(SampleCall const & oc, long aa)
{
  // Old call should be multi-allelic
  assert(oc.coverage.size() > 2);
  assert(aa + 1 < static_cast<long>(oc.coverage.size()));

  SampleCall c;
  c.coverage.reserve(2);
  c.coverage.push_back(oc.coverage[0]);
  c.ambiguous_depth = oc.ambiguous_depth;
  c.ref_total_depth = oc.ref_total_depth;
  c.alt_total_depth = oc.alt_total_depth;
  c.alt_proper_pair_depth = oc.alt_proper_pair_depth;

  // Re-calculate ambiguous_depth_alt
  int32_t ambiguous_depth_alt = c.coverage[0] + c.ambiguous_depth - c.ref_total_depth;

  // Can happen if the depth overflew..
  ambiguous_depth_alt = std::min(static_cast<int32_t>(c.ambiguous_depth), ambiguous_depth_alt);

  // What was before ambiguous_depth_alt is no longer ambiguous
  assert(c.ambiguous_depth >= ambiguous_depth_alt);
  c.ambiguous_depth -= ambiguous_depth_alt;
  long cov_aa = c.alt_total_depth - c.ambiguous_depth;
  //c.coverage.push_back(c.alt_total_depth - c.ambiguous_depth);

  // Reduce coverage of reads that aligned uniquely to other alleles
  for (long a = 1; a < static_cast<long>(oc.coverage.size()); ++a)
  {
    if (a == aa + 1)
      continue; // Skip the target alternative allele

    cov_aa -= static_cast<long>(oc.coverage[a]);
  }

  c.coverage.push_back(std::max(cov_aa, 0l));

  int32_t const alt_not_proper = c.coverage[1] > c.alt_proper_pair_depth ?
                                 c.coverage[1] - c.alt_proper_pair_depth :
                                 0;

  int32_t const alt_proper = c.coverage[1] - alt_not_proper;

  c.phred.resize(3, 0); // 0/0, 0/1, 1/1
  assert(c.alt_total_depth >= c.ambiguous_depth);
  assert(c.ref_total_depth >= c.coverage[0]);

  // Only output 3 phred scores (0/0, 0/1 and 1/1)
  uint64_t constexpr ERROR_PHRED_PROPER = 25;
  uint64_t constexpr ERROR_PHRED_NOT_PROPER = 10;
  uint64_t const gt_00 = alt_proper * ERROR_PHRED_PROPER + alt_not_proper * ERROR_PHRED_NOT_PROPER;
  uint64_t const gt_01 = static_cast<uint64_t>(3 * (c.coverage[0] + c.coverage[1]));
  uint64_t const gt_11 = c.coverage[0] * ERROR_PHRED_PROPER;

  uint64_t min_gt = std::min(gt_00, std::min(gt_01, gt_11));
  c.phred[0] = std::min(static_cast<uint64_t>(0xFFu), gt_00 - min_gt);
  c.phred[1] = std::min(static_cast<uint64_t>(0xFFu), gt_01 - min_gt);
  c.phred[2] = std::min(static_cast<uint64_t>(0xFFu), gt_11 - min_gt);
  return c;
}


} // namespce gyper
