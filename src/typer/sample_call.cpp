#include <algorithm>
#include <iostream>
#include <numeric>

#include <cereal/archives/binary.hpp>
#include <graphtyper/utilities/logging.hpp>
#include <cereal/types/vector.hpp>

#include <graphtyper/graph/absolute_position.hpp>
#include <graphtyper/graph/haplotype.hpp> // AlleleCoverage
#include <graphtyper/graph/reference_depth.hpp>
#include <graphtyper/graph/sv.hpp> // SV
#include <graphtyper/utilities/graph_help_functions.hpp> // to_index
#include <graphtyper/typer/sample_call.hpp>


namespace
{

uint16_t
get_uint16(long const val)
{
  return static_cast<uint16_t>(
    std::max(static_cast<long>(std::numeric_limits<uint16_t>::min()),
             std::min(static_cast<long>(std::numeric_limits<uint16_t>::max()),
                      val
                      )
             )
    );
}


template <typename Tvec>
typename Tvec::value_type
median(Tvec & vec)
{
  std::nth_element(vec.begin(), vec.begin() + vec.size() / 2, vec.end());
  return vec[vec.size() / 2];
}


} // anon namespace


namespace gyper
{

SampleCall::SampleCall(std::vector<uint8_t> && _phred,
                       std::vector<uint16_t> && _coverage,
                       uint8_t const _ambiguous_depth,
                       uint8_t const _ambiguous_depth_alt,
                       uint8_t const _alt_proper_pair_depth) noexcept
  : phred(std::forward<std::vector<uint8_t> >(_phred))
  , coverage(std::forward<std::vector<uint16_t> >(_coverage))
  , ambiguous_depth(_ambiguous_depth)
  , alt_proper_pair_depth(_alt_proper_pair_depth)
{
  assert(coverage.size() > 1);
  assert(ambiguous_depth >= _ambiguous_depth_alt);

  // Only (ambiguous_depth - _ambiguous_depth_alt) covers the reference
  uint32_t const ref_depth = coverage[0] + ambiguous_depth - _ambiguous_depth_alt;

  assert(ref_depth >= coverage[0]);

  ref_total_depth = static_cast<uint16_t>(
    std::min(static_cast<uint32_t>(0xFFFFu),
             static_cast<uint32_t>(ref_depth)
             )
    );

  // All ambiguous depth supports alt
  uint32_t const alt_depth = std::accumulate(coverage.begin() + 1, coverage.end(), 0u) + ambiguous_depth;

  alt_total_depth = static_cast<uint16_t>(
    std::min(static_cast<uint32_t>(0xFFFFu),
             static_cast<uint32_t>(alt_depth))
    );

  assert(alt_total_depth >= ambiguous_depth);
  assert(ref_total_depth >= coverage[0]);
}


uint32_t
SampleCall::get_depth() const
{
  return std::accumulate(coverage.begin(), coverage.end(), static_cast<uint32_t>(ambiguous_depth));
}


uint32_t
SampleCall::get_unique_depth() const
{
  return std::accumulate(coverage.begin(), coverage.end(), static_cast<uint32_t>(0));
}


uint32_t
SampleCall::get_alt_depth() const
{
  return std::accumulate(coverage.begin() + 1, coverage.end(), static_cast<uint32_t>(ambiguous_depth));
}


std::pair<uint16_t, uint16_t>
SampleCall::get_gt_call() const
{
  if (phred.size() == 0)
    return std::make_pair<uint16_t, uint16_t>(0, 0);

  std::size_t i = 0;

  for (std::size_t y = 0; y < coverage.size(); ++y)
  {
    for (std::size_t x = 0; x <= y; ++x, ++i)
    {
      assert(i < phred.size());

      if (phred[i] == 0)
      {
        return std::make_pair<uint16_t, uint16_t>(static_cast<uint16_t>(x),
                                                  static_cast<uint16_t>(y)
                                                  );
      }
    }
  }

  print_log(log_severity::warning, __HERE__, " No phred==0");

  for (auto p : phred)
    std::cerr << static_cast<long>(p) << ",";

  std::cerr << std::endl;
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
    else if (p < next_lowest_phred)
    {
      next_lowest_phred = p;
    }
  }

  return next_lowest_phred;
}


uint8_t
SampleCall::get_lowest_phred_not_with(uint16_t allele) const
{
  long i{0};
  uint8_t min_phred{255};

  for (long y{0}; y < static_cast<long>(coverage.size()); ++y)
  {
    if (y == allele)
    {
      i += y + 1;
      continue;
    }

    for (long x = 0; x <= y; ++x, ++i)
    {
      if (x == allele)
        continue;

      if (phred[i] < min_phred)
        min_phred = phred[i];
    }
  }

  return min_phred;
}


int8_t
SampleCall::check_filter(long gq) const
{
  if (filter < 0)
  {
    if (gq > 40)
      filter = 0;
    else if (gq > 20)
      filter = 1;
    else if (gq > 0)
      filter = 2;
    else
      filter = 3;
  }

  return filter;
}


template <typename Archive>
void
SampleCall::serialize(Archive & ar, unsigned const int /*version*/)
{
  ar & phred;
  ar & coverage;
  ar & ref_total_depth;
  ar & alt_total_depth;
  ar & ambiguous_depth;
  ar & alt_proper_pair_depth;
}


template void SampleCall::serialize<cereal::BinaryInputArchive>(cereal::BinaryInputArchive &,
                                                                     const unsigned int);
template void SampleCall::serialize<cereal::BinaryOutputArchive>(cereal::BinaryOutputArchive &,
                                                                     const unsigned int);


SampleCall
make_bi_allelic_call(SampleCall const & oc, long aa)
{
  if (oc.coverage.size() == 2)
    return oc;

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
  int cov_aa = c.alt_total_depth - c.ambiguous_depth;

  // Reduce coverage of reads that aligned uniquely to other alleles
  for (long a = 1; a < static_cast<long>(oc.coverage.size()); ++a)
  {
    if (a == aa + 1)
      continue; // Skip the target alternative allele

    cov_aa -= static_cast<int>(oc.coverage[a]);
    c.alt_total_depth = static_cast<uint16_t>(
      std::max(0, static_cast<int>(c.alt_total_depth) - static_cast<int>(oc.coverage[a])));

    c.alt_proper_pair_depth = static_cast<uint8_t>(
      std::max(0, static_cast<int>(c.alt_proper_pair_depth) - static_cast<int>(oc.coverage[a])));
  }

  c.coverage.push_back(static_cast<unsigned short>(std::max(cov_aa, 0)));

  int32_t const alt_not_proper = c.coverage[1] > c.alt_proper_pair_depth ?
                                 c.coverage[1] - c.alt_proper_pair_depth :
                                 0;

  int32_t const alt_proper = c.coverage[1] - alt_not_proper;

  c.phred.resize(3, 0); // 0/0, 0/1, 1/1
  assert(c.alt_total_depth >= c.ambiguous_depth);
  assert(c.ref_total_depth >= c.coverage[0]);

  // Only output 3 phred scores (0/0, 0/1 and 1/1)
  uint64_t constexpr ERROR_PHRED_PROPER = 24;
  uint64_t constexpr ERROR_PHRED_NOT_PROPER = 12;
  uint64_t const gt_00 = alt_proper * ERROR_PHRED_PROPER + alt_not_proper * ERROR_PHRED_NOT_PROPER;
  uint64_t const gt_01 = static_cast<uint64_t>(3ul * (c.coverage[0] + static_cast<uint64_t>(c.coverage[1])));
  uint64_t const gt_11 = c.coverage[0] * ERROR_PHRED_PROPER;

  uint64_t min_gt = std::min(gt_00, std::min(gt_01, gt_11));
  c.phred[0] = static_cast<uint8_t>(std::min(static_cast<uint64_t>(0xFFu), gt_00 - min_gt));
  c.phred[1] = static_cast<uint8_t>(std::min(static_cast<uint64_t>(0xFFu), gt_01 - min_gt));
  c.phred[2] = static_cast<uint8_t>(std::min(static_cast<uint64_t>(0xFFu), gt_11 - min_gt));
  return c;
}


SampleCall
make_call_based_on_coverage(long pn_index, SV const & sv, ReferenceDepth const & reference_depth)
{
  SampleCall call;
  long abs_begin = absolute_pos.get_absolute_position(sv.chrom, sv.begin);
  long abs_end;

  if (sv.size < 190000)
    abs_end = absolute_pos.get_absolute_position(sv.chrom, sv.end);
  else
    abs_end = abs_begin + 190000;

  long constexpr N = 101; // Number of points, should be an odd number
  long constexpr M = 20; // Number of bp to jump outside of the deletion
  std::vector<uint16_t> depths_in;
  std::vector<uint16_t> depths_out;
  depths_in.reserve(N);

  // Inside the SV
  {
    long const size = abs_end - abs_begin;
    long N_in = std::min(N, size - 2 * M);

    if (N_in % 2 == 0)
      --N_in; // Make it an odd number

    for (long i = 1; i <= N_in; ++i)
    {
      auto abs_pos = (i * (size - 2 * M)) / (N_in + 1) + abs_begin + M;
      assert(abs_pos >= 0);
      depths_in.push_back(reference_depth.get_read_depth(abs_pos, pn_index));
    }
  }

  // Before the SV
  for (long i = 1; i <= N / 2 + 1; ++i)
  {
    auto abs_pos = std::max(abs_begin - i * M, 0l);
    depths_out.push_back(reference_depth.get_read_depth(abs_pos, pn_index));
  }

  // After of the SV (if the SV is not enormously huge)
  if (sv.size < 190000)
  {
    for (long i = 1; i <= N / 2; ++i)
    {
      auto abs_pos = std::max(abs_end + i * M, 0l);
      depths_out.push_back(reference_depth.get_read_depth(abs_pos, pn_index));
    }
  }

  long median_out = median(depths_out);
  long median_in = median(depths_in);

  // Only output 3 phred scores (0/0, 0/1 and 1/1)
  uint64_t ERROR = 12;
  uint64_t gt_00 = 0;
  uint64_t gt_01 = 0;
  uint64_t gt_11 = 0;

  if (sv.type == DEL || sv.type == DEL_ALU)
  {
    // Set "coverage"
    call.coverage.push_back(get_uint16(median_in));
    call.coverage.push_back(get_uint16(median_out - median_in));
  }
  else if (sv.type == DUP || sv.type == INV)
  {
    // Set "coverage"
    call.coverage.push_back(get_uint16(2 * median_out - std::max(median_in, median_out)));
    call.coverage.push_back(get_uint16(median_out - static_cast<long>(call.coverage[0])));
  }
  else
  {
    assert(false); // Only deletions and duplications expected in calling based on coverage
  }

  gt_00 = call.coverage[1] * ERROR;
  gt_01 = static_cast<uint64_t>(3 * (call.coverage[0] + call.coverage[1]));
  gt_11 = call.coverage[0] * ERROR;

  /// Normalize the genotype likelihoods such that at least one of them is 0
  {
    uint64_t const min_gt = std::min(gt_00, std::min(gt_01, gt_11));

    gt_00 -= min_gt;
    gt_01 -= min_gt;
    gt_11 -= min_gt;
  }

  // Coverage model should work more reliably for larger SVs
  if (sv.size <= 100)
  {
    gt_00 = (gt_00 * 2) / 3;
    gt_01 = (gt_01 * 2) / 3;
    gt_11 = (gt_11 * 2) / 3;
  }
  else if (sv.size > 1000)
  {
    gt_00 = (gt_00 * 3) / 2;
    gt_01 = (gt_01 * 3) / 2;
    gt_11 = (gt_11 * 3) / 2;
  }
  else if (sv.size > 10000)
  {
    gt_00 = gt_00 * 2;
    gt_01 = gt_01 * 2;
    gt_11 = gt_11 * 2;
  }

  // Set PHRED scores
  call.phred.resize(3, 0); // 0/0, 0/1, 1/1
  call.phred[0] = static_cast<uint8_t>(std::min(static_cast<uint64_t>(0xFFu), gt_00));
  call.phred[1] = static_cast<uint8_t>(std::min(static_cast<uint64_t>(0xFFu), gt_01));
  call.phred[2] = static_cast<uint8_t>(std::min(static_cast<uint64_t>(0xFFu), gt_11));
  return call;
}


} // namespce gyper
