#pragma once

#include <array>
#include <cmath>

namespace
{

double constexpr LOGF_INTERCEPT{-29.28908};
double constexpr LOGF_ABHom{23.12909};
double constexpr LOGF_CRBySeqDepth{-10.22658};
double constexpr LOGF_MQ{0.01024};
double constexpr LOGF_PASS_ratio{0.85320};
double constexpr LOGF_GTYield{4.91178};
double constexpr LOGF_QD{0.23215};

std::array<double, 11l> constexpr LOGF_ABHet{{-6.03446,
  -6.03446,
  -1.35948,
  -0.84956,
  -0.28956,
  0.0,
  -1.05013,
  -1.35024,
  -1.34475,
  -3.74512,
  -3.74512}};

std::array<double, 11l> constexpr LOGF_SBAlt{{-0.32486,
  -0.32486,
  -0.25342,
  -0.32696,
  0.02442,
  0.0,
  -0.33522,
  -0.41332,
  -0.74043,
  -1.60844,
  -1.60844
}};

} // anon namespace


namespace gyper
{

inline double
get_logf(double ab_hom,
         double cr_by_seqdepth,
         double mq,
         double pass_ratio,
         double gt_yield,
         double qd,
         long ab_het_bin,
         long sbalt_bin)
{
  assert(ab_het_bin <= 10);
  assert(ab_het_bin >= 0);
  assert(sbalt_bin <= 10);
  assert(sbalt_bin >= 0);

  //BOOST_LOG_TRIVIAL(info) << __HERE__ << " " << ab_hom << " " << cr_by_seqdepth << " " << mq << " " << pass_ratio
  //                        << " " << gt_yield << " " << qd << " " << ab_het_bin << " " << sbalt_bin;

  double const pwr = LOGF_INTERCEPT + ab_hom * LOGF_ABHom + cr_by_seqdepth * LOGF_CRBySeqDepth +
                     mq * LOGF_MQ + pass_ratio * LOGF_PASS_ratio + gt_yield * LOGF_GTYield + qd * LOGF_QD +
                     LOGF_ABHet[ab_het_bin] + LOGF_SBAlt[sbalt_bin];
  //BOOST_LOG_TRIVIAL(info) << __HERE__ << " log model pwr = " << pwr;
  return 1.0 / (1.0 + std::exp(-pwr));
}


} // namespace gyper
