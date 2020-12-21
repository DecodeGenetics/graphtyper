#pragma once

#include <array>
#include <cmath>

namespace
{

// logf
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
get_logf(double abhom,
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

  double const pwr = LOGF_INTERCEPT + abhom * LOGF_ABHom + cr_by_seqdepth * LOGF_CRBySeqDepth +
                     mq * LOGF_MQ + pass_ratio * LOGF_PASS_ratio + gt_yield * LOGF_GTYield + qd * LOGF_QD +
                     LOGF_ABHet[ab_het_bin] + LOGF_SBAlt[sbalt_bin];
  //BOOST_LOG_TRIVIAL(info) << __HERE__ << " log model pwr = " << pwr;
  double const _exp = std::max(0.0, std::exp(-pwr));
  return 1.0 / (1.0 + _exp);
}


inline double
get_aa_score(double abhom, double sb, double mm, long sd, double qd, double cr, long mq)
{
  double constexpr AA_INTERCEPT{-6.347426707};
  double constexpr AA_SB{-0.25233400};
  double constexpr AA_MM{-0.04129973};
  double constexpr AA_SD{0.014572295};
  double constexpr AA_QD{0.065221319};
  double constexpr AA_CR{-0.01934834};
  double constexpr AA_MQ{0.055973424};

  std::array<double, 5l> constexpr AA_ABHom{{0.0,
    1.304140117,
    1.681221065,
    2.214801195,
    3.930106559}};

  long abhom_bin{4};

  if (abhom <= 0.85)
    abhom_bin = 0;
  else if (abhom <= 0.94)
    abhom_bin = 1;
  else if (abhom <= 0.98)
    abhom_bin = 2;
  else if (abhom <= 0.99)
    abhom_bin = 3;

  if (mq > 60)
    mq = 60;

  double const pwr = AA_INTERCEPT +
                     AA_ABHom[abhom_bin] +
                     sb * AA_SB +
                     mm * AA_MM +
                     sd * AA_SD +
                     qd * AA_QD +
                     cr * AA_CR +
                     mq * AA_MQ;

  /* TODO: check if this is ok
  if (pwr < 0)
  {
    double const _exp = std::exp(pwr);
    return _exp / (1.0 + _exp);
  }
  */

  double const _exp = std::exp(-pwr);
  assert(_exp >= 0.0);
  return 1.0 / (1.0 + _exp);
}


} // namespace gyper
