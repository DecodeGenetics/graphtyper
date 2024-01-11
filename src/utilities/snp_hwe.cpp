#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>

#include <graphtyper/utilities/logging.hpp>
#include <graphtyper/utilities/snp_hwe.hpp>

/*
// Source: http://csg.sph.umich.edu/abecasis/Exact/snp_hwe.c
//
// This code implements an exact SNP test of Hardy-Weinberg Equilibrium as described in
// Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of
// Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76: 000 - 000
//
// Written by Jan Wigginton
*/
namespace gyper
{
double p_hwe_excess_het(int obs_hets, int obs_hom1, int obs_hom2)
{
  if (obs_hom1 < 0 || obs_hom2 < 0 || obs_hets < 0)
  {
    print_error(
      "FATAL ERROR - SNP-HWE: Current genotype configuration (%d  %d %d ) includes a"
      " negative count",
      obs_hets,
      " ",
      obs_hom1,
      " ",
      obs_hom2);

    exit(EXIT_FAILURE);
  }

  if (obs_hets == 0 && (obs_hom1 == 0 || obs_hom2 == 0))
    return 1.0;

  int obs_homc = obs_hom1 < obs_hom2 ? obs_hom2 : obs_hom1;
  int obs_homr = obs_hom1 < obs_hom2 ? obs_hom1 : obs_hom2;

  int rare_copies = 2 * obs_homr + obs_hets;
  int genotypes = obs_hets + obs_homc + obs_homr;

  std::vector<double> het_probs(rare_copies + 1, 0.0);

  // start at midpoint
  double mid_double = static_cast<double>(rare_copies) * static_cast<double>(2 * genotypes - rare_copies) /
                      static_cast<double>(2 * genotypes);

  assert(mid_double < static_cast<double>(std::numeric_limits<int>::max()));
  int mid = std::lround(mid_double);

  // check to ensure that midpoint and rare alleles have same parity
  if ((rare_copies & 1) ^ (mid & 1))
    mid++;

  int curr_hets = mid;
  int curr_homr = (rare_copies - mid) / 2;
  int curr_homc = genotypes - curr_hets - curr_homr;

  assert(mid < static_cast<int>(het_probs.size()));
  het_probs[mid] = 1.0;
  double sum = het_probs[mid];

  for (/*curr_hets = mid*/; curr_hets > 1; curr_hets -= 2)
  {
    assert(curr_hets < static_cast<int>(het_probs.size()));
    assert(curr_hets > 1);

    het_probs[curr_hets - 2] =
      het_probs[curr_hets] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));

    sum += het_probs[curr_hets - 2];

    // 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote
    curr_homr++;
    curr_homc++;
  }

  curr_hets = mid;
  curr_homr = (rare_copies - mid) / 2;
  curr_homc = genotypes - curr_hets - curr_homr;

  for (curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2)
  {
    assert(curr_hets < static_cast<int>(het_probs.size()));
    assert((curr_hets + 2) < static_cast<int>(het_probs.size()));

    het_probs[curr_hets + 2] =
      het_probs[curr_hets] * 4.0 * curr_homr * curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0));

    sum += het_probs[curr_hets + 2];

    // add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote
    curr_homr--;
    curr_homc--;
  }

  for (int i = 0; i <= rare_copies; i++)
    het_probs[i] /= sum;

  // alternate p-value calculation for p_hi
  assert(obs_hets < static_cast<int>(het_probs.size()));
  double p_hi = 0.0;

  for (int i = obs_hets; i <= rare_copies; i++)
    p_hi += het_probs[i];

  return p_hi > 1.0 ? 1.0 : p_hi;
}

} // namespace gyper
