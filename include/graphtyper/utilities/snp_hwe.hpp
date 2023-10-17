#pragma once

namespace gyper
{
/*
// Source: http://csg.sph.umich.edu/abecasis/Exact/snp_hwe.c
//
// This code implements an exact SNP test of Hardy-Weinberg Equilibrium as described in
// Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of
// Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76: 000 - 000
//
// Written by Jan Wigginton
*/
double p_hwe_excess_het(int obs_hets, int obs_hom1, int obs_hom2);

} // namespace gyper
