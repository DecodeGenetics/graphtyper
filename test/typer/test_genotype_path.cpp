#include <catch.hpp>

#include <stdio.h>
#include <climits>
#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>

#include <graphtyper/graph/graph.hpp>
#include <graphtyper/index/kmer_label.hpp>
#include <graphtyper/utilities/type_conversions.hpp> // to_uint64()
#include <graphtyper/typer/path.hpp>
#include <graphtyper/typer/genotype_paths.hpp>


/*
TEST_CASE("Genotype paths")
{
  seqan::IupacString const test_read = "ACGTTGCA";

  gyper::KmerLabel label1(13, 44, 0);
  label1.variant_order = 39;
  label1.variant_num = 0;

  gyper::KmerLabel label2(13, 44, 1);
  label2.variant_order = 16;
  label2.variant_num = 0;

  gyper::GenotypePaths geno(0, seqan::length(test_read));
  geno.add_next_kmer_labels({label1, label2}, 0, 31);

  REQUIRE(geno.paths.size() == 1);
  REQUIRE(geno.paths[0].size() == 32);
  REQUIRE(geno.paths[0].start == 13);
  REQUIRE(geno.paths[0].end == 44);
  REQUIRE(geno.paths[0].var_order.size() == 2);
  REQUIRE(geno.paths[0].var_order[0] == 39);
  REQUIRE(geno.paths[0].var_order[1] == 16);
  REQUIRE(geno.paths[0].nums.size() == 2);
  REQUIRE(geno.paths[0].nums[0].to_ulong() == 1);
  REQUIRE(geno.paths[0].nums[1].to_ulong() == 1);

  gyper::KmerLabel label3(44, 75, 3);
  label3.variant_order = 69;
  label3.variant_num = 1;

  geno.add_next_kmer_labels({label3}, 31, 62, 0); // mismatches

  REQUIRE(geno.paths.size() == 1);
  REQUIRE(geno.paths[0].size() == 63);
  REQUIRE(geno.paths[0].start == 13);
  REQUIRE(geno.paths[0].end == 75);
  REQUIRE(geno.paths[0].read_start_index == 0);
  REQUIRE(geno.paths[0].var_order.size() == 3);
  REQUIRE(geno.paths[0].var_order[0] == 69);
  REQUIRE(geno.paths[0].var_order[1] == 39);
  REQUIRE(geno.paths[0].var_order[2] == 16);
  REQUIRE(geno.paths[0].nums.size() == 3);
  REQUIRE(geno.paths[0].nums[0].to_ulong() == (1 << 1));
  REQUIRE(geno.paths[0].nums[1].to_ulong() == (1 << 0));
  REQUIRE(geno.paths[0].nums[2].to_ulong() == (1 << 0));

  geno.add_next_kmer_labels({gyper::KmerLabel(75, 167, 4, 10, 136),
                             gyper::KmerLabel(75, 167, 5, 0, 121),
                             gyper::KmerLabel(75, 137, 6, 0, 136),
                             gyper::KmerLabel(75, 137, 7, 0, 121)}, 62, 93
                            );

  REQUIRE(geno.paths.size() == 2);
  REQUIRE(geno.paths[0].size() == 94);
  REQUIRE(geno.paths[0].start == 13);
  REQUIRE(geno.paths[0].end == 167);
  REQUIRE(geno.paths[0].read_start_index == 0);
  REQUIRE(geno.paths[0].var_order.size() == 5);
  REQUIRE(geno.paths[0].var_order[0] == 136);
  REQUIRE(geno.paths[0].var_order[1] == 121);
  REQUIRE(geno.paths[0].var_order[2] == 69);
  REQUIRE(geno.paths[0].var_order[3] == 39);
  REQUIRE(geno.paths[0].var_order[4] == 16);
  REQUIRE(geno.paths[0].nums.size() == 5);
  REQUIRE(geno.paths[0].nums[0].to_ulong() == (1 << 10));
  REQUIRE(geno.paths[0].nums[1].to_ulong() == (1 << 0));
  REQUIRE(geno.paths[0].nums[2].to_ulong() == (1 << 1));
  REQUIRE(geno.paths[0].nums[3].to_ulong() == (1 << 0));
  REQUIRE(geno.paths[0].nums[4].to_ulong() == (1 << 0));
  REQUIRE(geno.paths[1].size() == 94);
  REQUIRE(geno.paths[1].start == 13);
  REQUIRE(geno.paths[1].end_pos() == 137);
  REQUIRE(geno.paths[1].var_order.size() == 5);
  REQUIRE(geno.paths[1].var_order[0] == 136);
  REQUIRE(geno.paths[1].var_order[1] == 121);
  REQUIRE(geno.paths[1].var_order[2] == 69);
  REQUIRE(geno.paths[1].var_order[3] == 39);
  REQUIRE(geno.paths[1].var_order[4] == 16);
  REQUIRE(geno.paths[1].nums.size() == 5);
  REQUIRE(geno.paths[1].nums[0].to_ulong() == (1 << 0));
  REQUIRE(geno.paths[1].nums[1].to_ulong() == (1 << 0));
  REQUIRE(geno.paths[1].nums[2].to_ulong() == (1 << 1));
  REQUIRE(geno.paths[1].nums[3].to_ulong() == (1 << 0));
  REQUIRE(geno.paths[1].nums[4].to_ulong() == (1 << 0));

  // Kmer 0: CCCACCCTCCTGGGCTTCCCCTTGGCACTCCG
  // 43129113 43129144 43129139 0
  // 43129113 43129144 43129116 0
  // Kmer 1: GCTGTCACCCGTCTGGCCCCATTGCTGGGGCC
  // 43129144 43129175 43129169 1
  // Kmer 2: CTGCCCCGAGCCTGACTCCTCTGCTTTGCTCC
  // 43129175 43129206 43129192 0
  // 43129175 43129206 43129189 0
  // Kmer 3: CCACAGGTGTCCTGCGAACAGGTGCTGCTGGC
  // 43129206 43129267 43129236 36
  // 43129206 43129267 43129221 0
  // 43129206 43129237 43129236 35
  // 43129206 43129237 43129221 0
}

*/
