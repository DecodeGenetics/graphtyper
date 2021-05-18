#pragma once

#include <string>

namespace bamshrink
{

struct Options
{
  // Note: You must change defaults in gyper options as well
  int maxFragLen{1000};
  int minNumMatching{55};
  bool is_filtering_mapq0{true};
  bool no_filter_on_coverage{false};
  int minReadLen{75};
  int minReadLenMapQ0{94};
  int minUnpairedReadLen{94};
  long as_filter_threshold{40};

  double avgCovByReadLen{0.30000001};
  std::string bamPathIn{};
  std::string bamIndex{"<bamPathIn>.[bai,crai]"};
  std::string bamPathOut{"-"};
  std::string interval{};
  std::string intervalFile{};
  long SUPER_HI_DEPTH{2};
};

int main(Options & opts);

} // namespace bamshrink


namespace gyper
{

void
bamshrink(std::string const chrom,
          int begin,
          int end,
          std::string const path_in,
          std::string const sam_index_in,
          std::string const path_out,
          double const avg_cov_by_readlen,
          std::string const ref_fn);


void
bamshrink_multi(std::string const interval_fn,
                std::string const path_in,
                std::string const path_out,
                double const avg_cov_by_readlen,
                std::string const ref_fn);

} // namespace gyper
