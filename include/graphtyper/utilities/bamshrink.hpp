#pragma once

namespace bamshrink
{

struct Options
{
  int maxFragLen = 1000;
  int minNumMatching = 40;
  double avgCovByReadLen = 0.30000001;
  std::string bamPathIn;
  std::string bamIndex = "<bamPathIn>.[bai,crai]";
  std::string bamPathOut = "-";
  std::string interval;
  std::string intervalFile;
  long SUPER_HI_DEPTH = 2;
};

int main(Options & opts);

} // namespace bamshrink


namespace gyper
{

void
bamshrink(std::string const & chrom,
          int begin,
          int end,
          std::string const & path_in,
          std::string const & path_out,
          double const avg_cov_by_readlen,
          std::string const & ref_fn);


void
bamshrink_multi(std::string const & interval_fn,
                std::string const & path_in,
                std::string const & path_out,
                double const avg_cov_by_readlen,
                std::string const & ref_fn);

} // namespace gyper
