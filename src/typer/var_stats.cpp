#include <cassert> // assert()
#include <cmath>   // sqrt
#include <cstdint> // uint64_t
#include <cstdlib> // std::abs(int64_t)
#include <map>     // std::map<Key, Value>
#include <numeric> // std::accumulate
#include <sstream> // std::ostringstream
#include <string>  // std::string
#include <vector>  // std::vector

#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>

#include <graphtyper/graph/absolute_position.hpp> // gyper::absolute_pos
#include <graphtyper/typer/var_stats.hpp>         // gyper::VarStats
#include <graphtyper/typer/vcf.hpp>               // gyper::get_all_pos(line, delim)
#include <graphtyper/utilities/options.hpp>       // gyper::Options::instance()
#include <graphtyper/utilities/string.hpp>

namespace gyper
{
VarStats::VarStats(std::size_t const allele_count) noexcept
  //  : realignment_distance(allele_count)
  //  , realignment_count(allele_count)
  :
  per_allele(allele_count), read_strand(allele_count)
{
}

/*
void
VarStats::add_realignment_distance(uint8_t const allele_id,
                                   uint32_t const original_pos,
                                   uint32_t const new_pos)
{
  assert(allele_id < realignment_distance.size());
  assert(allele_id < realignment_count.size());
  ++realignment_count[allele_id];
  uint32_t const distance = std::abs(static_cast<int64_t>(original_pos) - static_cast<int64_t>(new_pos));

  // Check for overflow
  if (static_cast<uint64_t>(realignment_distance[allele_id]) + static_cast<uint64_t>(distance) < 0xFFFFFFFFull)
    realignment_distance[allele_id] += distance;
  else
    realignment_distance[allele_id] = 0xFFFFFFFFull;
}
*/

/**
 * CLASS INFORMATION
 */

void VarStats::write_stats(std::map<std::string, std::string> & infos) const
{
  assert(per_allele.size() == read_strand.size());
  long const num_allele = per_allele.size();

  if (num_allele <= 1)
    return;

  infos["CR"] = std::to_string(clipped_reads);
  infos["MQsquared"] = std::to_string(mapq_squared);
  write_read_strand_stats(infos);
  write_per_allele_stats(infos);
}

void VarStats::write_per_allele_stats(std::map<std::string, std::string> & infos) const
{
  long const num_allele = per_allele.size();
  assert(num_allele > 1);

  std::ostringstream cr_ss;
  std::ostringstream mq_ss;
  std::ostringstream sd_ss;
  std::ostringstream mm_ss;

  {
    VarStatsPerAllele const & stats_per_al = per_allele[0];
    cr_ss << stats_per_al.clipped_bp;
    mq_ss << stats_per_al.mapq_squared;
    sd_ss << stats_per_al.score_diff;
    mm_ss << stats_per_al.mismatches;
  }

  for (long i{1}; i < num_allele; ++i)
  {
    VarStatsPerAllele const & stats_per_al = per_allele[i];
    cr_ss << ',' << stats_per_al.clipped_bp;
    mq_ss << ',' << stats_per_al.mapq_squared;
    sd_ss << ',' << stats_per_al.score_diff;
    mm_ss << ',' << stats_per_al.mismatches;
  }

  infos["CRal"] = cr_ss.str();
  infos["MQSal"] = mq_ss.str();
  infos["SDal"] = sd_ss.str();
  infos["MMal"] = mm_ss.str();
}

void VarStats::write_read_strand_stats(std::map<std::string, std::string> & infos) const
{
  long const num_allele = read_strand.size();
  assert(num_allele > 1);

  std::ostringstream f_ss;
  std::ostringstream r_ss;
  std::ostringstream f1_ss;
  std::ostringstream f2_ss;
  std::ostringstream r1_ss;
  std::ostringstream r2_ss;

  {
    auto const & stats = read_strand[0];
    f_ss << (stats.r1_forward + stats.r2_forward);
    r_ss << (stats.r1_reverse + stats.r2_reverse);
    f1_ss << stats.r1_forward;
    f2_ss << stats.r2_forward;
    r1_ss << stats.r1_reverse;
    r2_ss << stats.r2_reverse;
  }

  for (long i{1}; i < num_allele; ++i)
  {
    auto const & stats = read_strand[i];
    f_ss << ',' << (stats.r1_forward + stats.r2_forward);
    r_ss << ',' << (stats.r1_reverse + stats.r2_reverse);
    f1_ss << ',' << stats.r1_forward;
    f2_ss << ',' << stats.r2_forward;
    r1_ss << ',' << stats.r1_reverse;
    r2_ss << ',' << stats.r2_reverse;
  }

  infos["SBF"] = f_ss.str();
  infos["SBR"] = r_ss.str();
  infos["SBF1"] = f1_ss.str();
  infos["SBF2"] = f2_ss.str();
  infos["SBR1"] = r1_ss.str();
  infos["SBR2"] = r2_ss.str();
}

void VarStats::add_stats(VarStats const & stats)
{
  assert(per_allele.size() == stats.per_allele.size());
  assert(read_strand.size() == stats.read_strand.size());

  clipped_reads += stats.clipped_reads;
  mapq_squared += stats.mapq_squared;
  n_genotyped += stats.n_genotyped;
  n_calls += stats.n_calls;
  n_passed_calls += stats.n_passed_calls;
  n_ref_ref += stats.n_ref_ref;
  n_ref_alt += stats.n_ref_alt;
  n_alt_alt += stats.n_alt_alt;
  n_max_alt_proper_pairs += stats.n_max_alt_proper_pairs;
  het_allele_depth.first += stats.het_allele_depth.first;
  het_allele_depth.second += stats.het_allele_depth.second;
  hom_allele_depth.first += stats.hom_allele_depth.first;
  hom_allele_depth.second += stats.hom_allele_depth.second;
  seqdepth += stats.seqdepth;

  for (long i{0}; i < static_cast<long>(per_allele.size()); ++i)
  {
    auto & new_a = per_allele[i];
    auto & old_a = stats.per_allele[i];
    auto & new_rs = read_strand[i];
    auto & old_rs = stats.read_strand[i];

    new_a.clipped_bp += old_a.clipped_bp;
    new_a.mapq_squared += old_a.mapq_squared;
    new_a.score_diff += old_a.score_diff;
    new_a.mismatches += old_a.mismatches;
    new_a.qd_qual += old_a.qd_qual;
    new_a.qd_depth += old_a.qd_depth;
    new_a.total_depth += old_a.total_depth;
    new_a.ac += old_a.ac;
    new_a.pass_ac += old_a.pass_ac;
    new_a.maximum_alt_support = std::max(new_a.maximum_alt_support, old_a.maximum_alt_support);
    new_a.maximum_alt_support_ratio = std::max(new_a.maximum_alt_support_ratio, old_a.maximum_alt_support_ratio);
    new_a.het_multi_allele_depth.first += old_a.het_multi_allele_depth.first;
    new_a.het_multi_allele_depth.second += old_a.het_multi_allele_depth.second;
    new_a.hom_multi_allele_depth.first += old_a.hom_multi_allele_depth.first;
    new_a.hom_multi_allele_depth.second += old_a.hom_multi_allele_depth.second;

    new_rs.r1_forward += old_rs.r1_forward;
    new_rs.r1_reverse += old_rs.r1_reverse;
    new_rs.r2_forward += old_rs.r2_forward;
    new_rs.r2_reverse += old_rs.r2_reverse;
  }
}

void VarStats::read_stats(std::map<std::string, std::string> const & infos)
{
  {
    auto find_it = infos.find("CR");

    if (find_it != infos.end())
      clipped_reads += std::strtoull(find_it->second.c_str(), NULL, 10);
  }

  {
    auto find_it = infos.find("MQsquared");

    if (find_it != infos.end())
      mapq_squared += std::strtoull(find_it->second.c_str(), NULL, 10);
  }

  {
    auto find_it = infos.find("SBF1");

    if (find_it != infos.end())
    {
      std::stringstream ss(find_it->second);
      long i{0};

      for (uint32_t num; ss >> num; ++i)
      {
        assert(i < static_cast<long>(read_strand.size()));
        read_strand[i].r1_forward += num;

        if (ss.peek() == ',')
          ss.ignore();
      }
    }
  }

  {
    auto find_it = infos.find("SBF2");

    if (find_it != infos.end())
    {
      std::stringstream ss(find_it->second);
      long i{0};

      for (uint32_t num; ss >> num; ++i)
      {
        assert(i < static_cast<long>(read_strand.size()));
        read_strand[i].r2_forward += num;

        if (ss.peek() == ',')
          ss.ignore();
      }
    }
  }

  {
    auto find_it = infos.find("SBR1");

    if (find_it != infos.end())
    {
      std::stringstream ss(find_it->second);
      long i{0};

      for (uint32_t num; ss >> num; ++i)
      {
        assert(i < static_cast<long>(read_strand.size()));
        read_strand[i].r1_reverse += num;

        if (ss.peek() == ',')
          ss.ignore();
      }
    }
  }

  {
    auto find_it = infos.find("SBR2");

    if (find_it != infos.end())
    {
      std::stringstream ss(find_it->second);
      long i{0};

      for (uint32_t num; ss >> num; ++i)
      {
        assert(i < static_cast<long>(read_strand.size()));
        read_strand[i].r2_reverse += num;

        if (ss.peek() == ',')
          ss.ignore();
      }
    }
  }

  {
    auto find_it = infos.find("CRal");

    if (find_it != infos.end())
    {
      std::stringstream ss(find_it->second);
      long i{0};

      for (uint32_t num; ss >> num; ++i)
      {
        assert(i < static_cast<long>(per_allele.size()));
        per_allele[i].clipped_bp += num;

        if (ss.peek() == ',')
          ss.ignore();
      }
    }
  }

  {
    auto find_it = infos.find("MQSal");

    if (find_it != infos.end())
    {
      std::stringstream ss(find_it->second);
      long i{0};

      for (uint32_t num; ss >> num; ++i)
      {
        assert(i < static_cast<long>(per_allele.size()));
        per_allele[i].mapq_squared += num;

        if (ss.peek() == ',')
          ss.ignore();
      }
    }
  }

  {
    auto find_it = infos.find("SDal");

    if (find_it != infos.end())
    {
      std::stringstream ss(find_it->second);
      long i{0};

      for (uint32_t num; ss >> num; ++i)
      {
        assert(i < static_cast<long>(per_allele.size()));
        per_allele[i].score_diff += num;

        if (ss.peek() == ',')
          ss.ignore();
      }
    }
  }

  {
    auto find_it = infos.find("MMal");

    if (find_it != infos.end())
    {
      std::stringstream ss(find_it->second);
      long i{0};

      for (uint32_t num; ss >> num; ++i)
      {
        assert(i < static_cast<long>(per_allele.size()));
        per_allele[i].mismatches += num;

        if (ss.peek() == ',')
          ss.ignore();
      }
    }
  }
}

/** Non-member functions */
template <class T>
std::string join_strand_bias(std::vector<T> const & bias)
{
  if (bias.size() == 0)
    return std::string(".");

  std::stringstream ss;
  ss << bias[0];

  for (long i = 1; i < static_cast<long>(bias.size()); ++i)
    ss << ',' << bias[i];

  return ss.str();
}

// Explicit instantiation
template std::string join_strand_bias(std::vector<uint32_t> const & bias);
template std::string join_strand_bias(std::vector<uint64_t> const & bias);

std::string join_strand_bias(std::vector<uint32_t> const & r1bias, std::vector<uint32_t> const & r2bias)
{
  assert(r1bias.size() == r2bias.size());

  if (r1bias.size() == 0)
    return std::string(".");

  std::stringstream ss;
  ss << (r1bias[0] + r2bias[0]);

  for (std::size_t i = 1; i < r1bias.size(); ++i)
    ss << ',' << (r1bias[i] + r2bias[i]);

  return ss.str();
}

std::vector<uint32_t> split_bias_to_numbers(std::string const & bias)
{
  auto splitted = split_on_delim(bias, ',');
  std::vector<uint32_t> nums;

  if (splitted.size() >= 1)
  {
    nums.reserve(splitted.size());

    for (auto const & spl : splitted)
      nums.push_back(stoi64(spl));
  }

  return nums;
}

std::vector<uint32_t> get_strand_bias(std::map<std::string, std::string> const & infos, std::string const & bias)
{
  auto find_it = infos.find(bias);

  if (find_it != infos.end())
    return split_bias_to_numbers(find_it->second);
  else
    return std::vector<uint32_t>(0);
}

long get_accumulated_strand_bias(std::map<std::string, std::string> const & infos, std::string const & bias)
{
  std::vector<uint32_t> strand_bias = get_strand_bias(infos, bias);
  return std::accumulate(strand_bias.begin(), strand_bias.end(), 0l);
}

long get_accumulated_alt_strand_bias(std::map<std::string, std::string> const & infos, std::string const & bias)
{
  std::vector<uint32_t> strand_bias = get_strand_bias(infos, bias);

  if (strand_bias.size() == 0)
    return 0;

  return std::accumulate(strand_bias.begin() + 1, strand_bias.end(), 0l);
}

std::vector<uint16_t> get_list_of_uncalled_alleles(std::string const & ac)
{
  std::vector<uint16_t> uncalled_alleles;

  {
    std::vector<std::size_t> all_ac_commas = get_all_pos(ac, ',');
    std::string zero("0");

    for (std::size_t i = 0; i < all_ac_commas.size() - 1; ++i)
    {
      if (get_string_at_tab_index(ac, all_ac_commas, i) == zero)
        uncalled_alleles.push_back(i + 1); // Reference does not have AC, and we never want to delete it anyway
    }
  }

  return uncalled_alleles;
}

template <typename Archive>
void VarStats::serialize(Archive & ar, unsigned const int /*version*/)
{
  ar & per_allele;
  ar & read_strand;
  ar & clipped_reads;
  ar & mapq_squared;
  ar & n_genotyped;
  ar & n_calls;
  ar & n_passed_calls;
  ar & n_ref_ref;
  ar & n_ref_alt;
  ar & n_alt_alt;
  ar & n_max_alt_proper_pairs;
  ar & seqdepth;
  ar & het_allele_depth;
  ar & hom_allele_depth;
}

template void VarStats::serialize<cereal::BinaryInputArchive>(cereal::BinaryInputArchive &, const unsigned int);
template void VarStats::serialize<cereal::BinaryOutputArchive>(cereal::BinaryOutputArchive &, const unsigned int);

} // namespace gyper
