#include <cassert> // assert()
#include <cmath> // sqrt
#include <cstdlib> // std::abs(int64_t)
#include <cstdint> // uint64_t
#include <map> // std::map<Key, Value>
#include <numeric> // std::accumulate
#include <string> // std::string
#include <sstream> // std::ostringstream
#include <vector> // std::vector

#include <boost/algorithm/string/split.hpp> // boost::split

#include <graphtyper/graph/absolute_position.hpp> // gyper::absolute_pos
#include <graphtyper/typer/var_stats.hpp> // gyper::VarStats
#include <graphtyper/typer/vcf.hpp> // gyper::get_all_pos(line, delim)
#include <graphtyper/utilities/options.hpp> // gyper::Options::instance()


namespace gyper
{

VarStats::VarStats(uint16_t allele_count) noexcept
  : mapq_allele_root_total(allele_count)
  , mapq_allele_counts(allele_count)
  , originally_clipped(allele_count)
  , realignment_distance(allele_count)
  , realignment_count(allele_count)
  , r1_strand_forward(allele_count)
  , r1_strand_reverse(allele_count)
  , r2_strand_forward(allele_count)
  , r2_strand_reverse(allele_count)
{}


void
VarStats::add_mapq(uint8_t allele_id, uint8_t const new_mapq)
{
  // Ignore MapQ == 255, that means MapQ is unavailable
  if (new_mapq < 255u)
  {
    ++mapq_count;

    if (new_mapq == 0u)
      ++mapq_zero_count;
    else
      mapq_root_total += static_cast<uint64_t>(new_mapq) * static_cast<uint64_t>(new_mapq); // Mapping quality squared

    // Per allele stats
    assert(allele_id < mapq_allele_counts.size());
    assert(allele_id < mapq_allele_root_total.size());

    ++mapq_allele_counts[allele_id];
    mapq_allele_root_total[allele_id] += static_cast<uint64_t>(new_mapq) * static_cast<uint64_t>(new_mapq);
  }
}


void
VarStats::add_realignment_distance(uint8_t const allele_id, uint32_t const original_pos, uint32_t const new_pos)
{
  assert (allele_id < realignment_distance.size());
  assert (allele_id < realignment_count.size());
  ++realignment_count[allele_id];
  uint32_t const distance = std::abs(static_cast<int64_t>(original_pos) - static_cast<int64_t>(new_pos));

  // Check for overflow
  if (static_cast<uint64_t>(realignment_distance[allele_id]) + static_cast<uint64_t>(distance) < 0xFFFFFFFFull)
    realignment_distance[allele_id] += distance;
  else
    realignment_distance[allele_id] = 0xFFFFFFFFull;
}

/**
 * CLASS INFORMATION
 */

uint8_t
VarStats::get_rms_mapq() const
{
  if (mapq_count > 0)
    return static_cast<uint64_t>(sqrt(static_cast<double>(mapq_root_total) / static_cast<double>(mapq_count)));
  else
    return 255u;
}


std::string
VarStats::get_rms_mapq_per_allele() const
{
  std::vector<uint64_t> rms_mapq_per_allele(mapq_allele_counts.size(), 255u);

  for (std::size_t i = 0; i < rms_mapq_per_allele.size(); ++i)
  {
    if (mapq_allele_counts[i] > 0)
    {
      rms_mapq_per_allele[i] =
        static_cast<uint64_t>(
          sqrt(static_cast<double>(mapq_allele_root_total[i]) / static_cast<double>(mapq_allele_counts[i]))
        );
    }
  }

  return join_strand_bias(rms_mapq_per_allele);
}


std::string
VarStats::get_originally_clipped_reads() const
{
  return join_strand_bias(this->originally_clipped);
}


std::string
VarStats::get_realignment_count() const
{
  return join_strand_bias(realignment_count); // TODO: Rename join_strand_bias()
}


std::string
VarStats::get_realignment_distance() const
{
  return join_strand_bias(realignment_distance); // TODO: Rename join_strand_bias()
}


std::string
VarStats::get_forward_strand_bias() const
{
  //return join_strand_bias(strand_forward);
  return join_strand_bias(r1_strand_forward, r2_strand_forward);
}


std::string
VarStats::get_reverse_strand_bias() const
{
  return join_strand_bias(r1_strand_reverse, r2_strand_reverse);
}


std::string
VarStats::get_unaligned_count() const
{
  std::ostringstream ss;
  ss << this->unaligned_reads;
  return ss.str();
}


std::string
VarStats::get_r1_forward_strand_bias() const
{
  return join_strand_bias(r1_strand_forward);
}


std::string
VarStats::get_r2_forward_strand_bias() const
{
  return join_strand_bias(r2_strand_forward);
}


std::string
VarStats::get_r1_reverse_strand_bias() const
{
  return join_strand_bias(r1_strand_reverse);
}


std::string
VarStats::get_r2_reverse_strand_bias() const
{
  return join_strand_bias(r2_strand_reverse);
}


/** Non-member functions */
template<class T>
std::string
join_strand_bias(std::vector<T> const & bias)
{
  if (bias.size() == 0)
    return std::string(".");

  std::stringstream ss;
  ss << bias[0];

  for (std::size_t i = 1; i < bias.size(); ++i)
    ss << ',' << bias[i];

  return ss.str();
}

// Explicit instantiation
template std::string join_strand_bias(std::vector<uint32_t> const & bias);
template std::string join_strand_bias(std::vector<uint64_t> const & bias);


std::string
join_strand_bias(std::vector<uint32_t> const & r1bias, std::vector<uint32_t> const & r2bias)
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


std::vector<std::string>
split_bias_to_strings(std::string const & bias)
{
  std::vector<std::string> splitted;
  boost::split(splitted, bias, [](char c){return c == ',';});
  return splitted;
}


std::vector<uint32_t>
split_bias_to_numbers(std::string const & bias)
{
  std::vector<std::string> splitted = split_bias_to_strings(bias);
  std::vector<uint32_t> nums;

  if (splitted.size() >= 1)
  {
    nums.reserve(splitted.size());

    for (auto const & spl : splitted)
      nums.push_back(std::strtoull(spl.c_str(), NULL, 10));
  }

  return nums;
}


std::vector<uint32_t>
get_strand_bias(std::map<std::string, std::string> const & infos, std::string const & bias)
{
  auto find_it = infos.find(bias);

  if (find_it != infos.end())
    return split_bias_to_numbers(find_it->second);
  else
    return std::vector<uint32_t>(0);
}


uint32_t
get_accumulated_strand_bias(std::map<std::string, std::string> const & infos, std::string const & bias)
{
  std::vector<uint32_t> strand_bias = get_strand_bias(infos, bias);
  return std::accumulate(strand_bias.begin(), strand_bias.end(), static_cast<uint32_t>(0u));
}


std::vector<uint16_t>
get_list_of_uncalled_alleles(std::string const & ac)
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


} // namespace gyper
