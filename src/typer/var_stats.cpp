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
//  : realignment_distance(allele_count)
//  , realignment_count(allele_count)
  : read_strand(allele_count)
//  , r1_strand_forward(allele_count)
//  , r1_strand_reverse(allele_count)
//  , r2_strand_forward(allele_count)
//  , r2_strand_reverse(allele_count)
{}


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

//std::string
//VarStats::get_realignment_count() const
//{
//  return join_strand_bias(realignment_count); // TODO: Rename join_strand_bias()
//}
//
//
//std::string
//VarStats::get_realignment_distance() const
//{
//  return join_strand_bias(realignment_distance); // TODO: Rename join_strand_bias()
//}


std::string
VarStats::get_forward_strand_bias() const
{
  if (read_strand.size() == 0)
    return ".";

  std::stringstream ss;
  ss << (read_strand[0].r1_forward + read_strand[0].r2_forward);

  for (long i = 1; i < static_cast<long>(read_strand.size()); ++i)
    ss << ',' << (read_strand[i].r1_forward + read_strand[i].r2_forward);

  return ss.str();
}


std::string
VarStats::get_reverse_strand_bias() const
{
  if (read_strand.size() == 0)
    return ".";

  std::stringstream ss;
  ss << (read_strand[0].r1_reverse + read_strand[0].r2_reverse);

  for (long i = 1; i < static_cast<long>(read_strand.size()); ++i)
    ss << ',' << (read_strand[i].r1_reverse + read_strand[i].r2_reverse);

  return ss.str();
}


std::string
VarStats::get_r1_forward_strand_bias() const
{
  if (read_strand.size() == 0)
    return ".";

  std::stringstream ss;
  ss << read_strand[0].r1_forward;

  for (long i = 1; i < static_cast<long>(read_strand.size()); ++i)
    ss << ',' << read_strand[i].r1_forward;

  return ss.str();
}


std::string
VarStats::get_r2_forward_strand_bias() const
{
  if (read_strand.size() == 0)
    return ".";

  std::stringstream ss;
  ss << read_strand[0].r2_forward;

  for (long i = 1; i < static_cast<long>(read_strand.size()); ++i)
    ss << ',' << read_strand[i].r2_forward;

  return ss.str();
}


std::string
VarStats::get_r1_reverse_strand_bias() const
{
  if (read_strand.size() == 0)
    return ".";

  std::stringstream ss;
  ss << read_strand[0].r1_reverse;

  for (long i = 1; i < static_cast<long>(read_strand.size()); ++i)
    ss << ',' << read_strand[i].r1_reverse;

  return ss.str();
}


std::string
VarStats::get_r2_reverse_strand_bias() const
{
  if (read_strand.size() == 0)
    return ".";

  std::stringstream ss;
  ss << read_strand[0].r2_reverse;

  for (long i = 1; i < static_cast<long>(read_strand.size()); ++i)
    ss << ',' << read_strand[i].r2_reverse;

  return ss.str();
}


void
VarStats::add_mapq(uint8_t const new_mapq)
{
  // Ignore MapQ == 255, that means MapQ is unavailable
  if (new_mapq < 255u)
    mapq_squared += static_cast<uint64_t>(new_mapq) * static_cast<uint64_t>(new_mapq);
}


/** Non-member functions */
template <class T>
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
  boost::split(splitted, bias, [](char c){
      return c == ',';
    });
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
get_accumulated_strand_bias(std::map<std::string, std::string> const & infos,
                            std::string const & bias)
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
