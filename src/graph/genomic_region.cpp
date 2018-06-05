#include <ostream>
#include <sstream>
#include <utility>

#include <seqan/basic.h>
#include <seqan/bam_io.h>

#include <graphtyper/graph/absolute_position.hpp>
#include <graphtyper/graph/genomic_region.hpp>
#include <graphtyper/graph/var_record.hpp>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/log/trivial.hpp>


namespace
{


bool
prefix_match(std::vector<char> const & seq1, std::vector<char> const & seq2)
{
  auto it1 = seq1.begin();
  auto it2 = seq2.begin();

  while (it1 != seq1.end() && it2 != seq2.end())
  {
    if (*it1 != *it2)
      return false;

    ++it1;
    ++it2;
  }

  return true;
}


bool
has_matching_longest_prefix(std::vector<char> const & ref, std::vector<std::vector<char> > const & alts)
{
  // Check reference
  for (uint32_t a = 0; a < alts.size(); ++a)
  {
    if (prefix_match(ref, alts[a]))
      return true;
  }

  assert (alts.size() > 0);

  for (uint32_t i = 0; i < alts.size() - 1; ++i)
  {
    for (uint32_t j = i + 1; j < alts.size(); ++j)
    {
      if (prefix_match(alts[i], alts[j]))
      {
        // Detect if there are duplicate alt alleles
        if (alts[i] == alts[j])
        {
          BOOST_LOG_TRIVIAL(fatal) << "Duplicated alt. alleles detected. "
                                   << "Aborting constructing graph.";
          std::exit(1);
        }

        return true;
      }
    }
  }

  return false;
}


} // anon namespace


namespace gyper
{

GenomicRegion::GenomicRegion(uint32_t _region_to_refnode)
  : rID(0), chr("chr1"), begin(0), end(CHR01_LENGTH - 1), region_to_refnode(_region_to_refnode)
{ }


GenomicRegion::GenomicRegion(uint16_t && r, std::string && c, uint32_t && b, uint32_t && e, uint32_t _region_to_refnode)
  : rID(std::forward<uint16_t>(r))
  , chr(std::forward<std::string>(c))
  , begin(std::forward<uint32_t>(b))
  , end(std::forward<uint32_t>(e))
  , region_to_refnode(_region_to_refnode)
{}


GenomicRegion::GenomicRegion(std::string region, uint32_t _region_to_refnode)
  : rID(0), chr("chr1"), begin(0), end(AS_LONG_AS_POSSIBLE), region_to_refnode(_region_to_refnode)
{
  if (region == std::string("."))
    return;

  assert(std::count(region.begin(), region.end(), ':') <= 1);
  assert(std::count(region.begin(), region.end(), '-') <= 1);

  if (std::count(region.begin(), region.end(), ':') == 0)
  {
    GenomicRegion::chr = region;
  }
  else
  {
    std::size_t colon = region.find(':');

    GenomicRegion::chr = region.substr(0, colon);

    if (std::count(region.begin(), region.end(), '-') == 0)
    {
      GenomicRegion::begin = stoi(region.substr(colon + 1, region.size() - colon));
    }
    else
    {
      std::size_t dash = region.find('-');

      GenomicRegion::begin = stoi(region.substr(colon + 1, dash - colon));
      GenomicRegion::end = stoi(region.substr(dash + 1, region.size() - dash));
    }
  }


  if (GenomicRegion::begin != 0)
  {
    // Switch to 0-based indexing
    --GenomicRegion::begin;
  }
}


std::string
GenomicRegion::to_string() const
{
  std::ostringstream ss;

  if (GenomicRegion::end == AS_LONG_AS_POSSIBLE)
  {
    if (GenomicRegion::begin == 0)
    {
      ss << chr;
    }
    else
    {
      ss << chr << ":" << (GenomicRegion::begin + 1);
    }
  }
  else
  {
    ss << chr << ":" << (GenomicRegion::begin + 1) << "-" << GenomicRegion::end;
  }

  return ss.str();
}


uint32_t
GenomicRegion::get_absolute_begin_position() const
{
  return absolute_pos.get_absolute_position(chr, begin + 1);
}


uint32_t
GenomicRegion::get_absolute_end_position() const
{
  return absolute_pos.get_absolute_position(chr, end + 1);
}


uint32_t
GenomicRegion::get_absolute_position(std::string const & chromosome, uint32_t contig_position) const
{
  return absolute_pos.get_absolute_position(chromosome, contig_position);
}


uint32_t
GenomicRegion::get_absolute_position(uint32_t contig_position) const
{
  return absolute_pos.get_absolute_position(chr, contig_position);
}


std::pair<std::string, uint32_t>
GenomicRegion::get_contig_position(uint32_t absolute_position) const
{
  return absolute_pos.get_contig_position(absolute_position);
}


void
GenomicRegion::check_if_var_records_match_reference_genome(std::vector<VarRecord> const & var_records, std::vector<char> const & reference)
{
  for (auto const & record : var_records)
  {
    uint32_t const pos = record.pos;

    if (pos >= GenomicRegion::begin + reference.size())
      continue;

    assert(record.ref.size() > 0);
    assert(record.ref[0] != '.');

    auto start_it = reference.begin() + (pos - GenomicRegion::begin);
    std::size_t const ref_len = std::min(static_cast<std::size_t>(std::distance(start_it, reference.end())),
                                         record.ref.size()
                                         );

    std::vector<char> const reference_ref(start_it, start_it + ref_len);

    if (reference_ref.size() > 0 && !prefix_match(record.ref, reference_ref))
    {
      BOOST_LOG_TRIVIAL(warning) << "[graphtyper::genomic_region] Record @ position " << pos
                                 << " (REF = " << std::string(record.ref.begin(), record.ref.end()) << ')'
                                 << " did not match the reference genome ("
                                 << std::string(reference_ref.begin(), reference_ref.end()) << ")";
      //std::exit(1);
    }
  }
}


void
GenomicRegion::add_reference_to_record_if_they_have_a_matching_prefix(VarRecord & var_record, std::vector<char> const & reference)
{
  if (var_record.is_sv)
    return;

  auto start_it = reference.begin() + (var_record.pos - GenomicRegion::begin) + var_record.ref.size();
  assert(std::distance(reference.begin(), start_it) <= static_cast<int64_t>(reference.size()));

  while (start_it != reference.end() && *start_it != 'N' &&
    has_matching_longest_prefix(var_record.ref, var_record.alts)
    )
  {
    // Add base to back
    var_record.ref.push_back(*start_it);

    for (auto & alt : var_record.alts)
      alt.push_back(*start_it);

    ++start_it;
  }
}


template <typename Archive>
void
GenomicRegion::serialize(Archive & ar, unsigned const /*version*/)
{
  ar & rID;
  ar & chr;
  ar & begin;
  ar & end;
  ar & region_to_refnode;
}


template void GenomicRegion::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive &, const unsigned int);
template void GenomicRegion::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive &, const unsigned int);

} // namespace gyper

BOOST_CLASS_VERSION(gyper::GenomicRegion, 1)
