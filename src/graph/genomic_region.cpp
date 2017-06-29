#include <ostream>
#include <sstream>
#include <utility>

#include <seqan/basic.h>
#include <seqan/bam_io.h>

#include <graphtyper/graph/genomic_region.hpp>

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
        return true;
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

  assert (std::count(region.begin(), region.end(), ':') <= 1);
  assert (std::count(region.begin(), region.end(), '-') <= 1);

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
  for (auto record_it = var_records.cbegin(); record_it != var_records.cend(); ++record_it)
  {
    uint32_t const pos = record_it->pos;
    assert(pos < GenomicRegion::begin + reference.size());

    std::vector<char> vcf_ref(record_it->ref);
    assert(vcf_ref.size() > 0);
    assert(vcf_ref[0] != '.');

    auto start_it = reference.begin() + (pos - GenomicRegion::begin);
    auto end_it = start_it + vcf_ref.size();
    std::vector<char> reference_ref(start_it, end_it);

    if (reference_ref != vcf_ref)
    {
      BOOST_LOG_TRIVIAL(error) << "[graphtyper::genomic_region] Record @ position " << std::to_string(pos)
                               << " (REF = " << std::string(vcf_ref.begin(), vcf_ref.end())
                               << ") did not match the reference genome (" << std::string(reference_ref.begin(), reference_ref.end()) << ")";
      std::exit(1);
    }
  }
}


void
GenomicRegion::add_reference_to_record_if_they_have_a_matching_prefix(VarRecord & var_record, std::vector<char> const & reference)
{
  auto start_it = reference.begin() + (var_record.pos - GenomicRegion::begin) + var_record.ref.size();
  assert(std::distance(reference.begin(), start_it) <= static_cast<int64_t>(reference.size()));

  while (start_it != reference.end() && *start_it != 'N' && has_matching_longest_prefix(var_record.ref, var_record.alts))
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
GenomicRegion::serialize(Archive & ar, const unsigned int version)
{
  ar & rID;
  ar & chr;
  ar & begin;
  ar & end;

  if (version > 0)
  {
    ar & region_to_refnode;
  }
  else
  {
    region_to_refnode = 0;
  }
}


template void GenomicRegion::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive &, const unsigned int);
template void GenomicRegion::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive &, const unsigned int);

} // namespace gyper

BOOST_CLASS_VERSION(gyper::GenomicRegion, 1)
