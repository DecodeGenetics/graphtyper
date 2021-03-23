#include <sstream> // std::ostringstream
#include <string> // std::string, std::stol

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/absolute_position.hpp>
#include <graphtyper/graph/genomic_region.hpp>
#include <graphtyper/graph/graph.hpp>
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
has_matching_longest_prefix(std::vector<char> const & ref, std::vector<gyper::Alt> const & alts)
{
  // Check reference
  for (uint32_t a = 0; a < alts.size(); ++a)
  {
    if (prefix_match(ref, alts[a].seq))
      return true;
  }

  assert(alts.size() > 0);

  for (uint32_t i = 0; i < alts.size() - 1; ++i)
  {
    for (uint32_t j = i + 1; j < alts.size(); ++j)
    {
      if (prefix_match(alts[i].seq, alts[j].seq))
      {
        // Detect if there are duplicate alt alleles
        if (alts[i].seq == alts[j].seq)
        {
          BOOST_LOG_TRIVIAL(error) << __HERE__ << " Duplicated alt. alleles detected. "
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

GenomicRegion::GenomicRegion()
  : chr("N/A"), begin(0), end(AS_LONG_AS_POSSIBLE)
{}


GenomicRegion::GenomicRegion(std::string const & region)
  : chr("N/A"), begin(0), end(AS_LONG_AS_POSSIBLE)
{
  if (region.size() == 0 || (region.size() == 1 && region[0] == '.'))
    return;

  assert(std::count(region.begin(), region.end(), ':') <= 1);
  //assert(std::count(region.begin(), region.end(), '-') <= 1);

  if (std::count(region.begin(), region.end(), ':') == 0)
  {
    chr = region;
  }
  else
  {
    std::size_t const colon = region.find(':');
    chr = region.substr(0, colon);

    if (std::count(region.begin(), region.end(), '-') == 0)
    {
      begin = std::stol(region.substr(colon + 1, region.size() - colon));
    }
    else
    {
      std::size_t const dash = region.find('-', colon + 1);

      begin = std::stol(region.substr(colon + 1, dash - colon));
      end = std::stol(region.substr(dash + 1, region.size() - dash));
    }
  }


  if (begin != 0)
  {
    // Switch to 0-based indexing
    --begin;
  }
}


GenomicRegion::GenomicRegion(std::string const & chrom, long _begin, long _end)
  : chr(chrom)
  , begin(_begin)
  , end(_end)
{
  if (begin != 0)
  {
    // Switch to 0-based indexing
    --begin;
  }
}


void
GenomicRegion::clear()
{
  chr.clear();
  begin = 0;
  end = 0;
}


void
GenomicRegion::pad(long bases)
{
  begin = std::max(static_cast<long>(begin) - bases, 0l);
  end += bases;
}


void
GenomicRegion::pad_end(long bases)
{
  end += bases;
}


std::string
GenomicRegion::to_string() const
{
  std::ostringstream ss;

  if (end == AS_LONG_AS_POSSIBLE)
  {
    ss << chr << ":" << (begin + 1);
  }
  else
  {
    ss << chr << ":" << (begin + 1) << "-" << end;
  }

  return ss.str();
}


long
GenomicRegion::get_absolute_begin_position() const
{
  return absolute_pos.get_absolute_position(chr, begin + 1);
}


long
GenomicRegion::get_absolute_end_position() const
{
  return absolute_pos.get_absolute_position(chr, end + 1);
}


long
GenomicRegion::get_absolute_position(std::string const & chromosome, long contig_position) const
{
  return absolute_pos.get_absolute_position(chromosome, contig_position);
}


long
GenomicRegion::get_absolute_position(long contig_position) const
{
  return absolute_pos.get_absolute_position(chr, contig_position);
}


std::pair<std::string, long>
GenomicRegion::get_contig_position(long absolute_position, Graph const & graph) const
{
  return absolute_pos.get_contig_position(absolute_position, graph.contigs);
}


void
GenomicRegion::check_if_var_records_match_reference_genome(std::vector<VarRecord> const & var_records,
                                                           std::vector<char> const & reference)
{
  for (auto const & record : var_records)
  {
    long const pos = record.pos;

    if (pos >= static_cast<long>(GenomicRegion::begin + reference.size()) || GenomicRegion::begin > pos)
      continue;

    assert(record.ref.seq.size() > 0);
    assert(record.ref.seq[0] != '.');

    auto start_it = reference.begin() + (pos - GenomicRegion::begin);
    std::size_t const ref_len = std::min(static_cast<std::size_t>(std::distance(start_it, reference.end())),
                                         record.ref.seq.size());

    std::vector<char> const reference_ref(start_it, start_it + ref_len);

    if (reference_ref.size() > 0 && !prefix_match(record.ref.seq, reference_ref))
    {
      BOOST_LOG_TRIVIAL(warning) << __HERE__ << " Record @ position " << pos
                                 << " (REF = "
                                 << std::string(record.ref.seq.begin(), record.ref.seq.end()) << ')'
                                 << " did not match the reference genome ("
                                 << std::string(reference_ref.begin(), reference_ref.end()) << ")";
    }
  }
}


void
GenomicRegion::add_reference_to_record_if_they_have_a_matching_prefix(VarRecord & var_record,
                                                                      std::vector<char> const & reference)
{
  if (var_record.is_sv)
    return;

  auto start_it = reference.begin() + (var_record.pos - GenomicRegion::begin) + var_record.ref.seq.size();
  assert(std::distance(reference.begin(), start_it) <= static_cast<int64_t>(reference.size()));

  while (start_it != reference.end() &&
         *start_it != 'N' &&
         has_matching_longest_prefix(var_record.ref.seq, var_record.alts))
  {
    // Add base to back
    var_record.ref.seq.push_back(*start_it);

    for (auto & alt : var_record.alts)
      alt.seq.push_back(*start_it);

    ++start_it;
  }
}


template <typename Archive>
void
GenomicRegion::serialize(Archive & ar, unsigned const /*version*/)
{
  ar & chr;
  ar & begin;
  ar & end;
}


template void GenomicRegion::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive &,
                                                                        const unsigned int);
template void GenomicRegion::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive &,
                                                                        const unsigned int);

} // namespace gyper

BOOST_CLASS_VERSION(gyper::GenomicRegion, 1)
