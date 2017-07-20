#include <iostream>
#include <sstream>
#include <fstream>

#include <boost/log/trivial.hpp>

#include <graphtyper/constants.hpp>
#include <graphtyper/utilities/options.hpp>
#include <graphtyper/utilities/sam_reader.hpp>

/*
namespace
{

// Deletes all tags, except Read Group, and ADapters (RG, AD)
void
delete_tags(seqan::BamAlignmentRecord & record)
{
  seqan::BamTagsDict tags_dict(record.tags);
  seqan::buildIndex(tags_dict);
  seqan::StringSet<seqan::CharString> keys_to_remove;

  for (std::size_t i = 0; i < seqan::length(tags_dict); ++i)
  {
    seqan::CharString key = tags_dict[i];
    seqan::resize(key, 2);

    if (key != "RG" && key != "AD")
      appendValue(keys_to_remove, key);
  }

  for (auto it = seqan::begin(keys_to_remove); it != seqan::end(keys_to_remove); ++it)
  {
    seqan::eraseTag(tags_dict, *it);
  }
}


} // namespace anon
*/


namespace gyper
{

SamReader::SamReader(std::string const & hts_path, std::vector<std::string> const & _regions)
  : hts_file(hts_path.c_str())
  , regions(_regions)
{
  if (!(regions.size() == 1 && regions[0] == std::string(".")) && !seqan::loadIndex(hts_file))
  {
    // Try to build it if we cannot open it
    seqan::buildIndex(hts_file);
    seqan::loadIndex(hts_file);
  }

  assert(regions.size() > 0);
  seqan::setRegion(hts_file, regions[r].c_str());
}


TReads
SamReader::read_N_reads(std::size_t const N)
{
  TReads reads;
  seqan::BamAlignmentRecord record;

  while (true)
  {
    if (seqan::readRegion(record, hts_file))
    {
      getContigName(record.rID, hts_file);
      // delete_tags(record); // Not tested

      // Remove Ns from front. Note that if the entire sequence is N this won't erase anything. However, the sequence will be erased in the next step.
      for (std::size_t i = 0; i < seqan::length(record.seq); ++i)
      {
        if (record.seq[i] != seqan::Iupac('N'))
        {
          if (i > 0)
          {
            seqan::erase(record.seq, 0, i);
            seqan::erase(record.qual, 0, i);
            seqan::clear(record.cigar); // Clear CIGAR as it is incorrect now
          }

          break;
        }
      }

      // The minimum read length is 2 overlapping k-mers
      std::size_t const MIN_READ_LENGTH = 2 * K - 1;

      // Remove Ns from back
      while (seqan::length(record.seq) >= MIN_READ_LENGTH && seqan::back(record.seq) == seqan::Iupac('N'))
      {
        seqan::eraseBack(record.seq);
        seqan::eraseBack(record.qual);
        seqan::clear(record.cigar);
      }

      // Require the read to be at least MIN_READ_LENGTH, otherwise skip it
      if (seqan::length(record.seq) < MIN_READ_LENGTH)
        continue;

      assert(seqan::length(record.seq) == seqan::length(record.qual));
      insert_reads(reads, std::move(record));

      if (reads.size() >= N)
      {
        return reads;
      }
    }
    else if (!second_file && (r + 1) < regions.size())
    {
      ++r;
      seqan::setRegion(hts_file, regions[r].c_str());
    }
    else
    {
      break;
    }
  }

  // If reads.size() == 0 this is the last iteration
  if (reads.size() == 0)
  {
    reads = std::move(unpaired_reads);

    for (auto it = reads_first.begin(); it != reads_first.end(); ++it)
    {
      seqan::BamAlignmentRecord empty_record;
      reads.push_back({it->second, empty_record});
    }

    reads_first.clear();
    unpaired_reads.clear();
  }

  return reads;
}


void
SamReader::insert_reads(TReads & reads, seqan::BamAlignmentRecord && record)
{
  if ((hts_file.hts_record->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FSUPPLEMENTARY | BAM_FDUP)) != 0u)
    return;

  if (!seqan::hasFlagMultiple(record))
  {
    // Read is not paired
    seqan::BamAlignmentRecord empty_record;

    // Make it the first read, if it isn't it already
    //if (!seqan::hasFlagFirst(record))
    //  seqan::toggleFlagFirst(record);

    // This read is read 1
    if (seqan::hasFlagRC(record))
    {
      // Read 1 is reversed
      seqan::reverseComplement(record.seq);
      seqan::reverse(record.qual);
    }

    unpaired_reads.push_back({record, empty_record});
    return;
  }

  if (seqan::hasFlagFirst(record))
  {
    // This read is read 1
    if (seqan::hasFlagRC(record))
    {
      // Read 1 is reversed
      seqan::reverseComplement(record.seq);
      seqan::reverse(record.qual);
    }
  }
  else if (!seqan::hasFlagRC(record))
  {
    // Read 2 is not reversed
    seqan::reverseComplement(record.seq);
    seqan::reverse(record.qual);
  }

  auto results_it = reads_first.find(record.qName);

  if (results_it != reads_first.end())
  {
    if (seqan::hasFlagFirst(record))
      reads.push_back({record, results_it->second});
    else
      reads.push_back({results_it->second, record});

    reads_first.erase(results_it);
  }
  else
  {
    reads_first[record.qName] = std::move(record);
  }
}


} // namespace gyper
