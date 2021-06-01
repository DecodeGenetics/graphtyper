#include <string> // std::string
#include <vector> // std::vector

#include <seqan/hts_io.h>

#include <graphtyper/constants.hpp>
#include <graphtyper/utilities/logging.hpp>
#include <graphtyper/utilities/options.hpp>
#include <graphtyper/utilities/sam_reader.hpp>

namespace gyper
{
SamReader::SamReader(std::string const & hts_path, std::vector<std::string> const & _regions) :
  r(0), hts_file(hts_path.c_str()), regions(_regions)
{
  if (!(regions.size() == 1 && regions[0] == std::string(".")) && !seqan::loadIndex(hts_file))
  {
    // Try to build it if we cannot open it
    seqan::buildIndex(hts_file);
    seqan::loadIndex(hts_file);
  }

  if (regions.size() == 0)
    regions.push_back(".");

  seqan::setRegion(hts_file, regions[r].c_str());
}

TReads SamReader::read_N_reads(std::size_t const N)
{
  TReads reads;
  seqan::BamAlignmentRecord record;

  while (true)
  {
    if (seqan::readRegion(record, hts_file))
    {
      // The minimum read length is 2 overlapping k-mers
      // std::size_t const MIN_READ_LENGTH = 2 * K - 1;
      //
      // Require the read to be at least MIN_READ_LENGTH, otherwise skip it
      // if (seqan::length(record.seq) < MIN_READ_LENGTH)
      //  continue;
      assert(seqan::length(record.seq) >= 2 * K - 1);
      assert(seqan::length(record.seq) == seqan::length(record.qual));
      insert_reads(reads, std::move(record));

      if (reads.size() >= N)
        return reads;
    }
    else if (r + 1 < regions.size())
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

void SamReader::insert_reads(TReads & reads, seqan::BamAlignmentRecord && record)
{
  if ((hts_file.hts_record->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FSUPPLEMENTARY | BAM_FDUP)) != 0u)
    return;

  if (!seqan::hasFlagMultiple(record))
  {
    // Read is not paired
    seqan::BamAlignmentRecord empty_record;

    // Make it the first read, if it isn't it already
    // if (!seqan::hasFlagFirst(record))
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
    {
      if (seqan::hasFlagFirst(results_it->second))
      {
        print_log(log_severity::warning,
                  "[graphtyper::sam_reader] Two first-in-pair reads with identical read name found. "
                  "Are you reading the same region more than once?");
      }
      else
      {
        reads.push_back({record, results_it->second});
      }
    }
    else
    {
      if (!seqan::hasFlagFirst(results_it->second))
      {
        print_log(log_severity::warning,
                  "[graphtyper::sam_reader] Two second-in-pair reads with identical read name found. "
                  "Are you reading the same region more than once?");
      }
      else
      {
        reads.push_back({results_it->second, record});
      }
    }

    reads_first.erase(results_it);
  }
  else
  {
    reads_first[record.qName] = std::move(record);
  }
}

} // namespace gyper
