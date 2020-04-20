#include <iostream>
#include <limits>
#include <map>
#include <unordered_map>
#include <string>

#include <htslib/sam.h>

#include <seqan/file.h>
#include <seqan/bam_io.h>
#include <seqan/align.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/store.h>

#include <paw/parser.hpp>

#include <boost/log/trivial.hpp>

#include <graphtyper/utilities/bamshrink.hpp>
#include <graphtyper/utilities/options.hpp>


namespace bamshrink
{

#ifndef NDEBUG
bool constexpr CHANGE_READ_NAMES = false;
#else
bool constexpr CHANGE_READ_NAMES = true;
#endif // NDEBUG

} // namespace bamshrink


namespace
{

void
removeHardClipped(seqan::String<seqan::CigarElement<> > & cigar)
{
  long n_cigar = seqan::length(cigar);

  if (n_cigar >= 1 && cigar[0].operation == 'H')
  {
    seqan::erase(cigar, 0);
    --n_cigar;
  }

  if (n_cigar >= 2 && cigar[n_cigar - 1].operation == 'H')
    seqan::eraseBack(cigar);
}


void
binarizeQual(seqan::CharString & qual)
{
  for (auto & q : qual)
    q = (q - 33 >= 24) ? '?' : ',';
}


bool
is_clipped_both_ends(seqan::String<seqan::CigarElement<> > const & cigar, long const min_clip = 15)
{
  return seqan::length(cigar) >= 1 &&
         cigar[0].operation == 'S' &&
         cigar[seqan::length(cigar) - 1].operation == 'S' &&
         static_cast<long>(cigar[0].count + cigar[seqan::length(cigar) - 1].count) >= min_clip;
}


bool
is_one_end_clipped(seqan::String<seqan::CigarElement<> > const & cigar, long const min_clip = 0)
{
  return seqan::length(cigar) == 0 ||
         (cigar[0].operation == 'S' && cigar[0].count >= min_clip) ||
         (cigar[seqan::length(cigar) - 1].operation == 'S' && cigar[seqan::length(cigar) - 1].count >= min_clip);
}


// returns true if the alignment is good
bool
process_tags(seqan::BamAlignmentRecord & record, seqan::CharString & new_tags, bamshrink::Options const & opts)
{
  unsigned i = 0;
  int64_t as = -1;
  int64_t xs = -1;

  while (i < seqan::length(record.tags))
  {
    //std::string curr_tag(2, record.tags[i]);
    //curr_tag[1] = record.tags[i+1];
    i += 3;
    char type = record.tags[i - 1];
    bool is_as = false;
    bool is_xs = false;

    // Check for AS and XS
    if (record.tags[i - 2] == 'S')
    {
      if (record.tags[i - 3] == 'A')
        is_as = true;
      else if (record.tags[i - 3] == 'X')
        is_xs = true;
    }

    switch (type)
    {
    case 'A':
    {
      // A printable character
      ++i;
      break;
    }

    case 'Z':
    {
      // A string! Let's loop it until qNULL
      // Check for RG
      if (record.tags[i - 3] == 'R' && record.tags[i - 2] == 'G')
      {
        auto begin_it = begin(record.tags) + i - 3;

        while (record.tags[i] != '\0' && record.tags[i] != '\n')
          ++i;

        ++i;
        auto end_it = begin(record.tags) + i;
        //unsigned const old_length = length(new_tags);
        resize(new_tags, std::distance(begin_it, end_it), seqan::Exact());
        seqan::arrayCopyForward(begin_it, end_it, begin(new_tags));
      }
      else
      {
        while (record.tags[i] != '\0' && record.tags[i] != '\n')
          ++i;

        ++i;
      }

      break;
    }

    case 'c':
    {
      if (is_as)
      {
        int8_t num;
        memcpy(&num, &record.tags[i], sizeof(int8_t));
        as = num;
      }
      else if (is_xs)
      {
        int8_t num;
        memcpy(&num, &record.tags[i], sizeof(int8_t));
        xs = num;
      }

      i += sizeof(int8_t);
      break;
    }

    case 'C':
    {
      if (is_as)
      {
        uint8_t num;
        memcpy(&num, &record.tags[i], sizeof(uint8_t));
        as = num;
      }
      else if (is_xs)
      {
        uint8_t num;
        memcpy(&num, &record.tags[i], sizeof(uint8_t));
        xs = num;
      }

      i += sizeof(uint8_t);
      break;
    }

    case 's':
    {
      if (is_as)
      {
        int16_t num;
        memcpy(&num, &record.tags[i], sizeof(int16_t));
        as = num;
      }
      else if (is_xs)
      {
        int16_t num;
        memcpy(&num, &record.tags[i], sizeof(int16_t));
        xs = num;
      }

      i += sizeof(int16_t);
      break;
    }

    case 'S':
    {
      if (is_as)
      {
        uint16_t num;
        memcpy(&num, &record.tags[i], sizeof(uint16_t));
        as = num;
      }
      else if (is_xs)
      {
        uint16_t num;
        memcpy(&num, &record.tags[i], sizeof(uint16_t));
        xs = num;
      }

      i += sizeof(uint16_t);
      break;
    }

    case 'i':
    {
      if (is_as)
      {
        int32_t num;
        memcpy(&num, &record.tags[i], sizeof(int32_t));
        as = num;
      }
      else if (is_xs)
      {
        int32_t num;
        memcpy(&num, &record.tags[i], sizeof(int32_t));
        xs = num;
      }

      i += sizeof(int32_t);
      break;
    }

    case 'I':
    {
      if (is_as)
      {
        uint32_t num;
        memcpy(&num, &record.tags[i], sizeof(uint32_t));
        as = num;
      }
      else if (is_xs)
      {
        uint32_t num;
        memcpy(&num, &record.tags[i], sizeof(uint32_t));
        xs = num;
      }

      i += sizeof(uint32_t);
      break;
    }

    case 'f':
    {
      i += sizeof(float);
      break;
    }

    default:
    {
      i = length(record.tags);       // Unkown tag, stop
      break;
    }
    }
  }

  if (as != -1 && xs != -1 && (!seqan::hasFlagMultiple(record) || seqan::hasFlagNextUnmapped(record)))
  {
    if (as <= xs)
      return false;

    long matches = 0;
    long indels = 0;

    for (auto const & c : record.cigar)
    {
      if (c.operation == 'M')
        matches += c.count;
      else if (c.operation == 'D' || c.operation == 'I')
        indels += c.count + 2; // Extra 2 for each event
    }

    if ((as + opts.as_filter_threshold) <= (matches - indels))
      return false;
  }

  return true;
}


} // anon namespace


namespace seqan
{

bool
operator<(seqan::BamAlignmentRecord const & a, seqan::BamAlignmentRecord const & b)
{
  return a.beginPos < b.beginPos;
}


}


namespace std
{

template <>
struct hash<seqan::String<char> >
{
  std::size_t
  operator()(seqan::String<char> const & s) const
  {
    std::size_t seed = 42;

    for (char c : s)
      seed ^= c + 0x9e3779b9 + (c << 6) + (c >> 2);

    return seed;
  }


};

} // namespace std

using namespace seqan;


namespace bamshrink
{

void
makeUnpaired(BamAlignmentRecord & record)
{
  //record.tLen = 0; // Removed for insert size distributions
  record.pNext = record.INVALID_POS;
  record.rNextId = record.INVALID_REFID;
  //unset: FlagNextUnmapped, FlagAllProper,FlagMultiple, FlagNextRC:
  record.flag &= ~BAM_FLAG_NEXT_UNMAPPED;
  record.flag &= ~BAM_FLAG_ALL_PROPER;
  record.flag &= ~BAM_FLAG_MULTIPLE;
  record.flag &= ~BAM_FLAG_NEXT_RC;
}


long
countMatchingBases(String<CigarElement<> > const & cigarString)
{
  long numOfMatches = 0;

  for (auto const & cigar : cigarString)
  {
    if (cigar.operation == 'M')
      numOfMatches += cigar.count;
  }

  return numOfMatches;
}


#ifndef NDEBUG
bool
cigarAndSeqMatch(BamAlignmentRecord & record)
{
  unsigned counter = 0;

  for (auto const & cigar : record.cigar)
  {
    if (cigar.operation != 'D')
      counter += cigar.count;
  }

  return length(record.seq) == counter;
}


#endif // ifndef NDEBUG


void
resetCigarStringEnd(String<CigarElement<> > & cigarString, unsigned nRemoved)
{
  if (length(cigarString) == 0)
    return;

  if (cigarString[length(cigarString) - 1].operation == 'D')
  {
    eraseBack(cigarString);

    if (length(cigarString) == 0)
      return;
  }

  auto & cigar_end = cigarString[length(cigarString) - 1];

  if (cigar_end.count > nRemoved)
  {
    cigar_end.count -= nRemoved;
  }
  else if (cigar_end.count == nRemoved)
  {
    eraseBack(cigarString);

    if (length(cigarString) > 0 && cigarString[length(cigarString) - 1].operation == 'D')
      eraseBack(cigarString);
  }
  else
  {
    unsigned nLeft = nRemoved - cigar_end.count;
    eraseBack(cigarString);
    resetCigarStringEnd(cigarString, nLeft);
  }
}


// Returns the amount of reference bases (cigars D or M) removed from cigar
unsigned
resetCigarStringBegin(String<CigarElement<> > & cigarString, unsigned nRemoved)
{
  if (length(cigarString) == 0)
    return 0;

  unsigned removed;

  if (cigarString[0].operation == 'D')
  {
    removed = cigarString[0].count;
    erase(cigarString, 0);

    if (length(cigarString) == 0)
      return removed;
  }
  else
  {
    removed = 0;
  }

  if (cigarString[0].count > nRemoved)
  {
    cigarString[0].count -= nRemoved;

    if (cigarString[0].operation == 'M')
      removed += nRemoved;
  }
  else if (cigarString[0].count == nRemoved)
  {
    if (cigarString[0].operation == 'M')
      removed += cigarString[0].count;

    erase(cigarString, 0);

    if (length(cigarString) == 0)
      return removed;

    if (cigarString[0].operation == 'D')
    {
      removed += cigarString[0].count;
      erase(cigarString, 0);
    }
  }
  else
  {
    if (cigarString[0].operation == 'M')
      removed += cigarString[0].count;

    unsigned nLeft = nRemoved - cigarString[0].count;
    erase(cigarString, 0);

    if (length(cigarString) == 0)
      return removed;
    else
      return removed + resetCigarStringBegin(cigarString, nLeft);
  }

  return removed;
}


bool
removeSoftClipped(BamAlignmentRecord & record, Options const & opts)
{
  long n_cigar = seqan::length(record.cigar);

  if (n_cigar >= 1)
  {
    auto const & first_cigar = record.cigar[0];

    if (first_cigar.operation == 'S')
    {
      erase(record.seq, 0, first_cigar.count);
      erase(record.qual, 0, first_cigar.count);
      erase(record.cigar, 0);
      --n_cigar;
    }

    if (n_cigar >= 2)
    {
      auto const & last_cigar = record.cigar[n_cigar - 1];

      if (last_cigar.operation == 'S')
      {
        long const sequence_length = length(record.seq);
        resize(record.seq, sequence_length - last_cigar.count);
        resize(record.qual, sequence_length - last_cigar.count);
        eraseBack(record.cigar);
      }
    }
  }

  if (static_cast<long>(length(record.seq)) < opts.minReadLen ||
      (record.mapQ < 25 && static_cast<long>(length(record.seq)) < opts.minReadLenMapQ0))
  {
    return false;
  }

  return true;
}


bool
removeNsAtEnds(BamAlignmentRecord & record, Options const & opts)
{
  int nOfNs = 0;

  if (record.seq[0] == 'N')
  {
    ++nOfNs;
    int idx = 1;

    while (record.seq[idx] == 'N' && idx < static_cast<long>(length(record.seq) - 1))
    {
      ++nOfNs;
      ++idx;
    }

    // Remove ns from beginning of sequence and qual fields
    erase(record.seq, 0, nOfNs);
    erase(record.qual, 0, nOfNs);

    if (!hasFlagUnmapped(record)) // Only have to fix CIGAR, beginPos and fragLen if the read is mapped
    {
      unsigned shift = resetCigarStringBegin(record.cigar, nOfNs);
      record.beginPos += shift;
    }
  }

  if (static_cast<long>(length(record.seq)) < opts.minReadLen ||
      (record.mapQ < 25 && static_cast<long>(length(record.seq)) < opts.minReadLenMapQ0))
  {
    return false;
  }

  nOfNs = 0;

  if (record.seq[length(record.seq) - 1] == 'N')
  {
    ++nOfNs;
    int idx = length(record.seq) - 2;

    while (record.seq[idx] == 'N' && idx > 0)
    {
      ++nOfNs;
      --idx;
    }

    // Remove ns from end of sequence and qual fields:
    erase(record.seq, length(record.seq) - nOfNs, length(record.seq));
    erase(record.qual, length(record.qual) - nOfNs, length(record.qual));

    if (!hasFlagUnmapped(record)) // Only have to fix CIGAR, beginPos and fragLen if the read is mapped
      resetCigarStringEnd(record.cigar, nOfNs);
  }

  if (static_cast<long>(length(record.seq)) < opts.minReadLen ||
      (record.mapQ < 25 && static_cast<long>(length(record.seq)) < opts.minReadLenMapQ0))
  {
    return false;
  }

  return true;
}


Pair<int>
findNum2Clip(BamAlignmentRecord & recordReverse, int forwardStartPos)
{
  int num2clip = 0;
  int num2shift = 0;
  unsigned cigarIndex = 0;
  long reverseStartPos = recordReverse.beginPos;
  unsigned n = 0;

  if (recordReverse.cigar[cigarIndex].operation == 'S')
  {
    num2clip = recordReverse.cigar[cigarIndex].count;
    ++cigarIndex;
  }

  while (cigarIndex < length(recordReverse.cigar))
  {
    char cigarOperation = recordReverse.cigar[cigarIndex].operation;
    n = 0;

    while (reverseStartPos < forwardStartPos && n < recordReverse.cigar[cigarIndex].count)
    {
      if (cigarOperation != 'D')
        ++num2clip;

      if (cigarOperation != 'I')
        ++reverseStartPos;

      ++n;
    }

    if (reverseStartPos == forwardStartPos)
      break;

    ++cigarIndex;
  }

  if (recordReverse.cigar[cigarIndex].operation == 'D')
    num2shift = recordReverse.cigar[cigarIndex].count - n;

  return Pair<int>(num2clip, num2shift);
}


bool
removeAdapters(BamAlignmentRecord & recordForward,
               BamAlignmentRecord & recordReverse,
               Options const & opts)
{
  //Check for soft clipped bases at beginning of forward record.
  if (removeSoftClipped(recordForward, opts) && removeSoftClipped(recordReverse, opts))
    return false;

  //delStats.nAdapterReads += 2;
  int startPosDiff = recordForward.beginPos - recordReverse.beginPos;

  if (startPosDiff < 0)
    return true;

  Pair<int> clipAndShift = findNum2Clip(recordReverse, recordForward.beginPos);
  int index = clipAndShift.i1;
  int shift = clipAndShift.i2;

  //erase from reverse read bases 0 to index
  erase(recordReverse.seq, 0, index);
  erase(recordReverse.qual, 0, index);
  resetCigarStringBegin(recordReverse.cigar, index);
  //erase from forward read bases from length(reverse.seq) to end

  if (length(recordForward.seq) > length(recordReverse.seq) && index > 0)
  {
    int forwardClip = length(recordForward.seq) - length(recordReverse.seq);
    erase(recordForward.seq, length(recordReverse.seq), length(recordForward.seq));
    erase(recordForward.qual, length(recordReverse.qual), length(recordForward.qual));
    resetCigarStringEnd(recordForward.cigar, forwardClip);
  }

  recordReverse.beginPos = recordForward.beginPos;

  if (shift > 0)
    recordReverse.beginPos += shift;

  recordForward.pNext = recordReverse.beginPos;

#ifndef NDEBUG
  if (!cigarAndSeqMatch(recordForward))
  {
    BOOST_LOG_TRIVIAL(warning) << __HERE__ << " The cigar string and sequence length don't match "
                               << "for the forward read.";
  }

  if (!cigarAndSeqMatch(recordReverse))
  {
    BOOST_LOG_TRIVIAL(warning) << __HERE__ << " The cigar string and sequence length don't match "
                               << "for the reverse read!";
  }
#endif // ifndef NDEBUG

  if (static_cast<long>(length(recordForward.seq)) < opts.minReadLen ||
      (recordForward.mapQ < 25 && static_cast<long>(length(recordForward.seq)) < opts.minReadLenMapQ0))
  {
    return false;
  }
  else
  {
    return true;
  }
}


void
qualityFilterSlice2(Options const & opts,
                    Triple<CharString, int, int> chr_start_end, // cannot be const& due to some seqan issue
                    BamFileIn & bamFileIn,
                    BamFileOut & bamFileOut,
                    long & read_num,
                    bool const is_single_contig)
{
  if (!loadIndex(bamFileIn, opts.bamIndex.c_str()))
  {
    BOOST_LOG_TRIVIAL(error) << __HERE__ << " Could not read index file " << opts.bamIndex;
    std::exit(1);
  }

  if (!setRegion(bamFileIn,
                 toCString(chr_start_end.i1),
                 chr_start_end.i2 - opts.maxFragLen,
                 chr_start_end.i3 + opts.maxFragLen
                 )
      )
  {
    BOOST_LOG_TRIVIAL(error) << __HERE__ << " Could not set region to "
                             << chr_start_end.i1 << ":"
                             << (chr_start_end.i2 + 1) << "-"
                             << (chr_start_end.i3 + 1)
                             << " when using index " << opts.bamIndex;
    std::exit(1);
  }

  bool is_paired_reads_with_no_mate{false};
  std::multiset<BamAlignmentRecord> read_set;
  std::unordered_map<seqan::String<char>, seqan::BamAlignmentRecord> read_first;
  std::vector<seqan::BamAlignmentRecord> unpaired_reads;
  std::vector<std::pair<seqan::BamAlignmentRecord, seqan::BamAlignmentRecord> > paired_reads;
  seqan::BamAlignmentRecord record;
  long first_pos = -1;
  std::vector<uint32_t> bin_counts;
  long const max_bin_sum = opts.no_filter_on_coverage ?
                           (std::numeric_limits<int>::max() / 10) :
                           static_cast<long>(opts.avgCovByReadLen * 50.0 * 2.5);

  long const max_fragment_length = opts.maxFragLen;

  auto filter_unpaired =
    [&chr_start_end, &opts](BamAlignmentRecord const & rec) -> bool
    {
      // Unpaired read that does not overlap the target region
      if (static_cast<long>(rec.beginPos + seqan::length(rec.seq)) < static_cast<long>(chr_start_end.i2) ||
          rec.beginPos > chr_start_end.i3)
      {
        return false;
      }

      if (rec.mapQ < 25 ||
          static_cast<long>(length(rec.seq)) < opts.minUnpairedReadLen ||
          is_one_end_clipped(rec.cigar, 12) ||
          is_clipped_both_ends(rec.cigar, 5) ||
          countMatchingBases(rec.cigar) < opts.minNumMatching + 5)
      {
        return false;
      }

      return true;
    };

  auto filter_paired =
    [&chr_start_end, &opts](BamAlignmentRecord const & rec) -> bool
    {
      if (opts.is_filtering_mapq0 && rec.mapQ == 0)
      {
        return false;
      }

      // Paired read that does not overlap the target region
      if (static_cast<long>(rec.beginPos + seqan::length(rec.seq)) < static_cast<long>(chr_start_end.i2) &&
          (static_cast<long>(rec.beginPos + rec.tLen) < static_cast<long>(chr_start_end.i2)))
      {
        return false;
      }

      if (static_cast<long>(rec.beginPos) > static_cast<long>(chr_start_end.i3) &&
          (static_cast<long>(rec.beginPos + rec.tLen - seqan::length(rec.seq)) > static_cast<long>(chr_start_end.i3)))
      {
        return false;
      }

      // Allow unmapped reads with mapped mates
      if (hasFlagUnmapped(rec))
      {
        //std::cerr << rec.qName << "\n";
        return true;
      }

      // Filter for paired reads
      if (static_cast<long>(length(rec.seq)) < opts.minReadLen ||
          (rec.mapQ < 30 && is_clipped_both_ends(rec.cigar, 12)) ||
          (rec.mapQ < 5 && is_one_end_clipped(rec.cigar, length(rec.seq) / 4)) ||
          is_clipped_both_ends(rec.cigar, length(rec.seq) / 2) ||
          countMatchingBases(rec.cigar) < opts.minNumMatching)
      {
        return false;
      }

      return true;
    };

  auto post_process_unpaired =
    [&bin_counts, &first_pos, &max_bin_sum, &read_num, &read_set, &opts](BamAlignmentRecord && rec) -> void
    {
      // Filter for unpaired reads
      seqan::CharString new_tags;
      bool is_good_alignment = process_tags(rec, new_tags, opts);

      if (!is_good_alignment)
        return;

      if (!removeNsAtEnds(rec, opts))
        return;

      rec.tags = std::move(new_tags);
      long const bin = (rec.beginPos - first_pos) / 50l;

      if (bin >= static_cast<long>(bin_counts.size()))
      {
        bin_counts.resize(bin + 1, 0u);
      }
      else if (bin_counts[bin] >= (max_bin_sum / 3l))
      {
        ++bin_counts[bin];
        return;
      }

      binarizeQual(rec.qual);
      removeHardClipped(rec.cigar);

      if (CHANGE_READ_NAMES)
      {
        std::ostringstream ss;
        ss << std::hex << read_num;
        rec.qName = ss.str();
        ++read_num;
      }

      ++bin_counts[bin];
      read_set.insert(std::move(rec));
    };

  auto post_process_paired =
    [&opts](BamAlignmentRecord & rec, long const read_num) -> bool
    {
      seqan::CharString new_tags;
      bool is_good_alignment = process_tags(rec, new_tags, opts);

      if (!is_good_alignment)
        return false;

      if (!removeNsAtEnds(rec, opts))
        return false;

      rec.tags = std::move(new_tags);
      binarizeQual(rec.qual);
      removeHardClipped(rec.cigar);

      if (CHANGE_READ_NAMES)
      {
        std::ostringstream ss;
        ss << std::hex << read_num;
        rec.qName = ss.str();
      }

      return true;
    };

  while (readRegion(record, bamFileIn))
  {
    if (((record.flag & gyper::Options::const_instance()->sam_flag_filter) != 0) ||
        (record.tLen != 0 && abs(record.tLen) < opts.minReadLen))
    {
      continue;
    }

    if (first_pos < 0)
    {
      if (record.beginPos < 0)
        continue;

      first_pos = record.beginPos;
    }

    if (read_first.size() > 0 &&
        static_cast<long>(record.beginPos) >
        static_cast<long>(2 * max_fragment_length + read_first.begin()->second.beginPos))
    {
      auto it = read_first.begin();

      while (it != read_first.end() &&
             (record.beginPos > max_fragment_length + it->second.beginPos + 151) &&
             (record.qName != it->first))
      {
        makeUnpaired(it->second);
        is_paired_reads_with_no_mate = true;

        if (filter_unpaired(it->second))
          post_process_unpaired(std::move(it->second));

        ++it;
      }

      read_first.erase(read_first.begin(), it);
    }

    if (read_set.size() > 0 &&
        static_cast<long>(record.beginPos) >
        static_cast<long>(3 * max_fragment_length + read_set.begin()->beginPos))
    {
      auto it = read_set.begin();

      while (it != read_set.end() && (record.beginPos > max_fragment_length + it->beginPos + 151))
      {
        long const bin1 = (it->beginPos - first_pos) / 50l;
        long const bin2 = (it->pNext - first_pos) / 50l;

        if (bin_counts[bin1] < (opts.SUPER_HI_DEPTH * max_bin_sum) ||
            (hasFlagMultiple(*it) && bin_counts[bin2] < (opts.SUPER_HI_DEPTH * max_bin_sum)))
        {
          writeRecord(bamFileOut, *it);
        }

        ++it;
      }

      read_set.erase(read_set.begin(), it);
    }

    // If there is only one interval the header was changed
    if (is_single_contig)
    {
      if (record.rNextId == record.rID)
      {
        record.rID = 0;
        record.rNextId = 0;
      }
      else
      {
        record.rID = 0;
        record.rNextId = 1; // Such that they are not the same
      }
    }

    if ((hasFlagUnmapped(record) || hasFlagNextUnmapped(record)) && hasFlagRC(record) == hasFlagNextRC(record))
    {
      seqan::reverseComplement(record.seq);
      seqan::reverse(record.qual);
      record.flag ^= 16;
    }

    // determine which reads are unpaired
    if (record.rID != record.rNextId ||
        hasFlagRC(record) == hasFlagNextRC(record) ||
        abs(record.tLen) > max_fragment_length ||
        (record.tLen > 0 && hasFlagRC(record)) ||
        (record.tLen < 0 && !hasFlagRC(record)))
    {
      makeUnpaired(record);
    }

    if (!hasFlagMultiple(record))
    {
      // Unpaired read
      if (filter_unpaired(record))
        post_process_unpaired(std::move(record));

      continue;
    }

    // Paired reads
    if (!filter_paired(record))
      continue;

    auto find_it = read_first.find(record.qName);

    if (find_it == read_first.end())
    {
      // No point in waiting if this read comes record
      assert(record.pNext != 0);

      if (record.pNext >= record.beginPos)
        read_first[record.qName] = std::move(record);

      continue;
    }

    long const bin1 = (record.beginPos - first_pos) / 50l;
    long const bin2 = (find_it->second.beginPos - first_pos) / 50l;

    {
      long const max_bin = std::max(bin1, bin2);

      if (max_bin >= static_cast<long>(bin_counts.size()))
        bin_counts.resize(max_bin + 1, 0u);
    }

    ++bin_counts[bin1];
    ++bin_counts[bin2];

    if (bin_counts[bin1] < max_bin_sum)
    {
      if (bin_counts[bin2] < max_bin_sum)
      {
        bool is_ok;

        if (record.tLen == 0 ||
            abs(record.tLen) > static_cast<long>(std::max(length(record.seq), length(find_it->second.seq))))
          is_ok = true;
        else if (hasFlagRC(record))
          is_ok = removeAdapters(find_it->second, record, opts);
        else
          is_ok = removeAdapters(record, find_it->second, opts);

        if (is_ok && post_process_paired(record, read_num) && post_process_paired(find_it->second, read_num))
        {
          if ((!hasFlagUnmapped(record) && !hasFlagUnmapped(find_it->second)) ||
              (hasFlagUnmapped(record) && filter_unpaired(find_it->second)) ||
              (hasFlagUnmapped(find_it->second) && filter_unpaired(record)))
          {
            ++read_num; // Only increase read_num if we actually add the reads
            read_set.insert(std::move(record));
            read_set.insert(std::move(find_it->second));
          }
        }
      }
      else if (bin_counts[bin1] < (max_bin_sum / 3))
      {
        makeUnpaired(record);

        if (filter_unpaired(record))
          post_process_unpaired(std::move(record));
      }
    }
    else if (bin_counts[bin2] < (max_bin_sum / 3))
    {
      makeUnpaired(find_it->second);

      if (filter_unpaired(find_it->second))
        post_process_unpaired(std::move(find_it->second));
    }

    read_first.erase(find_it);
  } //... while readRegion


  if (read_first.size() > 0)
  {
    for (auto && rec : read_first)
    {
      makeUnpaired(rec.second);

      if (filter_unpaired(rec.second))
        post_process_unpaired(std::move(rec.second));
    }

    read_first.clear();
  }

  // Write remaining reads
  for (auto & rec : read_set)
  {
    long const bin1 = (rec.beginPos - first_pos) / 50l;
    long const bin2 = (rec.pNext - first_pos) / 50l;

    if (bin_counts[bin1] < (opts.SUPER_HI_DEPTH * max_bin_sum) ||
        (hasFlagMultiple(rec) && bin_counts[bin2] < (opts.SUPER_HI_DEPTH * max_bin_sum)))
    {
      writeRecord(bamFileOut, rec);
    }
  }

  // Check if there are any paired reads which their mate was never found
  if (is_paired_reads_with_no_mate)
  {
    BOOST_LOG_TRIVIAL(info) << __HERE__ << " Some reads had flags indicating that they "
                            << "are paired with a mate read in the region but the mate was never found.";
  }
}


String<Triple<CharString, int, int> >
readIntervals(Options const & opts)
{
  // String of intervals to return
  String<Triple<CharString, int, int> > intervalString;
  std::string chrStr;
  Triple<CharString, int, int> chr_start_end;
  std::ifstream intFile(opts.intervalFile.c_str());

  if (intFile.fail())
  {
    BOOST_LOG_TRIVIAL(error) << "Unable to locate interval file at: " << opts.intervalFile << "\n";
    std::exit(1);
  }

  //Read first interval and add to string
  intFile >> chrStr;
  chr_start_end.i1 = chrStr;
  intFile >> chr_start_end.i2;
  --chr_start_end.i2;
  intFile >> chr_start_end.i3;
  --chr_start_end.i3;
  append(intervalString, chr_start_end);

  // If file only contains one interval, return.
  if (intFile.eof())
  {
    intFile.close();
    return intervalString;
  }

  while (!intFile.eof())
  {
    intFile >> chrStr;
    chr_start_end.i1 = chrStr;
    intFile >> chr_start_end.i2;
    --chr_start_end.i2;
    intFile >> chr_start_end.i3;
    --chr_start_end.i3;

    if (intFile.eof())
      break;

    if (chr_start_end.i1 == intervalString[length(intervalString) - 1].i1 &&
        chr_start_end.i2 < intervalString[length(intervalString) - 1].i2)
    {
      BOOST_LOG_TRIVIAL(error) << "The input intervals are not sorted.";
      std::exit(1);
    }

    // If beginning of interval is closer than 2*maxFragLen bases to the previous interval we merge them.
    // Otherwise we cannot ensure sorting of reads.
    if (chr_start_end.i2 - intervalString[length(intervalString) - 1].i3 <= 2 * opts.maxFragLen &&
        chr_start_end.i1 == intervalString[length(intervalString) - 1].i1)
    {
      intervalString[length(intervalString) - 1].i3 = chr_start_end.i3;
    }
    else
    {
      append(intervalString, chr_start_end);
    }
  }

  intFile.close();
  return intervalString;
}


int
main(bamshrink::Options & opts)
{
  if (opts.bamIndex == std::string("<bamPathIn>.[bai,crai]"))
  {
    if (opts.bamPathIn.size() > 5 && std::string(opts.bamPathIn.rbegin(), opts.bamPathIn.rbegin() + 5) == "marc.")
      opts.bamIndex = opts.bamPathIn + std::string(".crai");
    else
      opts.bamIndex = opts.bamPathIn + std::string(".bai");
  }

  String<Triple<CharString, int, int> > intervalString;

  if (opts.intervalFile.size() > 0)
  {
    intervalString = readIntervals(opts);

    if (length(intervalString) == 0)
    {
      BOOST_LOG_TRIVIAL(error) << "The interval file \"" << opts.intervalFile << "\" contained no intervals!";
      return 1;
    }
  }

  if (opts.interval.size() > 0)
  {
    assert(opts.interval.size() > 0);
    auto begin_it = opts.interval.begin();
    auto end_it = opts.interval.end();
    auto find_colon_it = std::find(begin_it, end_it, ':');
    auto find_dash_it = std::find(begin_it, end_it, '-');

    if (find_colon_it == end_it or find_dash_it == end_it)
    {
      BOOST_LOG_TRIVIAL(error) << "Could not parse interval '" << opts.interval << "'";
      return 1;
    }

    Triple<CharString, int, int> interval;
    interval.i1 = std::string(begin_it, find_colon_it).c_str();
    interval.i2 = std::stoi(std::string(find_colon_it + 1, find_dash_it)) - 1;
    interval.i3 = std::stoi(std::string(find_dash_it + 1, end_it)) - 1;
    appendValue(intervalString, interval);
  }

  BamFileIn bamFileIn;
  //
  open(bamFileIn, opts.bamPathIn.c_str());

  BamFileOut bamFileOut(opts.bamPathOut.c_str(), "wb");
  BamAlignmentRecord record;

  if (length(intervalString) == 0)
  {
    BOOST_LOG_TRIVIAL(error) << "[graphtyper::bamshrink] Some intervals are required to extract reads from.";
    std::exit(1);
  }

  // When there is only one contig, remove all other contigs from header to save space
  if (length(intervalString) == 1)
  {
    std::ostringstream ss;
    ss << "@SQ\tSN:" << intervalString[0].i1 << "\t";
    std::string const sq = ss.str();
    long const sq_len = sq.size();
    std::string const hd = "@HD\t";
    std::string const rg = "@RG\t";
    const char * t = bamFileIn.hdr->text;
    const char * t_end = bamFileIn.hdr->text + (bamFileIn.hdr->l_text - 4);
    std::ostringstream new_ss;

    while (t < t_end)
    {
      long line_size = std::distance(t, std::find(t, t + bamFileIn.hdr->l_text, '\n'));

      if (std::equal(hd.begin(), hd.begin() + 4, t) ||
          std::equal(rg.begin(), rg.begin() + 4, t) ||
          (line_size > sq_len && std::equal(sq.begin(), sq.begin() + sq_len, t))
          )
      {
        new_ss << std::string(t, line_size) << '\n';
      }

      t += line_size + 1;
    }

    std::string new_header = new_ss.str();
    bamFileOut.hdr = sam_hdr_parse(new_header.size(), new_header.c_str());
    bamFileOut.hdr->l_text = new_header.size();
    bamFileOut.hdr->text = static_cast<char *>(realloc(bamFileOut.hdr->text, sizeof(char) * new_header.size()));
    strncpy(bamFileOut.hdr->text, new_header.c_str(), new_header.size());
    writeHeader(bamFileOut);
    long read_num {
      0
    };
    qualityFilterSlice2(opts, intervalString[0], bamFileIn, bamFileOut, read_num, true); // is_single_contig
  }
  else
  {
    copyHeader(bamFileOut, bamFileIn);
    writeHeader(bamFileOut);
    long read_num {
      0
    };

    for (unsigned i = 0; i < length(intervalString); ++i)
    {
      qualityFilterSlice2(opts, intervalString[i], bamFileIn, bamFileOut, read_num, false); // is_single_contig
    }
  }

  return 0;
}


} // namespace bamshrink


namespace gyper
{

void
bamshrink(seqan::String<seqan::Triple<seqan::CharString, int, int> > const & intervals,
          std::string const & path_in,
          std::string const & path_out,
          double const avg_cov_by_readlen,
          const char * reference_genome)
{
  if (seqan::length(intervals) == 0)
  {
    BOOST_LOG_TRIVIAL(warning) << "No intervals to read regions from. Aborting bamshrink.";
    return;
  }

  bamshrink::Options opts;
  BamFileIn bamFileIn;
  open(bamFileIn, path_in.c_str(), reference_genome);

  BamFileOut bamFileOut(path_out.c_str(), "wb");

  if (path_in.size() > 5 && std::string(path_in.rbegin(), path_in.rbegin() + 5) == "marc.")
    opts.bamIndex = path_in + ".crai";
  else
    opts.bamIndex = path_in + ".bai";

  if (avg_cov_by_readlen > 0.0)
    opts.avgCovByReadLen = avg_cov_by_readlen;
  else
    opts.SUPER_HI_DEPTH = 1000; // Do not apply super hi depth if coverage is unknown

  auto const & copts = *(Options::const_instance());
  opts.maxFragLen = copts.bamshrink_max_fraglen;
  opts.minNumMatching = copts.bamshrink_min_matching;
  opts.is_filtering_mapq0 = !copts.bamshrink_is_not_filtering_mapq0;
  opts.minReadLen = copts.bamshrink_min_readlen;
  opts.minReadLenMapQ0 = copts.bamshrink_min_readlen_low_mapq;
  opts.minUnpairedReadLen = copts.bamshrink_min_unpair_readlen;
  opts.as_filter_threshold = copts.bamshrink_as_filter_threshold;
  opts.no_filter_on_coverage = copts.no_filter_on_coverage;

  // When there is only one contig, remove all other contigs from header to save space
  if (seqan::length(intervals) == 1)
  {
    auto const & interval = intervals[0];
    std::ostringstream ss;
    ss << "@SQ\tSN:" << interval.i1 << "\t";
    std::string const sq = ss.str();
    long const sq_len = sq.size();
    std::string const hd = "@HD\t";
    std::string const rg = "@RG\t";
    const char * t = bamFileIn.hdr->text;
    const char * t_end = bamFileIn.hdr->text + bamFileIn.hdr->l_text;
    std::ostringstream new_ss;

    while (t <= t_end)
    {
      long const line_size = std::distance(t, std::find(t, t_end, '\n'));

      if (line_size > 4 &&
          (std::equal(hd.begin(), hd.begin() + 4, t) ||
           std::equal(rg.begin(), rg.begin() + 4, t) ||
           (line_size > sq_len && std::equal(sq.begin(), sq.begin() + sq_len, t))))
      {
        new_ss << std::string(t, line_size) << '\n';
      }

      t += line_size + 1;
    }

    std::string new_header = new_ss.str();
    bamFileOut.hdr = sam_hdr_parse(new_header.size(), new_header.c_str());
    bamFileOut.hdr->l_text = new_header.size();
    bamFileOut.hdr->text = static_cast<char *>(realloc(bamFileOut.hdr->text, sizeof(char) * new_header.size()));
    strncpy(bamFileOut.hdr->text, new_header.c_str(), new_header.size());
    writeHeader(bamFileOut);
    long read_num {
      0
    };
    bamshrink::qualityFilterSlice2(opts, interval, bamFileIn, bamFileOut, read_num, true); // is_single_contig
  }
  else
  {
    copyHeader(bamFileOut, bamFileIn);
    writeHeader(bamFileOut);
    long read_num {
      0
    };

    for (auto const & interval : intervals)
    {
      bamshrink::qualityFilterSlice2(opts, interval, bamFileIn, bamFileOut, read_num, false); // is_single_contig
    }
  }
}


void
bamshrink(std::string const & chrom,
          int begin,
          int end,
          std::string const & path_in,
          std::string const & path_out,
          double const avg_cov_by_readlen,
          std::string const & ref_fn)
{
  seqan::String<seqan::Triple<seqan::CharString, int, int> > intervals;
  seqan::Triple<seqan::CharString, int, int> interval;
  interval.i1 = chrom.c_str();
  interval.i2 = begin;
  interval.i3 = end;
  seqan::appendValue(intervals, interval);
  BOOST_LOG_TRIVIAL(debug) << "Bamshrink is copying file " << path_in;
  bamshrink(intervals, path_in, path_out, avg_cov_by_readlen, ref_fn.c_str());
}


void
bamshrink_multi(std::string const & interval_fn,
                std::string const & path_in,
                std::string const & path_out,
                double const avg_cov_by_readlen,
                std::string const & ref_fn)
{
  bamshrink::Options opts;
  opts.intervalFile = interval_fn;
  seqan::String<seqan::Triple<seqan::CharString, int, int> > intervals = readIntervals(opts);
  bamshrink(intervals, path_in, path_out, avg_cov_by_readlen, ref_fn.c_str());
}


} // namespace gyper
