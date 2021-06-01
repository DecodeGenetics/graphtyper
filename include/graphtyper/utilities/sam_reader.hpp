#pragma once

#include <string>        // std::String
#include <unordered_map> // std::unordered_map
#include <utility>       // std::pair
#include <vector>        // std::vector

#include <seqan/basic.h>
#include <seqan/hts_io.h> // seqan::HtsFileIn, seqan::BamAlignmentRecord
#include <seqan/sequence.h>

#include <graphtyper/utilities/hash_seqan.hpp> // For std::unordered_map<seqan::String<char>, seqan::BamAlignmentRecord>

namespace gyper
{
using TReadPair = std::pair<seqan::BamAlignmentRecord, seqan::BamAlignmentRecord>;
using TReads = std::vector<TReadPair>;
using TReadsFirst = std::unordered_map<seqan::String<char>, seqan::BamAlignmentRecord>;

class SamReader
{
public:
  SamReader(std::string const & hts_path, std::vector<std::string> const & _regions);
  TReads read_N_reads(std::size_t const N);
  void insert_reads(TReads & reads, seqan::BamAlignmentRecord && record);

  std::size_t r = 0; /* Current region index */

private:
  seqan::HtsFileIn hts_file;
  std::vector<std::string> regions;

  TReadsFirst reads_first;
  TReads unpaired_reads;
};

} // namespace gyper
