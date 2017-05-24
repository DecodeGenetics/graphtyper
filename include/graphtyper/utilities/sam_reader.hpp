#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include <graphtyper/utilities/hash_seqan.hpp>


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
  std::size_t p = 0; /* Current pos in the region */

private:
  seqan::HtsFileIn hts_file;
  bool second_file = false;
  std::vector<std::string> regions;
  TReadsFirst reads_first;
  TReads unpaired_reads;

};

} // namespace gyper
