#pragma once

#include <vector>
#include <string>
#include <tuple>


#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <graphtyper/utilities/sam_reader.hpp>

namespace gyper
{

// constants
const char ADAPTER_TAG[3] = "AD";


template <typename Type>
class AdapterRemoval
{
public:
  AdapterRemoval();

  /*************
   * MODIFIERS *
   *************/
  // void add_adapters(seqan::Dna5String && adapter1, seqan::Dna5String && adapter2);
  std::tuple<uint16_t, uint16_t, uint16_t, uint16_t>
  remove_adapters_from_read(seqan::Dna5String & read1, seqan::Dna5String & read2);

  std::tuple<uint16_t, uint16_t, uint16_t, uint16_t>
  remove_adapters_from_read_read2_complemented(seqan::Dna5String const & read1, seqan::Dna5String const & read2);

  template <typename TSeq>
  void
  remove(std::tuple<uint16_t, uint16_t, uint16_t, uint16_t> const & bases_to_keep, TSeq & read1, TSeq & read2);


private:
  std::vector<std::pair<seqan::Dna5String, std::vector<seqan::Dna5String> > > adapters;
  std::vector<seqan::Dna5String> adapter_first_prefix;
  // std::vector<std::pair<unsigned, std::vector<unsigned> > > adapter_counter;
};

struct Illumina {};

using IlluminaAdapter = AdapterRemoval<Illumina>;

void remove_adapters_from_reads(TReads & reads);
void remove_adapters(TReads & reads);

}
