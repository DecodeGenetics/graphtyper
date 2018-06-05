#include <iostream>
#include <cstdlib>
#include <cassert>

#include <boost/log/trivial.hpp>
#include <boost/algorithm/string.hpp>

#include <seqan/align.h>
#include <seqan/stream.h>
#include <seqan/hts_io.h>

#include <graphtyper/constants.hpp>
#include <graphtyper/utilities/adapter_removal.hpp>


namespace
{

using TEraseTuple = std::tuple<uint16_t, uint16_t, uint16_t, uint16_t>;

} // anon namespace


namespace gyper
{

template <typename Type>
AdapterRemoval<Type>::AdapterRemoval()
  : adapters(0)
{}


template <typename Type>
template <typename TSeq>
void
AdapterRemoval<Type>::remove(TEraseTuple const & bases_to_keep, TSeq & read1, TSeq & read2)
{
  seqan::erase(read1, std::get<0>(bases_to_keep), std::get<1>(bases_to_keep));
  seqan::erase(read2, std::get<2>(bases_to_keep), std::get<3>(bases_to_keep));
}


template <typename Type>
TEraseTuple
AdapterRemoval<Type>::remove_adapters_from_read(seqan::Dna5String & read1, seqan::Dna5String & read2)
{
  seqan::reverseComplement(read2);
  auto results = remove_adapters_from_read_read2_complemented(read1, read2);
  seqan::reverseComplement(read2);

  // Reverse the bases to keep
  std::size_t const read_size = seqan::length(read2);
  uint16_t const tmp = std::get<2>(results);
  std::get<2>(results) = read_size - std::get<3>(results);
  std::get<3>(results) = read_size - tmp;
  return results;
}


template <typename Type>
TEraseTuple
AdapterRemoval<Type>::remove_adapters_from_read_read2_complemented(seqan::Dna5String const & read1,
                                                                   seqan::Dna5String const & read2
                                                                   )
{
  typedef seqan::Dna5String TSequence;
  typedef seqan::Align<TSequence, seqan::ArrayGaps> TAlign;

  static const unsigned ADAPTER_MISMATCHES_ALLOWED = 2;
  static const unsigned READ_MISMATCHES_ALLOWED = 5;
  static const unsigned MIN_PARTIAL_MATCHES = 10;

  bool found_entire_adapter_on_read1 = false;
  bool found_entire_adapter_on_read2 = false;

  // Best alignments (if found)
  TAlign best_align_case1; seqan::resize(seqan::rows(best_align_case1), 2);
  TAlign best_align_case2; seqan::resize(seqan::rows(best_align_case2), 2);

  // The scores and lengths of the best alignments
  int best_score_case1 = 10;
  unsigned best_adapter_length_case1 = 100;
  int best_score_case2 = 10;
  unsigned best_adapter_length_case2 = 100;

  // std::cout << "read1 = " << read1 << std::endl;
  // std::cout << "read2 = " << read2 << std::endl;
  std::vector<seqan::Dna5String> best_adapter2_from_1;
  seqan::Dna5String best_adapter2_from_2;

  TAlign align;
  resize(rows(align), 2);

  // unsigned best_first_i = 0; // For counter
  // unsigned best_second_i = 0; // For counter

  // Check prefix on first read
  bool found_prefix_on_first_read = false;
  for (unsigned i = 0; i < adapter_first_prefix.size(); ++i)
  {
    {
      seqan::assignSource(seqan::row(align, 0), read1);
      seqan::assignSource(seqan::row(align, 1), adapter_first_prefix[i]);

      int score = globalAlignment(align,
                                  seqan::Score<int, seqan::Simple>(1, 0, -1),
                                  seqan::AlignConfig<true, true, false, true>()
                                  );

      if (score >= static_cast<int>(seqan::length(adapter_first_prefix[i]) - ADAPTER_MISMATCHES_ALLOWED))
      {
        found_prefix_on_first_read = true;
        break;
      }
    }
  }

  // Case 1: Read 1
  if (found_prefix_on_first_read)
  {
    // std::cout << "Check case 1: Check if read 1 covers adapter 1" << std::endl;
    for (unsigned i = 0; i < adapters.size(); ++i)
    {
      // Check case 1: Check if read 1 covers adapter 1
      // ----------read1-------->
      //    %%%Adapter1%%%>
      {
        seqan::assignSource(seqan::row(align, 0), read1);
        seqan::assignSource(seqan::row(align, 1), adapters[i].first);

        int score = globalAlignment(align,
                                    seqan::Score<int, seqan::Simple>(1, 0, -1),
                                    seqan::AlignConfig<true, true, false, true>()
                                    );

        if (score > best_score_case1)
        {
          // std::cout << "possible adapter for read 1 = " << adapt_it->first << std::endl;
          best_score_case1 = score;
          best_align_case1 = align;
          best_adapter_length_case1 = seqan::length(adapters[i].first);

          best_adapter2_from_1 = adapters[i].second;
        }
      }
    }

    if (static_cast<unsigned>(best_score_case1) >= best_adapter_length_case1 - ADAPTER_MISMATCHES_ALLOWED)
    {
      // std::cout << "Entire match of read 1" << std::endl;
      found_entire_adapter_on_read1 = true;
    }
  }

  // Case 1: Read 2
  if (best_adapter2_from_1.size() != 0)
  {
    {
      // Check case 2: Check if read 2 covers adapter 2
      //     %%%Adapter2%%%>
      // ----------read2-------->
      for (unsigned i = 0; i < best_adapter2_from_1.size(); ++i)
      {
        seqan::assignSource(seqan::row(align, 0), best_adapter2_from_1[i]);
        seqan::assignSource(seqan::row(align, 1), read2);

        int score = seqan::globalAlignment(align,
                                           seqan::Score<int, seqan::Simple>(1, 0, -1),
                                           seqan::AlignConfig<false, true, true, true>()
                                           );

        if (score > best_score_case2)
        {
          // std::cout << "possible adapter for read 2 = " << best_adapter2_from_1[i] << std::endl;
          best_score_case2 = score;
          best_align_case2 = align;
          best_adapter_length_case2 = seqan::length(best_adapter2_from_1[i]);
          // best_second_i = i;
        }
      }
    }

    if (static_cast<unsigned>(best_score_case2) >= best_adapter_length_case2 - ADAPTER_MISMATCHES_ALLOWED)
    {
      found_entire_adapter_on_read2 = true;
    }
  }

  if (found_entire_adapter_on_read1 and found_entire_adapter_on_read2)
  {
    // BOOST_LOG_TRIVIAL(debug) << "Case 1: Erase both adapters and everything else afterwards.";
    return std::make_tuple<uint16_t, uint16_t, uint16_t, uint16_t>
           (
      seqan::toViewPosition(row(best_align_case1, 1), 0),
      seqan::length(read1),
      0u,
      uint16_t(std::min(static_cast<uint16_t>(seqan::length(read2)),
                        static_cast<uint16_t>(seqan::toViewPosition(row(best_align_case2, 0), best_adapter_length_case2))
                        )
               )
           );
  }
  else if (found_entire_adapter_on_read1)
  {
    // Check if we can find a partial alignment to read 2, otherwise we don't remove any adaptors
    // BOOST_LOG_TRIVIAL(debug) << "Case 2: Check for a partial match of read 2.";
    int best_score = 9;
    unsigned best_adapter_length = 100;
    TAlign best_align;

    for (auto adapt_it_second = best_adapter2_from_1.cbegin(); adapt_it_second != best_adapter2_from_1.cend(); ++adapt_it_second)
    {
      // Check if read 2 partially matches adapter 2
      //               %%%Adapter2%%%>
      // ----------read2-------->
      {
        seqan::assignSource(seqan::row(align, 0), *adapt_it_second);
        seqan::assignSource(seqan::row(align, 1), read2);

        int score = globalAlignment(align,
                                    seqan::Score<int, seqan::Simple>(1, -5, -100),
                                    seqan::AlignConfig<true, false, true, false>()
                                    );

        if (score > best_score)
        {
          best_score = score;
          best_align = align;
          best_adapter_length = seqan::length(*adapt_it_second);
        }
      }
    }

    if (static_cast<unsigned>(best_score) >= MIN_PARTIAL_MATCHES)
    {
      // BOOST_LOG_TRIVIAL(debug) << "Case 2: Erased a partial match on read 2 and entire match on read 1";
      return std::make_tuple<uint16_t, uint16_t, uint16_t, uint16_t>
             (
        seqan::toViewPosition(row(best_align_case1, 1), 0),
        seqan::length(read1),
        0u,
        seqan::toViewPosition(row(best_align, 0), best_adapter_length) - seqan::toViewPosition(row(best_align, 1), 0)
             );
    }

  }
  else if (found_entire_adapter_on_read2)
  {
    // BOOST_LOG_TRIVIAL(debug) << "Case 3: Check for a partial match of read 1.";
    int best_score = 9;
    // unsigned best_adapter_length = 100;
    TAlign best_align;

    for (auto adapt_it = adapters.cbegin(); adapt_it != adapters.cend(); ++adapt_it)
    {
      // Check if read 1 covers adapter 1
      // ----------read1-------->
      //               %%%Adapter1%%%>
      {
        seqan::assignSource(seqan::row(align, 0), read1);
        seqan::assignSource(seqan::row(align, 1), adapt_it->first);

        int score = globalAlignment(align,
                                    seqan::Score<int, seqan::Simple>(1, -5, -100),
                                    seqan::AlignConfig<true, false, true, false>()
                                    );

        if (score > best_score)
        {
          best_score = score;
          best_align = align;
          // best_adapter_length = seqan::length(adapt_it->first);
        }
      }
    }

    if (static_cast<unsigned>(best_score) >= MIN_PARTIAL_MATCHES)
    {
      // BOOST_LOG_TRIVIAL(debug) << "Case 3: Erased for a partial match of read 1";
      return std::make_tuple<uint16_t, uint16_t, uint16_t, uint16_t>
             (
        seqan::toViewPosition(row(best_align, 1), 0),
        seqan::length(read1),
        0u,
        std::min(static_cast<std::size_t>(seqan::length(read2)),
                 static_cast<std::size_t>(seqan::toViewPosition(row(best_align_case2, 0), best_adapter_length_case2))
                 )
             );
    }
  }
  else
  {
    // BOOST_LOG_TRIVIAL(debug) << "Case 4: Check for two partial matches on both reads";
    int best_score = std::min(seqan::length(read1), seqan::length(read2)) - READ_MISMATCHES_ALLOWED;
    seqan::Dna5String best_adapter1(read1);
    seqan::Dna5String best_adapter2(read2);
    TAlign best_align;

    seqan::Dna5String old_prefix = "";

    for (auto adapt_it = adapters.cbegin(); adapt_it != adapters.cend(); ++adapt_it)
    {
      seqan::Dna5String prefix = "";
      seqan::resize(prefix, 8);
      seqan::arrayCopyForward(seqan::begin(adapt_it->first), seqan::begin(adapt_it->first) + 8, seqan::begin(prefix));

      if (prefix == old_prefix)
      {
        continue;
      }

      for (auto adapt2_it = adapt_it->second.begin(); adapt2_it != adapt_it->second.end(); ++adapt2_it)
      {
        // Check if read 1 covers adapter 1
        // ----------read1-------->%%%Adapter2%%%
        //               %%%Adapter1%%%>------read2------->
        {
          seqan::Dna5String sequence1(*adapt2_it);
          seqan::append(sequence1, read1);
          seqan::Dna5String sequence2(read2);
          seqan::append(sequence2, adapt_it->first);

          seqan::assignSource(seqan::row(align, 0), sequence1);
          seqan::assignSource(seqan::row(align, 1), sequence2);

          int score = seqan::globalAlignment(align,
                                             seqan::Score<int, seqan::Simple>(1, -1, -5),
                                             seqan::AlignConfig<true, false, true, false>()
                                             );

          if (score > best_score)
          {
            best_score = score;
            best_align = align;
            best_adapter1 = adapt_it->first;
            best_adapter2 = *adapt2_it;
          }
          else if (score < static_cast<int>(std::min(seqan::length(read1), seqan::length(read2)) - READ_MISMATCHES_ALLOWED - 16))
          {
            old_prefix = prefix;
            break;
          }
        }
      }
    }

    if (static_cast<unsigned>(best_score) > std::min(seqan::length(read1), seqan::length(read2)) - READ_MISMATCHES_ALLOWED)
    {
      uint16_t const read1_bases_to_erase_from = seqan::length(row(best_align, 1)) - seqan::length(best_adapter1) - seqan::length(best_adapter2);

      if (read1_bases_to_erase_from < seqan::length(read1))
      {
        return std::make_tuple<uint16_t, uint16_t, uint16_t, uint16_t>
               (
          uint16_t(read1_bases_to_erase_from),
          seqan::length(read1),
          0u,
          seqan::length(read2) < read1_bases_to_erase_from ? 0u : seqan::length(read2) - read1_bases_to_erase_from
               );
      }
    }
  }

  return std::make_tuple(0u, 0u, 0u, 0u);
}


/***************************
 * ILLUMINA SPECIALIZATION *
 ***************************/
// ATCGGAAGAGCACACGTCTGAACTCCAGTCAC CTGAAGCT ATCTCGTATGCCGTCTTCTGCTTG
//
static const seqan::Dna5String hiseq_universal_string    = "ATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT";
static const seqan::Dna5String hiseq_universal_string_RC  = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGAT";
static const seqan::Dna5String hiseqx_universal_string    = "ATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTGCCTCTATGTGTAGATCTCGGTGGTCGCCGTATCATT";
static const seqan::Dna5String hiseqx_universal_string_RC = "AATGATACGGCGACCACCGAGATCTACACATAGAGGCACACTCTTTCCCTACACGACGCTCTTCCGAT";
// AATGATACGGCGACCACCGAGATCTACAC GGCTCTGA ACACTCTTTCCCTACACGACGCTCTTCCGATCT

static const seqan::Dna5String hiseq_prefix = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
static const seqan::Dna5String hiseqx_prefix = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
static const std::vector<seqan::Dna5String> barcode_reads2 =
{
  "AATGATACGGCGACCACCGAGATCTACAC" "TATAGCCT" "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
  "AATGATACGGCGACCACCGAGATCTACAC" "ATAGAGGC" "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
  "AATGATACGGCGACCACCGAGATCTACAC" "CCTATCCT" "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
  "AATGATACGGCGACCACCGAGATCTACAC" "GGCTCTGA" "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
  "AATGATACGGCGACCACCGAGATCTACAC" "AGGCGAAG" "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
  "AATGATACGGCGACCACCGAGATCTACAC" "TAATCTTA" "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
  "AATGATACGGCGACCACCGAGATCTACAC" "CAGGACGT" "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
  "AATGATACGGCGACCACCGAGATCTACAC" "GTACTGAC" "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
  hiseqx_universal_string_RC
};

template <>
AdapterRemoval<Illumina>::AdapterRemoval()
  : adapters
  (
  {
    // HiSeq
    {"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "ATCACG"   "ATCTCGTATGCCGTCTTCTGCTTG", {hiseq_universal_string_RC}},
    {"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "CGATGT"   "ATCTCGTATGCCGTCTTCTGCTTG", {hiseq_universal_string_RC}},
    {"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "TTAGGC"   "ATCTCGTATGCCGTCTTCTGCTTG", {hiseq_universal_string_RC}},
    {"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "TGACCA"   "ATCTCGTATGCCGTCTTCTGCTTG", {hiseq_universal_string_RC}},
    {"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "ACAGTG"   "ATCTCGTATGCCGTCTTCTGCTTG", {hiseq_universal_string_RC}},
    {"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "GCCAAT"   "ATCTCGTATGCCGTCTTCTGCTTG", {hiseq_universal_string_RC}},
    {"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "CAGATC"   "ATCTCGTATGCCGTCTTCTGCTTG", {hiseq_universal_string_RC}},
    {"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "ACTTGA"   "ATCTCGTATGCCGTCTTCTGCTTG", {hiseq_universal_string_RC}},
    {"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "GATCAG"   "ATCTCGTATGCCGTCTTCTGCTTG", {hiseq_universal_string_RC}},
    {"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "TAGCTT"   "ATCTCGTATGCCGTCTTCTGCTTG", {hiseq_universal_string_RC}},
    {"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "GGCTAC"   "ATCTCGTATGCCGTCTTCTGCTTG", {hiseq_universal_string_RC}},
    {"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "CTTGTA"   "ATCTCGTATGCCGTCTTCTGCTTG", {hiseq_universal_string_RC}},
    {"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "AGTCAACA" "ATCTCGTATGCCGTCTTCTGCTTG", {hiseq_universal_string_RC}},
    {"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "AGTTCCGT" "ATCTCGTATGCCGTCTTCTGCTTG", {hiseq_universal_string_RC}},
    {"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "ATGTCAGA" "ATCTCGTATGCCGTCTTCTGCTTG", {hiseq_universal_string_RC}},
    {"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "CCGTCCCG" "ATCTCGTATGCCGTCTTCTGCTTG", {hiseq_universal_string_RC}},
    {"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "GTCCGCAC" "ATCTCGTATGCCGTCTTCTGCTTG", {hiseq_universal_string_RC}},
    {"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "GTGAAACG" "ATCTCGTATGCCGTCTTCTGCTTG", {hiseq_universal_string_RC}},
    {"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "GTGGCCTT" "ATCTCGTATGCCGTCTTCTGCTTG", {hiseq_universal_string_RC}},
    {"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "GTTTCGGA" "ATCTCGTATGCCGTCTTCTGCTTG", {hiseq_universal_string_RC}},
    {"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "CGTACGTA" "ATCTCGTATGCCGTCTTCTGCTTG", {hiseq_universal_string_RC}},
    {"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "GAGTGGAT" "ATCTCGTATGCCGTCTTCTGCTTG", {hiseq_universal_string_RC}},
    {"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "ACTGATAT" "ATCTCGTATGCCGTCTTCTGCTTG", {hiseq_universal_string_RC}},
    {"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "ATTCCTTT" "ATCTCGTATGCCGTCTTCTGCTTG", {hiseq_universal_string_RC}},

    // HiSeqX
    {"GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "ATTACTCG" "ATCTCGTATGCCGTCTTCTGCTTG", barcode_reads2},
    {"GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "TCCGGAGA" "ATCTCGTATGCCGTCTTCTGCTTG", barcode_reads2},
    {"GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "CGCTCATT" "ATCTCGTATGCCGTCTTCTGCTTG", barcode_reads2},
    {"GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "GAGATTCC" "ATCTCGTATGCCGTCTTCTGCTTG", barcode_reads2},
    {"GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "ATTCAGAA" "ATCTCGTATGCCGTCTTCTGCTTG", barcode_reads2},
    {"GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "GAATTCGT" "ATCTCGTATGCCGTCTTCTGCTTG", barcode_reads2},
    {"GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "CTGAAGCT" "ATCTCGTATGCCGTCTTCTGCTTG", barcode_reads2},
    {"GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "TAATGCGC" "ATCTCGTATGCCGTCTTCTGCTTG", barcode_reads2},
    {"GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "CGGCTATG" "ATCTCGTATGCCGTCTTCTGCTTG", barcode_reads2},
    {"GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "TCCGCGAA" "ATCTCGTATGCCGTCTTCTGCTTG", barcode_reads2},
    {"GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "TCTCGCGC" "ATCTCGTATGCCGTCTTCTGCTTG", barcode_reads2},
    {"GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "AGCGATAG" "ATCTCGTATGCCGTCTTCTGCTTG", barcode_reads2}

    // TruSeq universal = AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
    // {"A", "T"}
  }


  )
  , adapter_first_prefix({hiseq_prefix, hiseqx_prefix})
{}


template class AdapterRemoval<Illumina>;

template void AdapterRemoval<Illumina>::remove<seqan::Dna5String>(TEraseTuple const & bases_to_keep,
                                                                  seqan::Dna5String & read1,
                                                                  seqan::Dna5String & read2
                                                                  );
template void AdapterRemoval<Illumina>::remove<seqan::CharString>(TEraseTuple const & bases_to_keep,
                                                                  seqan::CharString & read1,
                                                                  seqan::CharString & read2
                                                                  );
template void AdapterRemoval<Illumina>::remove<seqan::DnaString>(TEraseTuple const & bases_to_keep,
                                                                 seqan::DnaString & read1,
                                                                 seqan::DnaString & read2
                                                                 );
template void AdapterRemoval<Illumina>::remove<seqan::IupacString>(TEraseTuple const & bases_to_keep,
                                                                   seqan::IupacString & read1,
                                                                   seqan::IupacString & read2
                                                                   );

void
remove_adapters_from_reads(TReads & reads)
{
  gyper::AdapterRemoval<gyper::Illumina> ar;

  for (auto read_it = reads.begin(); read_it != reads.end(); ++read_it)
  {
    auto bases_to_erase = ar.remove_adapters_from_read_read2_complemented(read_it->first.seq, read_it->second.seq);
    std::size_t const SIZE_BEFORE_1 = seqan::length(read_it->first.seq);
    std::size_t const SIZE_BEFORE_2 = seqan::length(read_it->second.seq);
    ar.remove(bases_to_erase, read_it->first.seq, read_it->second.seq);
    ar.remove(bases_to_erase, read_it->first.qual, read_it->second.qual);

    if (seqan::length(read_it->first.seq) != SIZE_BEFORE_1)
      seqan::clear(read_it->first.cigar);

    if (seqan::length(read_it->second.seq) != SIZE_BEFORE_2)
      seqan::clear(read_it->second.cigar);
  }

  // Delete all short read pairs
  reads.erase(
    std::remove_if(
      reads.begin(),
      reads.end(),
      [](TReadPair const & rp)
    {
      return seqan::length(rp.first.seq) < (2 * K - 1) || seqan::length(rp.second.seq) < (2 * K - 1);
    }),
    reads.end()
    );
}


// Remove adapters from Bam records
void
remove_adapters(TReads & reads)
{
  gyper::AdapterRemoval<gyper::Illumina> ar;
  unsigned i = 0;

  for (auto read_it = reads.begin(); read_it != reads.end(); ++read_it, ++i)
  {
    seqan::BamTagsDict tags1(read_it->first.tags);

    // Check if adapter information is already known
    {
      unsigned tagIdx = 0;

      if (seqan::findTagKey(tagIdx, tags1, ADAPTER_TAG))
      {
        // The read pair has been checked for adapters. Leave it as it is.
        continue;
      }
    }

    auto bases_to_erase = ar.remove_adapters_from_read_read2_complemented(read_it->first.seq, read_it->second.seq);
    // BOOST_LOG_TRIVIAL(debug) << "Case 4: " << std::get<0>(bases_to_erase) << "-" << std::get<1>(bases_to_erase) << " "
    //                                        << std::get<2>(bases_to_erase) << "-" << std::get<3>(bases_to_erase);

    // lambda function which adds the the adapter removal tag to the
    auto set_tags = [](seqan::BamTagsDict & tags, uint16_t start, uint16_t end)
                    {
                      if (start == 0)
                      {
                        setTagValue(tags, ADAPTER_TAG, static_cast<short>(end));
                        // BOOST_LOG_TRIVIAL(debug) << "AD = " << start << " " << end;
                      }
                      else
                      {
                        setTagValue(tags, ADAPTER_TAG, static_cast<short>(start) - static_cast<short>(end));
                      }
                    };

    seqan::BamTagsDict tags2(read_it->second.tags);

    set_tags(tags1, std::get<0>(bases_to_erase), std::get<1>(bases_to_erase));
    set_tags(tags2, std::get<2>(bases_to_erase), std::get<3>(bases_to_erase));

    // if (std::get<1>(bases_to_erase) > std::get<0>(bases_to_erase))
    // {
    //   // Adapter at front
    //   assert (std::get<0>(bases_to_erase) == 0);
    //   setTagValue(tags1, ADAPTER_TAG, static_cast<short>(std::get<1>(bases_to_erase)));
    // }
    // else
    // {
    //   // Adapter at back (or no adapter if std::get<1>(bases_to_erase) == std::get<0>(bases_to_erase))
    //   setTagValue(tags1, ADAPTER_TAG, static_cast<short>(std::get<1>(bases_to_erase)) - static_cast<short>(std::get<0>(bases_to_erase)));
    // }
    //
    //
    //
    // if (std::get<3>(bases_to_erase) > std::get<2>(bases_to_erase))
    // {
    //   // Adapter at front
    //   assert (std::get<2>(bases_to_erase) == 0);
    //   setTagValue(tags2, ADAPTER_TAG, static_cast<short>(std::get<3>(bases_to_erase)));
    // }
    // else
    // {
    //   // Adapter at back (or no adapter if std::get<3>(bases_to_erase) == std::get<2>(bases_to_erase))
    //   setTagValue(tags2, ADAPTER_TAG, static_cast<short>(std::get<3>(bases_to_erase)) - static_cast<short>(std::get<2>(bases_to_erase)));
    // }

    // ar.remove(bases_to_erase, read_it->first.seq, read_it->second.seq);
    // ar.remove(bases_to_erase, read_it->first.qual, read_it->second.qual);
    // if (std::get<1>(bases_to_erase) - std::get<0>(bases_to_erase) == 0)
    // {
    //   BOOST_LOG_TRIVIAL(debug) << seqan::length(read_it->first.seq);
    // }
  }

  BOOST_LOG_TRIVIAL(info) << "[adapter_removal] I have checked " << reads.size() << " read pairs.";
}


} // namespace gyper
