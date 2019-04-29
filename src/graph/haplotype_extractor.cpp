#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/graph_msa.h>
#include <seqan/score.h>
#include <seqan/stream.h>

#include <boost/algorithm/string.hpp>
#include <boost/log/trivial.hpp>

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/graph.hpp>
#include <graphtyper/graph/graph_serialization.hpp>
#include <graphtyper/graph/haplotype_calls.hpp>
#include <graphtyper/graph/haplotype_extractor.hpp>
#include <graphtyper/typer/vcf_writer.hpp>
#include <graphtyper/typer/vcf.hpp>
#include <graphtyper/typer/variant.hpp>
#include <graphtyper/utilities/sequence_operations.hpp> // Remove prefix and suffix
#include <graphtyper/utilities/options.hpp>
#include <graphtyper/utilities/graph_help_functions.hpp>
#include <graphtyper/utilities/vcf_help_functions.hpp>
#include <graphtyper/utilities/type_conversions.hpp>


namespace
{

std::vector<std::vector<uint32_t> >
read_haplotype_calls_from_file(std::string const haps_path,
                               std::vector<std::vector<gyper::Genotype> > & gts
  )
{
  using namespace gyper;

  std::map<uint32_t, std::vector<uint32_t> > gt2hapcalls;

  // File with a list of haplotype files
  std::ifstream haps_file(haps_path.c_str());

  // Make sure the file could be opened
  if (not haps_file.is_open())
  {
    BOOST_LOG_TRIVIAL(error) << "[graphtyper::haplotype_extractor] Unable to open file '" << haps_path << "'.";
    std::exit(1);
  }

  // Each line in the haps file contains a location of a haplotype file
  std::string line;

  // Lambda function which processes each line of haps_file
  auto process_line =
    [&gts, &gt2hapcalls](std::string const & line, bool const FIRST_LINE)
    {
      THapCalls new_calls = load_calls(line.c_str());

      for (unsigned i = 0; i < new_calls.size(); ++i)
      {
        assert(not new_calls[i].gts.empty());

        // The first genotype in a haplotype is its ID
        uint32_t const gt_id = new_calls[i].gts[0].id;

        if (FIRST_LINE)
        {
          // Initialize all vectors with 0 (i.e. the reference haplotype).
          gt2hapcalls[gt_id] = {0u};
          gts.push_back(new_calls[i].gts);
        }

        //assert(new_calls[i].first.size() >= 2);

        for (auto const & hap_call : new_calls[i].calls)
          gt2hapcalls.at(gt_id).push_back(hap_call);
      }
    };

  // Process the first line seperately
  if (std::getline(haps_file, line))
  {
    process_line(line, true /*first line*/);
  }
  else
  {
    BOOST_LOG_TRIVIAL(info) << "[graphtyper::haplotype_extractor] No haplotype calls were extracted";
    std::exit(0);
  }

  uint32_t num_files = 1;

  // Process the rest of the lines
  while (std::getline(haps_file, line))
  {
    process_line(line, false /*first line*/);
    ++num_files;
  }

  // Close the file after all lines have been processed
  haps_file.close();

  BOOST_LOG_TRIVIAL(info) << "[graphtyper::haplotype_extractor] Processed all " << num_files << " files.";

  // Create the hap calls
  std::vector<std::vector<uint32_t> > hap_calls;

  for (auto map_it = gt2hapcalls.cbegin(); map_it != gt2hapcalls.cend(); ++map_it)
    hap_calls.push_back(map_it->second);

  assert(gt2hapcalls.size() == gts.size());
  assert(hap_calls.size() == gts.size());

  return hap_calls;
}


} // anon namespace


namespace gyper
{

/**
 * \brief * Aligns reference sequence to read and gets a gapped sequence of both
 * \param gapped_ref Output gapped reference sequence
 * \param gapped_alt Output gapped alternative sequence
 * \param ref Input reference sequence.
 * \param read Input read sequence.
 * \return True if we should use this alignment, False if we should throw it away
 */
bool
get_capped_strings(std::string & gapped_ref,
                   std::string & gapped_alt,
                   std::vector<char> const & ref,
                   seqan::Dna5String const & read
  )
{
  seqan::Align<seqan::Dna5String> alignment;
  seqan::resize(seqan::rows(alignment), 2);
  seqan::assignSource(seqan::row(alignment, 0), ref);
  seqan::assignSource(seqan::row(alignment, 1), read);

  {
    int constexpr SCORE_MATCH = 2;
    int score = seqan::globalAlignment(alignment,
                                       seqan::Score<int, seqan::Simple>(SCORE_MATCH /*match*/,
                                                                        -3 /*mismatch*/,
                                                                        -1 /*gap extend*/,
                                                                        -6 /*gap open*/),
                                       seqan::AlignConfig<true, false, false, true>(),
                                       -10 /* lower diagonal */,
                                       100 /* upper diagonal */
    );

    if (score == SCORE_MATCH * static_cast<int>(seqan::length(read)))
      return false; // Perfect match
    else if (score < 42)
      return false; // Very poor alignment
  }

  {
    std::ostringstream ref_ss;
    ref_ss << row(alignment, 0);
    gapped_ref = ref_ss.str();
  }

  {
    std::ostringstream alt_ss;
    alt_ss << row(alignment, 1);
    gapped_alt = alt_ss.str();
  }

  return true;
}


Variant
make_variant_of_gapped_strings(std::string & gapped_ref,
                               std::string & gapped_alt,
                               long pos,
                               long & ref_to_seq_offset
  )
{
  Variant new_var_c;
  new_var_c.abs_pos = 0;
  auto a_reference_begin = gapped_ref.begin();
  auto a_sequence_begin = gapped_alt.begin();
  auto a_reference_end = gapped_ref.end();
  auto a_sequence_end = gapped_alt.end();

  /*
  if (*a_sequence_begin != '-')
  {
    // BOOST_LOG_TRIVIAL(debug) << "I did not start on a gap, perhaps increasing the EXTRA_BASES can help.";
    return new_var_c;
  }
   */

  // Remove clipping prefix
  while ((*a_sequence_begin == '-' || *a_sequence_begin != *a_reference_begin) &&
         a_sequence_begin != a_sequence_end
    )
  {
    if (*a_reference_begin != '-')
      ++pos;

    ++a_sequence_begin;
    ++a_reference_begin;
  }

  // We require to find at least READ_END_MATCHES matches before starting new variant discovery
  std::size_t const READ_END_MATCHES = 2;

  {
    std::size_t start_matches = 0;

    while (start_matches < READ_END_MATCHES)
    {
      if (a_sequence_begin == a_sequence_end)
        return new_var_c;

      if (*a_sequence_begin == *a_reference_begin)
        ++start_matches;
      else
        start_matches = 0;

      if (*a_reference_begin != '-')
        ++pos;

      ++a_sequence_begin;
      ++a_reference_begin;
    }
  }

  // Descripes how to go from the reference position to the read position (r)
  ref_to_seq_offset = pos - 2;

  // ..and then remove common prefix
  while (*a_sequence_begin == *a_reference_begin)
  {
    if (a_sequence_begin == a_sequence_end)
      return new_var_c;

    if (*a_reference_begin != '-')
      ++pos;

    ++a_sequence_begin;
    ++a_reference_begin;
  }

  // Then move one back
  --a_reference_begin;
  --pos;

  // Remove clipping suffix
  while (*a_sequence_end == '-')
  {
    if (a_sequence_begin == a_sequence_end)
      return new_var_c;

    --a_sequence_end;
    --a_reference_end;
  }

  // Remove common suffix with at least 7 matches
  {
    std::size_t end_matches = 0;

    while (end_matches < READ_END_MATCHES)
    {
      if (a_sequence_begin == a_sequence_end)
        return new_var_c;

      if (*a_sequence_end == *a_reference_end)
        ++end_matches;
      else
        end_matches = 0;

      --a_sequence_end;
      --a_reference_end;
    }
  }

  ++a_sequence_end;
  ++a_reference_end;

  // Get first base
  char const first_base = *a_reference_begin;
  assert(first_base != '-');

  std::vector<char> reference(1, first_base);
  ++a_reference_begin;

  while (a_reference_begin != a_reference_end)
  {
    reference.push_back(*a_reference_begin);
    ++a_reference_begin;
  }

  if (reference.size() == 1)
    return new_var_c;

  std::vector<char> alt(1, first_base);

  while (a_sequence_begin != a_sequence_end)
  {
    alt.push_back(*a_sequence_begin);
    ++a_sequence_begin;
  }

  // Check if there are any variants
  if (reference == alt)
    return new_var_c;

  new_var_c.abs_pos = pos;
  new_var_c.seqs.push_back(std::move(reference));
  new_var_c.seqs.push_back(std::move(alt));
  return new_var_c;
}


std::vector<VariantCandidate>
find_variants_in_alignment(uint32_t const pos,
                           std::vector<char> const & ref,
                           seqan::Dna5String const & read,
                           std::vector<char> const & qual
                           )
{
  std::vector<VariantCandidate> new_var_candidates;

  assert(ref.size() > 1);
  assert(seqan::length(read) > 1);

  std::string gapped_ref;
  std::string gapped_alt;

  if(!get_capped_strings(gapped_ref, gapped_alt, ref, read))
    return new_var_candidates;

  long ref_to_seq_offset = 0;
  Variant new_var = make_variant_of_gapped_strings(gapped_ref, gapped_alt, pos, ref_to_seq_offset);

  if (new_var.abs_pos == 0)
    return new_var_candidates;

  std::vector<Variant> new_vars =
    extract_sequences_from_aligned_variant(std::move(new_var), SPLIT_VAR_THRESHOLD);

  new_var_candidates.resize(new_vars.size());

  for (unsigned i = 0; i < new_vars.size(); ++i)
  {
    assert(i < new_var_candidates.size());
    auto & new_var = new_vars[i];
    auto & new_var_candidate = new_var_candidates[i];
    assert(new_var.seqs.size() >= 2);
    new_var.normalize();

    // Check if high or low quality
    int64_t const r = new_var.abs_pos - ref_to_seq_offset;
    int64_t const r_end = r + new_var.seqs[1].size();
    ref_to_seq_offset += new_var.seqs[0].size() - new_var.seqs[1].size();

    // In case no qual is available
    if (r_end < static_cast<int>(qual.size()))
    {
      int const MAX_QUAL = *std::max_element(qual.begin() + r, qual.begin() + r_end) - 33;
      new_var_candidate.is_low_qual = MAX_QUAL < 25;
    }

    new_var_candidate.abs_pos = new_var.abs_pos;
    new_var_candidate.seqs = std::move(new_var.seqs);
  }

  return new_var_candidates;
}


void
extract_to_vcf(std::string const & graph_path,
               std::string const & haps_path,
               std::string const & output,
               std::string const & region
  )
{
  std::vector<std::vector<Genotype> > gts;
  load_graph(graph_path);
  std::vector<std::vector<uint32_t> > hap_calls = read_haplotype_calls_from_file(haps_path, gts);

  for (auto it = hap_calls.begin(); it != hap_calls.end(); ++it)
  {
    assert(it->size() > 0);
    assert((*it)[0] == 0);
    std::sort(it->begin(), it->end());

    // Count haplotype call occurences and remove
    std::map<uint32_t, std::size_t> occurences;

    for (auto const & call : *it)
    {
      auto find_it = occurences.find(call);

      if (find_it == occurences.end())
        occurences[call] = 1u;
      else
        ++find_it->second;
    }

    unsigned min_calls{1u};

    while (true)
    {
      std::size_t unique_calls{0u};

      for (auto map_it = occurences.cbegin(); map_it != occurences.cend(); ++map_it)
      {
        if (map_it->second >= min_calls)
          ++unique_calls;
      }

      if (unique_calls <= Options::instance()->max_extracted_haplotypes)
        break;

      BOOST_LOG_TRIVIAL(info) << "[graphtyper::haplotype_extractor] Found a haplotype with "
                              << unique_calls
                              << " different calls. I will require more occurences for each one"
                              << " of them.";
      ++min_calls;
    }

    it->erase(std::unique(it->begin(), it->end()), it->end());
    it->erase(std::remove_if(it->begin(), it->end(), [min_calls, &occurences](uint32_t const val){
        return val > 0 && occurences.at(val) < min_calls;
      }), it->end()
      );
  }

  Vcf hap_extract_vcf(WRITE_BGZF_MODE, output);
  hap_extract_vcf.add_haplotypes_for_extraction(gts, hap_calls);
  hap_extract_vcf.post_process_variants();
  hap_extract_vcf.write(region);
}


} // namespace gyper
