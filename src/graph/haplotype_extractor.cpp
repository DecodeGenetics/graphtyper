#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <unordered_map>

#include <paw/align.hpp>

#include <boost/log/trivial.hpp>

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/graph.hpp>
#include <graphtyper/graph/graph_serialization.hpp>
#include <graphtyper/graph/haplotype_calls.hpp>
#include <graphtyper/graph/haplotype_extractor.hpp>
#include <graphtyper/typer/vcf_writer.hpp>
#include <graphtyper/typer/vcf.hpp>
#include <graphtyper/typer/variant.hpp>
#include <graphtyper/typer/variant_candidate.hpp>
#include <graphtyper/utilities/sequence_operations.hpp> // Remove prefix and suffix
#include <graphtyper/utilities/options.hpp>
#include <graphtyper/utilities/graph_help_functions.hpp>
#include <graphtyper/utilities/type_conversions.hpp>


namespace
{

std::vector<gyper::HaplotypeCall>
read_haplotype_calls(std::vector<std::string> const & haps_paths)
{
  using namespace gyper;

  if (haps_paths.size() == 0)
  {
    BOOST_LOG_TRIVIAL(info) << "[graphtyper::haplotype_extractor] No haplotype calls were extracted";
    std::exit(0);
  }

  std::vector<gyper::HaplotypeCall> hap_calls{};

  // Read first line
  {
    std::string file_path = haps_paths[0] + ".hap";
    hap_calls = load_calls(file_path.c_str());
  }

  // Read remaining lines
  for (long i = 1; i < static_cast<long>(haps_paths.size()); ++i)
  {
    std::string file_path = haps_paths[i] + ".hap";
    std::vector<gyper::HaplotypeCall> new_hap_calls = load_calls(file_path.c_str());
    assert(hap_calls.size() == new_hap_calls.size());

    for (long h = 0; h < static_cast<long>(hap_calls.size()); ++h)
      hap_calls[h].merge_with(new_hap_calls[h]);
  }

  return hap_calls;
}


std::vector<gyper::HaplotypeCall>
read_haplotype_calls_from_file(std::string const & haps_path)
{
  std::vector<std::string> haps_paths;

  // File with a list of haplotype files
  std::ifstream haps_file(haps_path.c_str());

  // Make sure the file could be opened
  if (!haps_file.is_open())
  {
    BOOST_LOG_TRIVIAL(error) << "[graphtyper::haplotype_extractor] Unable to open file '" << haps_path << "'.";
    std::exit(1);
  }

  std::string line;

  while (std::getline(haps_file, line))
  {
    // remove '.hap'
    haps_paths.push_back(line.substr(0, line.size() - 4));
  }

  haps_file.close();
  return read_haplotype_calls(haps_paths);
}


} // anon namespace


namespace gyper
{

bool
get_gapped_strings(std::pair<std::string, std::string> & gapped_strings,
                   std::vector<char> const & ref,
                   std::vector<char> const & seq
                   )
{
  paw::AlignmentOptions<uint8_t> opts;
  opts.set_match(2).set_mismatch(4);
  opts.set_gap_open(6).set_gap_extend(1);
  opts.left_column_free = true;
  opts.right_column_free = true;
  paw::global_alignment(seq, ref, opts);

  paw::AlignmentResults<uint8_t> const & ar = *opts.get_alignment_results();

  if (ar.score == 2 * static_cast<int>(seq.size()))
    return false; // Perfect match
  else if (ar.score < 42)
    return false;

  gapped_strings = ar.get_aligned_strings(seq, ref);
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

  // Descripes how to go from the reference position to the read position (r)
  ref_to_seq_offset = pos;

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

  // Remove common prefix
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
  --a_reference_end;
  --a_sequence_end;

  if (a_sequence_begin == a_sequence_end)
    return new_var_c;

  while (*a_sequence_end == '-')
  {
    if (a_sequence_begin == a_sequence_end)
      return new_var_c;

    --a_sequence_end;
    --a_reference_end;
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
                           std::vector<char> const & seq,
                           std::vector<char> const & qual
                           )
{
  std::vector<VariantCandidate> new_var_candidates;

  assert(ref.size() > 1);
  assert(seq.size() > 1);

  std::pair<std::string, std::string> gapped_strings;

  if (!get_gapped_strings(gapped_strings, ref, seq))
    return new_var_candidates;

  //if (!get_gapped_strings(gapped_ref, gapped_alt, ref, seq))
  //  return new_var_candidates;

  //if (gapped_strings.first != gapped_alt)
  //  std::cerr << gapped_strings.first << " != " << gapped_alt << "\n";
  //else if (gapped_strings.second != gapped_ref)
  //  std::cerr << gapped_strings.second << " != " << gapped_ref << "\n";
  //else
  //{
  //  std::cerr << gapped_strings.first << " and\n" << gapped_strings.second << "\n";
  //}

  std::string & gapped_ref = gapped_strings.second;
  std::string & gapped_alt = gapped_strings.first;

  long ref_to_seq_offset = 0;
  Variant var = make_variant_of_gapped_strings(gapped_ref, gapped_alt, pos, ref_to_seq_offset);

  if (var.abs_pos == 0)
    return new_var_candidates;

  std::vector<Variant> new_vars =
    extract_sequences_from_aligned_variant(std::move(var), SPLIT_VAR_THRESHOLD);

  if (new_vars.size() == 0)
    return new_var_candidates;

  new_var_candidates.reserve(new_vars.size());

  for (long i = 0; i < static_cast<long>(new_vars.size()); ++i)
  {
    //assert(i < new_var_candidates.size());
    auto & new_var = new_vars[i];

#ifndef NDEBUG
    if (new_var.seqs.size() != 2 || new_var.seqs[0].size() == 0 || new_var.seqs[1].size() == 0)
    {
      BOOST_LOG_TRIVIAL(error) << "Error in discovery. " << new_var.seqs.size() << ","
                               << new_var.seqs[0].size() << ","
                               << new_var.seqs[1].size();
      std::exit(1);
    }
#endif // NDEBUG

    VariantCandidate new_var_candidate;
    assert(new_var.seqs.size() == 2);
    assert(new_var.seqs[0].size() > 0);
    assert(new_var.seqs[1].size() > 0);
    assert(new_var.is_normalized());

    // its a sign of a problem if the new variant candidate is not normalized, most likely this means the 50bp where not enough
    if (!new_var.is_normalized())
    {
      BOOST_LOG_TRIVIAL(debug) << "Removed variant candidate since it was not in normalized form " << new_var.print();
      continue;
    }

    // Check if high or low quality
    long const r = std::max(0l, new_var.abs_pos - ref_to_seq_offset - 50l);
    long const r_end = r + new_var.seqs[1].size();
    ref_to_seq_offset += new_var.seqs[0].size() - new_var.seqs[1].size();

    // In case no qual is available
    if (r < static_cast<long>(qual.size()))
    {
      std::vector<char>::const_iterator max_qual_it;

      if (r_end < static_cast<long>(qual.size()))
        max_qual_it = std::max_element(qual.begin() + r, qual.begin() + r_end);
      else
        max_qual_it = std::max_element(qual.begin() + r, qual.end());

      long MAX_QUAL = (*max_qual_it) - 33l;
      new_var_candidate.flags |= static_cast<uint16_t>(static_cast<bool>(MAX_QUAL <= 25l)) << IS_LOW_BASE_QUAL_SHIFT;
    }

    new_var_candidate.abs_pos = new_var.abs_pos;
    new_var_candidate.seqs = std::move(new_var.seqs);
    new_var_candidates.push_back(std::move(new_var_candidate));
  }

  return new_var_candidates;
}


void
post_process_hap_calls(std::vector<gyper::HaplotypeCall> & hap_calls)
{
  Options const & copts = *(Options::const_instance());

  ///*
  for (auto & hap_call : hap_calls)
  {
    auto & calls = hap_call.calls;
    assert(calls.size() > 0);
    assert(calls[0] == 0);
    std::sort(calls.begin(), calls.end());
    std::vector<long> occurences;
    occurences.resize(calls[calls.size() - 1] + 1);

    {
      long prev_call = 0; //hap_call.calls[0];
      long n = 1;

      for (long i = 1; i < static_cast<long>(calls.size()); ++i)
      {
        auto const call = calls[i];

        if (call == prev_call)
        {
          ++n;
          continue;
        }

        assert(prev_call < static_cast<long>(occurences.size()));
        occurences[prev_call] = n;
        prev_call = call;
        n = 1;
      }

      assert(prev_call < static_cast<long>(occurences.size()));
      occurences[prev_call] = n;
    }

    // Make the calls unique
    calls.erase(std::unique(calls.begin(), calls.end()), calls.end());

    if (static_cast<long>(occurences.size()) > copts.max_extracted_haplotypes)
    {
      long min_calls{2};

      while (true)
      {
        long unique_calls =
          std::count_if(occurences.cbegin(), occurences.cend(), [min_calls](long const a){
            return a >= min_calls;
          });

        if (unique_calls <= copts.max_extracted_haplotypes)
          break;


        ++min_calls;
      }

      BOOST_LOG_TRIVIAL(debug) << "Found a variant with a high number of different haplotype calls. "
                               << "I will require " << min_calls << " sample calls.";

      hap_call.calls.erase(std::remove_if(hap_call.calls.begin(),
                                          hap_call.calls.end(),
                                          [min_calls, &occurences](uint16_t const val){
          return val > 0 && occurences[val] < min_calls;
        }), hap_call.calls.end()
                           );
    }
  }
  //*/

  //*
  for (auto & hap_call : hap_calls)
  {
    assert(hap_call.calls.size() >= 1);

    if (hap_call.calls.size() == 1)
      continue;

    std::unordered_set<uint16_t> bad_calls;

    // Check for high impurity if we have multiple samples
    /*
    if (hap_call.num_samples >= 50 && copts.impurity_threshold < 0.25)
    {
      for (long i = 0; i < static_cast<long>(hap_call.haplotype_impurity.size()); ++i)
      {
        double const impurity = static_cast<double>(hap_call.haplotype_impurity[i]) /
                                static_cast<double>(hap_call.num_samples);
        assert(impurity > -0.001);
        assert(impurity < 0.251);

        if (impurity > copts.impurity_threshold)
        {
          BOOST_LOG_TRIVIAL(debug) << "Too high impurity of " << impurity << " on haplotype " << (i + 1);
          bad_calls.insert(i + 1);
        }
      }
    }*/

    auto const & gts = hap_call.gts;
    assert(gts.size() > 0);
    assert(gts.size() == hap_call.read_strand.size());
    uint32_t num{1};

    for (auto gt_it = gts.cbegin(); gt_it != gts.cend(); ++gt_it)
      num *= gt_it->num;

    uint32_t const cnum{num};

    for (uint16_t haplotype_call : hap_call.calls)
    {
      if (haplotype_call == 0)
        continue;

      uint32_t q = cnum;
      uint32_t call = haplotype_call;

      for (long i = 0; i < static_cast<long>(gts.size()); ++i)
      {
        auto const & rs = hap_call.read_strand[i];
        auto const & gt = gts[i];

        q /= gt.num;

        unsigned const div = call / q;
        assert(div < rs.size());

        if (copts.filter_on_read_bias && copts.filter_on_strand_bias && div > 0)
        {
          auto const & rs_num = rs[div];
          long const weight = rs_num.get_weight();
          long const max_bias = rs_num.get_max_bias();

          if (max_bias > (9 * weight / 10) ||
              (weight > 250 && max_bias > (4 * weight / 5)))
          {
            bad_calls.insert(haplotype_call);
#ifndef NDEBUG
            BOOST_LOG_TRIVIAL(debug) << "Heavy read/strand bias: " << max_bias << " / " << weight
                                     << " " << rs_num.str();
#endif // ifndef NDEBUG
            break;
          }
        }

        call %= q;
      }
    }

    if (bad_calls.size() > 0)
    {
      BOOST_LOG_TRIVIAL(debug) << "Number of haplotypes with bad read or strand bias were " << bad_calls.size();

      auto & ca = hap_call.calls;
      ca.erase(std::remove_if(ca.begin(), ca.end(), [&bad_calls](uint16_t val){
          return bad_calls.count(val) >= 1;
        }), ca.end());
    }
  }//*/
}


void
extract_to_vcf(std::string const & graph_path,
               std::string const & haps_path,
               std::string const & output_vcf,
               std::string const & region,
               bool const is_splitting_vars)
{
  load_graph(graph_path);
  std::vector<gyper::HaplotypeCall> hap_calls = read_haplotype_calls_from_file(haps_path);
  post_process_hap_calls(hap_calls);
  Vcf hap_extract_vcf(WRITE_BGZF_MODE, output_vcf);
  hap_extract_vcf.add_haplotypes_for_extraction(hap_calls, is_splitting_vars);

  for (Variant & var : hap_extract_vcf.variants)
    var.normalize();

  hap_extract_vcf.write(region);
}


void
extract_to_vcf(Vcf & hap_extract_vcf,
               std::vector<std::string> const & haps_paths,
               std::string const & output_vcf,
               bool const is_splitting_vars)
{
  std::vector<gyper::HaplotypeCall> hap_calls = read_haplotype_calls(haps_paths);
  post_process_hap_calls(hap_calls);
  hap_extract_vcf.open(WRITE_BGZF_MODE, output_vcf);
  hap_extract_vcf.add_haplotypes_for_extraction(hap_calls, is_splitting_vars);

  for (Variant & var : hap_extract_vcf.variants)
    var.normalize();
}


} // namespace gyper
