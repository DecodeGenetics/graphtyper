#include <cassert> // assert
#include <sstream> // std::ostringstream
#include <string> // std::string
#include <vector> // std::vector
#include <utility>

#include <boost/log/trivial.hpp>

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/absolute_position.hpp>
#include <graphtyper/graph/graph.hpp>
#include <graphtyper/graph/graph_serialization.hpp>
#include <graphtyper/graph/constructor.hpp>
#include <graphtyper/graph/var_record.hpp>
#include <graphtyper/utilities/options.hpp>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/vcf_io.h>


namespace
{

void
append_sv_tag_to_node(std::vector<char> & alt, char const end_with)
{
  std::ostringstream ss;
  ss << "<SV:" << std::setw(7) << std::setfill('0') << gyper::graph.SVs.size();
  std::string sv_id = ss.str();
  std::move(sv_id.c_str(), sv_id.c_str() + sv_id.size(), std::back_inserter(alt));
  alt.push_back(end_with);
}


void
open_tabix(seqan::Tabix & tabix_file, std::string const & tabix_filename, gyper::GenomicRegion const & genomic_region)
{
  seqan::open(tabix_file, tabix_filename.c_str());
  seqan::String<char> header;
  seqan::getHeader(header, tabix_file);

  std::string region_str = genomic_region.to_string();
  seqan::setRegion(tabix_file, region_str.c_str());
}


void
open_reference_genome(seqan::FaiIndex & fasta_index, std::string const & fasta_filename)
{
  if (!seqan::open(fasta_index, fasta_filename.c_str()))
  {
    if (!seqan::build(fasta_index, fasta_filename.c_str()))
    {
      BOOST_LOG_TRIVIAL(error) << "[graphtyper::constructor] Index could not be loaded or built.";
    }
    else if (!seqan::save(fasta_index))
    {
      BOOST_LOG_TRIVIAL(error) << "[graphtyper::constructor] Index could not be saved to disk.";
    }
  }
}


char
complement(char const c)
{
  switch(c)
  {
  case 'A': return 'T';
  case 'C': return 'G';
  case 'G': return 'C';
  case 'T': return 'A';
  default: return c;
  }
}

template<typename TIt>
int
longest_poly_repeat(TIt begin, TIt end, char const c, int & max_i)
{
  int max_repeats = 0;
  int repeats = 0;
  int const SIZE = static_cast<int>(std::distance(begin, end));

  for (int i = 0; begin != end; ++i, ++begin)
  {
    if (*begin == c)
    {
      ++repeats;
      continue;
    }
    else if (repeats > max_repeats)
    {
      max_repeats = repeats;
      max_i = i - 1;
    }

    repeats = 0;
  }

  if (repeats > max_repeats)
  {
    max_i = SIZE - 1;
    max_repeats = repeats;
  }

  return max_repeats;
}


template<typename Tseq>
void
combine_similar_sequences(std::vector<Tseq> & seqs)
{
  // Find the hamming distance between all pairs of ALU sequences
  using Tint = int8_t;
  std::size_t constexpr MAX_HAMMING_DIST_1 = 1;
  std::size_t constexpr MAX_HAMMING_DIST_2 = 2;
  std::vector<std::vector<Tint> > hamming_distances_1(seqs.size(), std::vector<Tint>(seqs.size()));
  std::vector<std::vector<Tint> > hamming_distances_2(seqs.size(), std::vector<Tint>(seqs.size()));

  for (std::size_t i = 1; i < seqs.size(); ++i)
  {
    for (std::size_t j = 0; j < i; ++j)
    {
      std::size_t hamming_dist = 0;

      for (std::size_t k = 0; k < seqs[0].size(); ++k)
      {
        hamming_dist += seqs[i][k] != seqs[j][k];

        if (hamming_dist > MAX_HAMMING_DIST_2)
          break;
      }

      // 1
      hamming_distances_1[i][j] = static_cast<Tint>(hamming_dist > MAX_HAMMING_DIST_1);
      hamming_distances_1[j][i] = hamming_distances_1[i][j];

      // 2
      hamming_distances_2[i][j] = static_cast<Tint>(hamming_dist > MAX_HAMMING_DIST_2);
      hamming_distances_2[j][i] = hamming_distances_2[i][j];
    }
  }

  std::set<std::size_t> selected_seq_indexes;

  {
    Tint minimum_sum_1 = static_cast<Tint>(hamming_distances_1.size());
    Tint minimum_sum_2 = static_cast<Tint>(hamming_distances_1.size());
    std::size_t minimum_sum_index = 0;
    std::set<std::size_t> covered_seq_indexes;

    while (covered_seq_indexes.size() < seqs.size())
    {
      for (std::size_t i = 0; i < hamming_distances_1.size(); ++i)
      {
        if (covered_seq_indexes.count(i) == 1)
          continue;

        Tint sum_1 = 0;
        Tint sum_2 = 0;

        for (std::size_t j = 0; j < hamming_distances_1[i].size(); ++j)
        {
          if (covered_seq_indexes.count(j) == 1)
            continue;

          sum_1 += hamming_distances_1[i][j];
          sum_2 += hamming_distances_2[i][j];
        }

        if (sum_2 < minimum_sum_2 || (sum_2 == minimum_sum_2 && sum_1 < minimum_sum_1))
        {
          minimum_sum_1 = sum_1;
          minimum_sum_2 = sum_2;
          minimum_sum_index = i;
        }
      }

      selected_seq_indexes.insert(minimum_sum_index);

      for (std::size_t j = 0; j < hamming_distances_2[minimum_sum_index].size(); ++j)
      {
        if (hamming_distances_2[minimum_sum_index][j] == 0)
        {
          covered_seq_indexes.insert(j);
        }
      }
    }
  }

  std::vector<Tseq> new_seqs;

  for (auto seq_index : selected_seq_indexes)
  {
    assert (seq_index < seqs.size());
    new_seqs.push_back(seqs[seq_index]);
  }

  std::swap(new_seqs, seqs);
}



std::vector<std::vector<char> > inline
read_reference_graph(std::vector<char> const & prefix,
                     gyper::Graph const & sv_graph,
                     std::string const & chr,
                     uint32_t const begin,
                     uint32_t length
                     )
{
  uint32_t const abs_begin = gyper::absolute_pos.get_absolute_position(chr, begin);
  return sv_graph.get_all_sequences_of_length(abs_begin, length, prefix);
}


void inline
read_reference_seq(std::vector<char> & reference_sequence,
                   seqan::FaiIndex const & fasta_index,
                   unsigned const chrom_idx,
                   uint32_t const begin,
                   uint32_t const length
                   )
{
  seqan::Dna5String ref_seq;

  seqan::readRegion(ref_seq,
                    fasta_index,
                    chrom_idx,
                    begin,
                    begin + length
                    );

  std::move(seqan::begin(ref_seq), seqan::end(ref_seq), std::back_inserter(reference_sequence));
}

std::vector<std::vector<char> > inline
read_reference(std::vector<char> const & prefix,
               seqan::FaiIndex const & fasta_index,
               unsigned const chrom_idx,
               std::string const & chr,
               uint32_t const begin,
               uint32_t const length,
               std::unique_ptr<gyper::Graph> const & sv_graph
               )
{
  if (sv_graph)
  {
    // Read from SV graph
    std::cerr << "DEBUG: Reading from reference graph.\n";
    std::vector<std::vector<char> > seqs = read_reference_graph(prefix, *sv_graph, chr, begin, length);
    std::cerr << "DEBUG seqs before " << seqs.size() << std::endl;

    if (seqs.size() < 100)
      combine_similar_sequences(seqs);

    std::cerr << "DEBUG seqs after " << seqs.size() << std::endl;
    return seqs;
  }
  else
  {
    // No SV graph provided, use the reference sequence
    std::cerr << "DEBUG: Reading from reference sequence.\n";
    std::vector<std::vector<char> > reference_sequences(1, prefix);
    read_reference_seq(reference_sequences[0], fasta_index, chrom_idx, begin, length);
    return reference_sequences;
  }
}


void
read_reference_genome(std::vector<char> & reference_sequence,
                      seqan::FaiIndex const & fasta_index,
                      gyper::GenomicRegion const & genomic_region
                      )
{
  unsigned chrom_idx = 0;

  if (!seqan::getIdByName(chrom_idx, fasta_index, genomic_region.chr.c_str()))
  {
    BOOST_LOG_TRIVIAL(warning) << "[graphtyper::constructor] FAI index has no entry for chromosome '"
                               << genomic_region.chr << "' (region = " << genomic_region.to_string() << ")";
    assert(false);
  }
  else
  {
    read_reference_seq(reference_sequence, fasta_index, chrom_idx, genomic_region.begin, genomic_region.end - genomic_region.begin);
  }
}


std::vector<char>
read_reference_genome_ends(seqan::FaiIndex const & fasta_index,
                           unsigned const chrom_idx,
                           uint32_t const begin,
                           uint32_t const end,
                           uint32_t const length
                           )
{
  // Read the ends of the reference sequence
  std::vector<char> ref;

  if ((end - begin) > 2 * length)
  {
    read_reference_seq(ref, fasta_index, chrom_idx, begin, length);
    read_reference_seq(ref, fasta_index, chrom_idx, end - length, length);
  }
  else
  {
    read_reference_seq(ref, fasta_index, chrom_idx, begin, end - begin);
  }

  return ref;
}

} // anon namespace


namespace gyper
{

void
get_alu_sequences(std::vector<std::vector<char> > & seqs_begin,
                  std::vector<std::vector<char> > & seqs_end,
                  SV & sv,
                  VarRecord const & var,
                  seqan::FaiIndex const & fasta_index,
                  unsigned const chrom_idx,
                  uint32_t const EXTRA_SEQUENCE_LENGTH
                  )
{
  assert (sv.or_end != -1);
  assert (sv.or_start != -1);

  if (sv.or_start >= sv.or_end)
    return;

  std::vector<char> alu_prefix; // Sometimes the ALUs have a prefix sequence inserted

  // Original ALU sequence had 140 Ns padded to them... so here is a little trick to fix the positions
  if (sv.or_start < 141)
  {
    // assume that the ALU prefix is the reference before the ALU
    int const n = 141 - sv.or_start;
    read_reference_seq(alu_prefix, fasta_index, chrom_idx, var.pos - n, static_cast<uint32_t>(n));
    sv.or_start = 141;
  }

  if (sv.or_end - sv.or_start < 32)
    return;

  assert(sv.or_end >= 141);
  assert(sv.or_start >= 141);
  sv.or_end -= 141; // Change to 0-based system instread of 141-based system
  sv.or_start -= 141;

  assert(sv.or_end > sv.or_start);
  assert(sv.or_end - sv.or_start >= 32);

  // Read reference sequence before and after to find binding site
  std::vector<char> ref_before;
  std::vector<char> ref_after;
  std::size_t constexpr BINDING_SITE_CHECK_LENGTH = 30;

  read_reference_seq(ref_before, fasta_index, chrom_idx, var.pos - BINDING_SITE_CHECK_LENGTH - 1, BINDING_SITE_CHECK_LENGTH);
  read_reference_seq(ref_after, fasta_index, chrom_idx, var.pos + 1, BINDING_SITE_CHECK_LENGTH);

  int max_i_before = -1;
  int max_i_after = -1;

  // poly(A) are expected to bind before the position, and poly(T) before
  int const binding_length_before = longest_poly_repeat(ref_before.rbegin(),
                                                        ref_before.rend(),
                                                        'A',
                                                        max_i_before);

  int const binding_length_after = longest_poly_repeat(ref_after.begin(),
                                                       ref_after.end(),
                                                       'T',
                                                       max_i_after);

  int constexpr MIN_BINDING_LENGTH = 3;

  if (binding_length_before < MIN_BINDING_LENGTH && binding_length_after < MIN_BINDING_LENGTH)
    return;

  std::vector<std::vector<char> > alus;

  for (auto const & alu : gyper::ALU_SEQUENCES) /*defined in graph/alu_sequences.hpp*/
  {
    if ((sv.or_end - 32) < static_cast<int>(alu.size()))
    {
      std::vector<char> new_alu(alu_prefix);
      std::copy(alu.begin() + sv.or_start, alu.end(), std::back_inserter(new_alu));
      alus.push_back(std::move(new_alu));
    }
  }

  // No ALUs where found with these coordinates
  if (alus.size() == 0)
    return;

  if (binding_length_before >= MIN_BINDING_LENGTH)
  {
    // Add ALUs in forward direction
    std::vector<std::vector<char> > f_alus;

    if (binding_length_after < MIN_BINDING_LENGTH)
      f_alus = std::move(alus);
    else
      f_alus = alus;

    // Add reference betweeen the ALU breakpoint
    for (auto & alu : f_alus)
    {
      if (alu.size() < EXTRA_SEQUENCE_LENGTH)
      {
        read_reference_seq(alu,
                           fasta_index,
                           chrom_idx,
                           var.pos + 1,
                           static_cast<uint32_t>(EXTRA_SEQUENCE_LENGTH - alu.size())
                           );

        assert(alu.size() == EXTRA_SEQUENCE_LENGTH);
      }
      else if (alu.size() > EXTRA_SEQUENCE_LENGTH)
      {
        alu.resize(EXTRA_SEQUENCE_LENGTH);
      }
    }

    // Reduce the amount of sequences to add by combining similar sequences
    combine_similar_sequences(f_alus);

    // Append to results
    std::move(f_alus.begin(), f_alus.end(), std::back_inserter(seqs_begin));

    // Read target site duplicates
    std::vector<char> tsd(EXTRA_SEQUENCE_LENGTH - max_i_before, 'A');
    read_reference_seq(tsd, fasta_index, chrom_idx, var.pos - max_i_before, max_i_before);
    seqs_end.push_back(std::move(tsd));
  }

  if (binding_length_after >= MIN_BINDING_LENGTH)
  {
    for (auto & alu : alus)
    {
      // Complemet the ALU sequence
      std::transform(alu.begin(), alu.end(), alu.begin(), complement);

      if (alu.size() < EXTRA_SEQUENCE_LENGTH)
      {
        std::size_t const n = EXTRA_SEQUENCE_LENGTH - alu.size();
        std::vector<char> ref;
        read_reference_seq(ref, fasta_index, chrom_idx, var.pos - n - 1, n);
        std::transform(ref.begin(), ref.end(), ref.begin(), complement);
        std::move(ref.begin(), ref.end(), std::back_inserter(alu));
      }
      else if (alu.size() > EXTRA_SEQUENCE_LENGTH)
      {
        alu.resize(EXTRA_SEQUENCE_LENGTH);
      }
    }

    combine_similar_sequences(alus);

    std::move(alus.rbegin(), alus.rend(), std::back_inserter(seqs_end));

    // Read target site duplication
    std::vector<char> tsd;
    tsd.reserve(EXTRA_SEQUENCE_LENGTH);
    read_reference_seq(tsd, fasta_index, chrom_idx, var.pos + 1, max_i_after);
    tsd.resize(EXTRA_SEQUENCE_LENGTH, 'T');
    seqs_begin.push_back(std::move(tsd));
  }
}


void
get_duplicated_sequences(std::vector<char> & dup_begin,
                         std::vector<char> & dup_end,
                         SV const & sv,
                         VarRecord const & var,
                         seqan::FaiIndex const & fasta_index,
                         unsigned const chrom_idx,
                         uint32_t const EXTRA_SEQUENCE_LENGTH
                         )
{
  if (sv.or_end == -1)
  {
    if (sv.or_start == -1)
    {
      // Case 1: Both breakpoints are known and the duplication is tandem
      // Duplication starts after: beginPos
      // Duplicated sequence is [beginPos + 1, beginPos + min(SVLEN, 150)]
      BOOST_LOG_TRIVIAL(debug) << "Case 1: Both breakpoints are known and the duplication is tandem.";

      // Read the duplicated sequence
      std::vector<char> dup = read_reference_genome_ends(fasta_index,
                                                         chrom_idx,
                                                         var.pos + 1,
                                                         var.pos + 1 + sv.length,
                                                         EXTRA_SEQUENCE_LENGTH
                                                         );

      // Both breakpoints are in the duplicated sequence if it is greater or equal to the extra_sequence_length
      if (dup.size() >= EXTRA_SEQUENCE_LENGTH)
      {
        // Read the beginning of the duplicated sequence
        std::copy(dup.begin(), dup.begin() + EXTRA_SEQUENCE_LENGTH, std::back_inserter(dup_begin));

        // Read the end of the duplicated sequence
        std::copy(dup.end() - EXTRA_SEQUENCE_LENGTH, dup.end(), std::back_inserter(dup_end));
      }
      else
      {
        std::size_t padding_size = EXTRA_SEQUENCE_LENGTH - dup.size();

        /// Read the beginning of the duplicated sequence
        std::copy(dup.begin(), dup.end(), std::back_inserter(dup_begin));
        read_reference_seq(dup_begin, fasta_index, chrom_idx, var.pos + 1, padding_size);

        /// Read the end of the duplicated sequence
        // Make sure we don't read before the chromosome starts!
        padding_size = std::min(padding_size, static_cast<std::size_t>(var.pos));
        read_reference_seq(dup_end, fasta_index, chrom_idx, var.pos - padding_size, padding_size);
        std::move(dup.begin(), dup.end(), std::back_inserter(dup_end));
      }
    }
    else
    {
      // Case 2: ORSTART given but OREND not. Only one of the breakpoints is known
      // Duplication starts after: beginPos
      // Duplicated sequence is [ORSTART, ORSTART + 150]
      BOOST_LOG_TRIVIAL(debug) << "Case 2: ORSTART="
                               << sv.or_start
                               << " but OREND=N/A. Only one of the breakpoints is known";
      read_reference_seq(dup_begin, fasta_index, chrom_idx, sv.or_start - 1, EXTRA_SEQUENCE_LENGTH);
    }
  }
  else
  {
    if (sv.or_start != -1)
    {
      // Case 3: OREND given but ORSTART not. Only one of the breakpoints is known
      // Duplication starts after: beginPos
      // Duplicated sequence is [OREND - 150, OREND]
      BOOST_LOG_TRIVIAL(debug) << "Case 3: OREND="
                               << sv.or_end
                               << " but ORSTART=N/A. Only one of the breakpoints is known";

      // Do not read before the chromosome starts!
      std::size_t const start_reading_at = std::min(static_cast<std::size_t>(EXTRA_SEQUENCE_LENGTH),
                                                    static_cast<std::size_t>(sv.or_end - 1)
                                                    );
      read_reference_seq(dup_begin, fasta_index, chrom_idx, start_reading_at, EXTRA_SEQUENCE_LENGTH);
    }
    else
    {
      // Case 4: Both ORSTART and OREND are given. Both breakpoints are known but the duplication is not tandem
      // Duplication starts after: beginPos
      // END=min(ORSTART+150,OREND)
      // Duplicated sequence is [ORSTART, END][beginPos + 1, beginPos + 150 - END - ORSTART]
      BOOST_LOG_TRIVIAL(debug) << "Case 4: Both ORSTART and OREND are given. Both breakpoints "
                                  "are known but the duplication is not tandem";
      /// Read the duplicated sequence
      std::vector<char> dup = read_reference_genome_ends(fasta_index,
                                                         chrom_idx,
                                                         sv.or_start - 1,
                                                         sv.or_end - 1,
                                                         EXTRA_SEQUENCE_LENGTH
                                                         );

      // Both breakpoints are in the duplicated sequence if it is greater or equal to the extra_sequence_length
      if (dup.size() >= EXTRA_SEQUENCE_LENGTH)
      {
        // Read the beginning of the duplicated sequence
        std::copy(dup.begin(), dup.begin() + EXTRA_SEQUENCE_LENGTH, std::back_inserter(dup_begin));

        // Read the end of the duplicated sequence
        std::copy(dup.end() - EXTRA_SEQUENCE_LENGTH, dup.end(), std::back_inserter(dup_end));
      }
      else
      {
        std::size_t padding_size = EXTRA_SEQUENCE_LENGTH - dup.size();

        /// Read the beginning of the duplicated sequence
        std::copy(dup.begin(), dup.end(), std::back_inserter(dup_begin));
        read_reference_seq(dup_begin, fasta_index, chrom_idx, var.pos + 1, padding_size);

        /// Read the end of the duplicated sequence
        // Make sure we don't read before the chromosome starts!
        padding_size = std::min(padding_size, static_cast<std::size_t>(var.pos));
        read_reference_seq(dup_end, fasta_index, chrom_idx, var.pos - padding_size, padding_size);
        std::move(dup.begin(), dup.end(), std::back_inserter(dup_end));
      }
    }
  }
}


template<typename Tint>
bool
parse_info(std::string const check_if_key,
           std::string const & key,
           std::string const & val,
           Tint & parsed_val
           )
{
  if (check_if_key == key)
  {
    std::istringstream is{val};
    is >> parsed_val;

    if (!is.eof())
    {
      BOOST_LOG_TRIVIAL(error) << "[graphtyper::constructor] VCF PARSE ERROR: Could not parse "
                               << check_if_key
                               << " from the INFO field (val = '"
                               << val
                               << "')";
      std::exit(1);
    }

    return true;
  }

  return false;
}


void
add_var_record(std::vector<VarRecord> & var_records,
               seqan::VcfRecord const & vcf_record,
               seqan::FaiIndex const & fasta_index,
               GenomicRegion genomic_region,
               std::unique_ptr<Graph> const & //sv_graph
  )
{
  assert(vcf_record.beginPos >= 0);
  VarRecord var(static_cast<uint32_t>(vcf_record.beginPos));
  var.is_sv = static_cast<char>(vcf_record.alt[0]) == '<';

  if (var.is_sv)
  {
    // Read the FASTA index of the chromosome (since it won't change)
    unsigned chrom_idx = 0;

    if (!seqan::getIdByName(chrom_idx, fasta_index, genomic_region.chr.c_str()))
    {
      BOOST_LOG_TRIVIAL(error) << "[graphtyper::constructor] FAI index has no entry for "
                                << "chromosome '"
                               << genomic_region.chr
                               << "' (region = " << genomic_region.to_string() << ")";

      std::exit(30);
    }

    SV sv;

    // Parse the INFO field
    seqan::StringSet<seqan::CharString> infos;
    seqan::strSplit(infos, vcf_record.info, seqan::EqualsChar<';'>());

    for (auto const & info : infos)
    {
      long const EQ_SIGN_POS = std::distance(seqan::begin(info),
                                             std::find(seqan::begin(info), seqan::end(info), '=')
                                             );

      std::string const key = std::string(seqan::begin(info), seqan::begin(info) + EQ_SIGN_POS);
      std::string const val = EQ_SIGN_POS < static_cast<int>(seqan::length(info))
                              ? std::string(seqan::begin(info) + EQ_SIGN_POS + 1, seqan::end(info))
                              : std::string("");



      if (key == "SVTYPE")
      {
        if (val == "DEL")
          sv.type = DEL;
        else if (val == "DEL:ME:ALU")
          sv.type = DEL_ALU;
        else if (val == "DUP")
          sv.type = DUP;
        else if (val == "INV")
          sv.type = INV;
        else if (val == "INS")
          sv.type = INS;
        else if (val == "INS:ME:ALU")
          sv.type = INS_ALU;
        else
          sv.type = OTHER;
      }
      else if (parse_info("END", key, val, sv.end)) {}
      else if (parse_info("SVSIZE", key, val, sv.size)) {}
      else if (parse_info("SVLEN", key, val, sv.length))
      {
        // Make SVLEN positive
        sv.length = std::abs(sv.length);
      }
      else if (parse_info("NCLUSTERS", key, val, sv.n_clusters)) {}
      else if (parse_info("ORSTART", key, val, sv.or_start)) {}
      else if (parse_info("OREND", key, val, sv.or_end)) {}
      else if (key == "SEQ" && val != ".")
      {
        sv.seq.clear();
        sv.seq.reserve(val.size());
        std::copy(val.begin(), val.end(), std::back_inserter(sv.seq));
      }
    }

    uint32_t constexpr EXTRA_SEQUENCE_LENGTH = 152;
    assert (EXTRA_SEQUENCE_LENGTH % 2 == 0);

    // Handle structural variants
    switch(sv.type)
    {
    case DEL:
    case DEL_ALU:
    {
      // Handle deletions
      BOOST_LOG_TRIVIAL(debug) << "A deletion of size "
                               << sv.size
                               << " @ "
                               << (vcf_record.beginPos + 1);

      // Read the first matching reference base
      read_reference_seq(var.ref, fasta_index, chrom_idx, var.pos, 1);

      // Read the alternative allele
      std::vector<char> alt1(var.ref);

      // If there is any inserted sequence, append it
      if (sv.seq.size() > 0 && sv.seq[0] != '.')
        std::copy(sv.seq.begin(), sv.seq.end(), std::back_inserter(alt1));

      // Add breakpoint
      if (alt1.size() < EXTRA_SEQUENCE_LENGTH + 1)
      {
        read_reference_seq(alt1,
                           fasta_index,
                           chrom_idx,
                           var.pos + sv.length + 1,
                           EXTRA_SEQUENCE_LENGTH + 1 - alt1.size()
                           );
      }

      // Append SV tag
      append_sv_tag_to_node(alt1, '<'); // end with
      var.alts.push_back(std::move(alt1));

      // Add the SV
      graph.SVs.push_back(std::move(sv));

      break;
    }


    case DUP:
    {
      BOOST_LOG_TRIVIAL(debug) << "A duplication of size "
                               << sv.size
                               << " @ "
                               << (vcf_record.beginPos + 1);

      // Read the first matching reference base
      read_reference_seq(var.ref, fasta_index, chrom_idx, var.pos, 1);

      std::vector<char> dup_begin; // Beginning of the duplication
      std::vector<char> dup_end; // End of the duplication

      if (sv.or_start == -1 && sv.or_end == -1)
      {
        // Move the duplication to beginPos + sv.length
        assert(vcf_record.beginPos >= 0);
        var.pos = static_cast<uint32_t>(vcf_record.beginPos) + sv.length;

        // Re-read the reference base with the new position
        var.ref.clear();
        read_reference_seq(var.ref, fasta_index, chrom_idx, var.pos, 1);
      }

      get_duplicated_sequences(dup_begin,
                               dup_end,
                               sv,
                               var,
                               fasta_index,
                               chrom_idx,
                               EXTRA_SEQUENCE_LENGTH
        );

      std::vector<char> alt;
      alt.push_back(var.ref.front()); // Add reference base to alternative allele
      std::move(dup_begin.begin(), dup_begin.end(), std::back_inserter(alt));

      if (dup_end.size() == 0)
      {
        append_sv_tag_to_node(alt, '<' /*end with*/);
      }
      else
      {
        append_sv_tag_to_node(alt, '>' /*end with*/);
        std::move(dup_end.begin(), dup_end.end(), std::back_inserter(alt));
        alt.push_back(var.ref.front());
      }

      var.alts.push_back(std::move(alt)); // Add the alternative path
      graph.SVs.push_back(std::move(sv));
      break;
    }

    case INV:
    {
      BOOST_LOG_TRIVIAL(debug) << "An inversion of size "
                               << sv.size
                               << " @ "
                               << (vcf_record.beginPos + 1);

      // Read the first matching reference base
      read_reference_seq(var.ref, fasta_index, chrom_idx, var.pos, 1);

      std::vector<char> inv_begin; // Beginning of the inversion
      std::vector<char> inv_end; // End of the inversion

      // Get the duplicated sequences, now we still need to reverse complement them
      uint32_t INV_EXTRA_LENGTH = 400;
      get_duplicated_sequences(inv_begin,
                               inv_end,
                               sv,
                               var,
                               fasta_index,
                               chrom_idx,
                               INV_EXTRA_LENGTH
        );

      BOOST_LOG_TRIVIAL(debug) << "inv_begin: "
                               << std::string(inv_begin.begin(), inv_begin.end());

      BOOST_LOG_TRIVIAL(debug) << "inv_end: "
                               << std::string(inv_end.begin(), inv_end.end());

      // Find the complement sequence
      std::transform(inv_begin.begin(), inv_begin.end(), inv_begin.begin(), complement);
      std::transform(inv_end.begin(), inv_end.end(), inv_end.begin(), complement);

      // Reverse the sequence
      BOOST_LOG_TRIVIAL(debug) << "inv_begin reverse comp: "
                               << std::string(inv_begin.rbegin(), inv_begin.rend());

      BOOST_LOG_TRIVIAL(debug) << "inv_end reverse comp: "
                               << std::string(inv_end.rbegin(), inv_end.rend());

      std::vector<char> alt;
      alt.push_back(var.ref.front()); // Add reference base to alternative allele
      std::move(inv_end.rbegin(), inv_end.rend(), std::back_inserter(alt));

      if (inv_begin.size() == 0)
      {
        append_sv_tag_to_node(alt, '<' /*end with*/);
      }
      else
      {
        append_sv_tag_to_node(alt, '>' /*end with*/);
        std::move(inv_begin.rbegin(), inv_begin.rend(), std::back_inserter(alt));
        alt.push_back(var.ref.front());
      }

      //alt.push_back(var.ref.front()); // Add reference base
      var.alts.push_back(std::move(alt)); // Add the alternative path
      graph.SVs.push_back(std::move(sv));

      break;
    }

    case INS:
    {
      BOOST_LOG_TRIVIAL(debug) << "An insertion of size "
                               << sv.size
                               << " @ " << (vcf_record.beginPos + 1);

      if (!sv.seq.empty())
      {
        BOOST_LOG_TRIVIAL(debug) << "Sequence given: "
                                 << std::string(sv.seq.begin(), sv.seq.end());

        // Read the first matching reference base
        read_reference_seq(var.ref, fasta_index, chrom_idx, var.pos, 1);
        std::vector<char> alt(var.ref); // Matching padding base

        // If the given SV sequence is more or equal in length as EXTRA_SEQUENCE_LENGTH we extract
        // only sequence at breakpoints.
        if (sv.seq.size() >= EXTRA_SEQUENCE_LENGTH)
        {
          std::copy(sv.seq.begin(),
                    sv.seq.begin() + EXTRA_SEQUENCE_LENGTH,
                    std::back_inserter(alt)
            );

          append_sv_tag_to_node(alt, '>' /*end with*/);

          std::copy(sv.seq.end() - EXTRA_SEQUENCE_LENGTH,
                    sv.seq.end(),
                    std::back_inserter(alt)
            );
        }
        else
        {
          // Padding is the length of  padding sequence before and after the SV breakpoint
          long padding_length = EXTRA_SEQUENCE_LENGTH - sv.seq.size();
          assert(padding_length > 0);
          std::copy(sv.seq.begin(), sv.seq.end(), std::back_inserter(alt));

          read_reference_seq(alt,
                             fasta_index,
                             chrom_idx,
                             var.pos + 1,
                             padding_length
            );

          append_sv_tag_to_node(alt, '>' /*end with*/);

          read_reference_seq(alt,
                             fasta_index,
                             chrom_idx,
                             var.pos - padding_length,
                             padding_length + 1
            );

          std::copy(sv.seq.begin(), sv.seq.end(), std::back_inserter(alt));
        }

        var.alts.push_back(std::move(alt));  // Add the alternative path
      }
      else
      {
        assert (sv.or_start != -1);
        assert (sv.or_end != -1);

        // Read the first matching reference base
        read_reference_seq(var.ref, fasta_index, chrom_idx, var.pos, 1);

        std::vector<char> ins_begin;  // Beginning of the insertion
        std::vector<char> ins_end;  // End of the insertion

        get_duplicated_sequences(ins_begin,
                                 ins_end,
                                 sv,
                                 var,
                                 fasta_index,
                                 chrom_idx,
                                 EXTRA_SEQUENCE_LENGTH
          );

        std::vector<char> alt;
        alt.push_back(var.ref.front());  // Add reference base to alternative allele
        std::move(ins_begin.begin(), ins_begin.end(), std::back_inserter(alt));
        append_sv_tag_to_node(alt, '>' /*end with*/);
        std::move(ins_end.begin(), ins_end.end(), std::back_inserter(alt));
        alt.push_back(var.ref.front());
        var.alts.push_back(std::move(alt));  // Add the alternative path
      }

      graph.SVs.push_back(std::move(sv));
      break;
    }

    case INS_ALU:
    {
      if (sv.or_start == -1 || sv.or_end == -1)
        break;

      BOOST_LOG_TRIVIAL(debug) << "An ALU insertion of size "
                               << sv.size
                               << " @ "
                               << (vcf_record.beginPos + 1);

      // Read the first matching reference base
      read_reference_seq(var.ref, fasta_index, chrom_idx, var.pos, 1);

      std::vector<std::vector<char> > ins_begin; // Beginning of the insertion
      std::vector<std::vector<char> > ins_end; // End of the insertion

      get_alu_sequences(ins_begin, ins_end, sv, var, fasta_index, chrom_idx, EXTRA_SEQUENCE_LENGTH);

      if (ins_begin.size() == 0 || ins_end.size() == 0)
        break;

      // Make sure these two sequence are the same size by inserting empty sequences
      if (ins_begin.size() > ins_end.size())
        ins_end.resize(ins_begin.size(), std::vector<char>(0));
      else
        ins_begin.resize(ins_end.size(), std::vector<char>(0));

      assert(ins_begin.size() == ins_end.size());

      for (std::size_t i = 0; i < ins_begin.size(); ++i)
      {
        auto && seq_begin = ins_begin[i];
        auto && seq_end = ins_end[i];

        std::vector<char> alt;
        alt.push_back(var.ref.front()); // Add reference base to alternative allele
        std::move(seq_begin.begin(), seq_begin.end(), std::back_inserter(alt));
        append_sv_tag_to_node(alt, '>' /*end with*/);
        std::move(seq_end.begin(), seq_end.end(), std::back_inserter(alt));
        alt.push_back(var.ref.front());
        var.alts.push_back(std::move(alt)); // Add the alternative path
      }

      // Add the variant
      graph.SVs.push_back(std::move(sv));

      break;
    }


    default:{return;} // Skip adding this variant
    }

    /*
    if (var.alts.size() > 0)
    {
      // Inefficient because it runs even when there is no logging, comment in release
      std::ostringstream ss;
      ss << var.pos << " ";
      ss << std::string(var.ref.begin(), var.ref.end());

      for (auto const & alt : var.alts)
        ss  << " " << std::string(alt.cbegin(), alt.cend());

      BOOST_LOG_TRIVIAL(debug) << ss.str();
    }
    */
  }
  else
  {
    // Handle typical small variants (SNPs, indels <50 bps)
    var.ref = std::vector<char>(seqan::begin(vcf_record.ref), seqan::end(vcf_record.ref));
    seqan::StringSet<seqan::CharString> seqan_alts;
    seqan::strSplit(seqan_alts, vcf_record.alt, seqan::EqualsChar<','>());
    var.alts.reserve(seqan::length(seqan_alts));

    for (unsigned i = 0; i < seqan::length(seqan_alts); ++i)
      var.alts.push_back(std::vector<char>(seqan::begin(seqan_alts[i]), seqan::end(seqan_alts[i])));

    // Remove alleles with Ns
    var.alts.erase(std::remove_if(var.alts.begin(),
                              var.alts.end(),
                              [](std::vector<char> const & alt)
                              {
                                return std::find(alt.begin(), alt.end(), 'N') != alt.end() || alt.size() == 0;
                              }
                     ), var.alts.end());


  }

  // Only add var if there are some alternative alleles
  if (var.alts.size() > 0)
  {
    //if (var.is_sv)
    //{
    //  BOOST_LOG_TRIVIAL(debug) << "Adding a SVvariant at position " << var.pos;
    //}

    var_records.push_back(std::move(var));
  }
}


void
construct_graph(std::string const & reference_filename,
                std::string const & vcf_filename,
                std::string const & region,
                bool const use_absolute_positions
                )
{
  graph = Graph(use_absolute_positions);

  if (region.substr(0, 3) != "chr" && region.substr(0, 1) != ".")
    graph.use_prefix_chr = false;

  BOOST_LOG_TRIVIAL(info) << "[graphtyper::constructor] Constructing graph for region " << region;
  GenomicRegion genomic_region(region);

  BOOST_LOG_TRIVIAL(info) << "[graphtyper::constructor] Reading FASTA file located at " << reference_filename;

  std::vector<char> reference_sequence;
  seqan::FaiIndex fasta_index;
  std::vector<VarRecord> var_records;

  // Load the reference genome
  open_reference_genome(fasta_index, reference_filename);
  read_reference_genome(reference_sequence, fasta_index, genomic_region);

  std::unique_ptr<Graph> sv_graph;

  // Load a optional SV graph, in case SVs are being called
  if (Options::instance()->sv_graph.size() > 0)
  {
    sv_graph = std::unique_ptr<Graph>(new Graph(load_secondary_graph(Options::instance()->sv_graph)));
    BOOST_LOG_TRIVIAL(debug) << "DEBUG: sv_graph is of size " << sv_graph->size();
  }

  if (vcf_filename.size() > 0)
  {
    BOOST_LOG_TRIVIAL(info) << "[graphtyper::constructor] INFO: Reading VCF file located at " << vcf_filename;

    // Load region using the tabix index
    seqan::Tabix tabix_file;
    open_tabix(tabix_file, vcf_filename, genomic_region);

    // Load VCF record in loaded region
    seqan::VcfRecord vcf_record;

    // Read records
    bool is_read_record = seqan::readRegion(vcf_record, tabix_file);

    while (is_read_record)
    {
      if (vcf_record.beginPos >= static_cast<int64_t>(genomic_region.begin) &&
          vcf_record.beginPos < static_cast<int64_t>(genomic_region.end))
      {
        add_var_record(var_records, vcf_record, fasta_index, genomic_region, sv_graph);
      }

      is_read_record = seqan::readRegion(vcf_record, tabix_file);
    }

    // Remove duplicate alternative alleles
    for (auto & var_record : var_records)
    {
      std::sort(var_record.alts.begin(), var_record.alts.end());
      var_record.alts.erase(std::unique(var_record.alts.begin(), var_record.alts.end()),
                            var_record.alts.end()
      );
    }

    genomic_region.check_if_var_records_match_reference_genome(var_records, reference_sequence);

    for (auto & var_record : var_records)
      genomic_region.add_reference_to_record_if_they_have_a_matching_prefix(var_record, reference_sequence);

#ifndef NDEBUG
    genomic_region.check_if_var_records_match_reference_genome(var_records, reference_sequence);
#endif
  }

  seqan::clear(fasta_index); // Close reference genome FASTA

  // Sort var_records by position in increasing order
  std::sort(var_records.begin(), var_records.end(), [](VarRecord const & a, VarRecord const & b){return a.pos < b.pos;});

  graph.add_genomic_region(std::move(reference_sequence), std::move(var_records), std::move(genomic_region));

  if (!graph.check())
  {
    BOOST_LOG_TRIVIAL(error) << "[graphtyper::graph] Problem creating graph.";
    std::exit(1);
  }

  BOOST_LOG_TRIVIAL(info) << "[graphtyper::constructor] Graph was successfully constructed.";

  // Create all specials positions
  graph.create_special_positions();
}


} // namespace gyper
