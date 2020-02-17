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
#include <graphtyper/utilities/gzstream.hpp>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/vcf_io.h>


namespace
{


template <typename Tint>
bool
parse_info_int(std::string const check_if_key,
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


bool
parse_info_str(std::string const check_if_key,
               std::string const & key,
               std::string const & val,
               std::vector<char> & parsed_val
               )
{
  if (check_if_key == key)
  {
    if (val.size() == 0 || val[0] == '.')
      return true;

    parsed_val.clear();
    parsed_val.reserve(val.size());
    std::copy(val.begin(), val.end(), std::back_inserter(parsed_val));
    return true;
  }

  return false;
}


bool
parse_info_sv_type(std::string const check_if_key,
                   std::string const & key,
                   std::string const & val,
                   gyper::SVTYPE & parsed_val
                   )
{
  using namespace gyper;

  if (key == check_if_key)
  {
    if (val == "DEL")
      parsed_val = DEL;
    else if (val == "DEL:ME:ALU")
      parsed_val = DEL_ALU;
    else if (val == "DUP")
      parsed_val = DUP;
    else if (val == "INV")
      parsed_val = INV;
    else if (val == "INS")
      parsed_val = INS;
    else if (val == "INS:ME:ALU")
      parsed_val = INS_ALU;
    else if (val == "BND")
      parsed_val = BND;
    else
      parsed_val = OTHER;

    return true;
  }

  return false;
}


bool
parse_info_inv_type(std::string const & key, gyper::INVTYPE & parsed_val)
{
  using namespace gyper;

  if (key == "INV3")
  {
    parsed_val = INV3;
    return true;
  }
  else if (key == "INV5")
  {
    parsed_val = INV5;
    return true;
  }

  return false;
}


unsigned
get_chrom_idx(seqan::FaiIndex const & fasta_index, std::string const & chrom)
{
  // Read the FASTA index of the chromosome (since it won't change)
  unsigned chrom_idx = 0;

  if (!seqan::getIdByName(chrom_idx, fasta_index, chrom.c_str()))
  {
    std::cerr << "[graphtyper::constructor] ERROR: FAI index has no entry for "
              << "contig/chromosome '" << chrom << "'" << std::endl;

    std::exit(30);
  }

  return chrom_idx;
}


void
append_sv_tag_to_node(std::vector<char> & alt)
{
  std::ostringstream ss;
  ss << "<SV:" << std::setw(7) << std::setfill('0') << gyper::graph.SVs.size() << ">";
  std::string sv_id = ss.str();
  std::move(sv_id.begin(), sv_id.end(), std::back_inserter(alt));
}


void
open_tabix(seqan::Tabix & tabix_file,
           std::string const & tabix_filename,
           gyper::GenomicRegion const & genomic_region
           )
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
  // Read contigs and add them to the graph
  {
    std::ifstream f(fasta_filename + std::string(".fai"));

    if (!f.is_open() || f.eof())
    {
      BOOST_LOG_TRIVIAL(error) << "[graphtyper::constructor] Failed to open FASTA index '"
                               << fasta_filename << ".fai'";
      std::exit(111);
    }

    for (std::string line; std::getline(f, line);)
    {
      std::istringstream ss{line};
      gyper::Contig new_contig;
      ss >> new_contig.name;
      ss >> new_contig.length;
      gyper::graph.contigs.push_back(std::move(new_contig));
    }
  }

  // Check if the total length of contigs is too large
  {
    uint64_t sum = 0;

    for (gyper::Contig const & contig : gyper::graph.contigs)
      sum += contig.length;

    if (sum > 0x00000000FFFFFFFFull)
    {
      BOOST_LOG_TRIVIAL(error) << "[graphtyper::construct] Graphtyper only supports references "
                               << "with a total length less than 4,294,967,295 (2^32 - 1)\n";
      std::exit(112);
    }
  }

  if (!seqan::open(fasta_index, fasta_filename.c_str()))
  {
    if (!seqan::build(fasta_index, fasta_filename.c_str()))
    {
      std::cerr << "[graphtyper::constructor] ERROR: Index could not be loaded or built.\n";
      std::exit(31);
    }
    else if (!seqan::save(fasta_index))
    {
      std::cerr << "[graphtyper::constructor] ERROR: Index could not be saved to disk.\n";
      std::exit(32);
    }
  }
}


char
complement(char const c)
{
  switch (c)
  {
  case 'A': return 'T';

  case 'C': return 'G';

  case 'G': return 'C';

  case 'T': return 'A';

  default: return c;
  }
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


void
read_reference_genome(std::vector<char> & reference_sequence,
                      seqan::FaiIndex const & fasta_index,
                      gyper::GenomicRegion const & genomic_region
                      )
{
  unsigned chrom_idx = get_chrom_idx(fasta_index, genomic_region.chr);
  read_reference_seq(reference_sequence,
                     fasta_index,
                     chrom_idx,
                     genomic_region.begin,
                     genomic_region.end - genomic_region.begin
                     );
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
add_sv_breakend(SV & sv,
                VarRecord & var,
                seqan::VcfRecord const & vcf_record,
                seqan::FaiIndex const & fasta_index,
                unsigned const chrom_idx,
                uint32_t const EXTRA_SEQUENCE_LENGTH
                )
{
  // Read the first matching reference base
  read_reference_seq(var.ref, fasta_index, chrom_idx, var.pos, 1);

  assert(var.ref.size() > 0);
  std::vector<char> alt(seqan::begin(vcf_record.alt), seqan::end(vcf_record.alt));
  sv.original_alt = alt;

  auto invalid_bnd_alt = [&var, &alt]()
                         {
                           std::cerr << "[graphtyper::constructor] ERROR: Invalid breakend alt allele: "
                                     << std::string(begin(alt), end(alt))
                                     << " at position " << var.pos << std::endl;
                           std::exit(2);
                         };

  auto parse_chromosome_name = [&](std::vector<char> const & seq, char const c) -> std::string
                               {
                                 auto find_it = std::find(seq.begin(), seq.end(), c);
                                 auto find_colon_it = std::find(find_it + 1, seq.end(), ':');

                                 if (find_colon_it == seq.end())
                                   invalid_bnd_alt();

                                 return std::string(find_it + 1, find_colon_it);
                               };

  auto parse_position = [&](std::vector<char> const & seq, char const c) -> long
                        {
                          auto find_colon_it = std::find(seq.begin() + 1, seq.end(), ':');
                          auto find_end_it = std::find(find_colon_it + 1, seq.end(), c);

                          if (find_end_it == seq.end())
                            invalid_bnd_alt();

                          long pos = 0;
                          std::istringstream ss{std::string(find_colon_it + 1, find_end_it)};
                          ss >> pos;

                          if (!ss.eof())
                            invalid_bnd_alt();

                          return pos;
                        };

  assert(std::count(begin(alt), end(alt), ',') == 0); // Not multi-allelic
  auto find_it = std::find(alt.begin(), alt.end(), '[');
  std::vector<char> bnd;

  if (find_it != alt.end())
  {
    std::string const chrom = parse_chromosome_name(alt, '[');
    auto const chrom_idx_2 = get_chrom_idx(fasta_index, chrom);
    long pos = parse_position(alt, '[');

    if (find_it != alt.begin())
    {
      // Case 1: S SNNN[chr:pos[ => Extending sequence right of chr:pos
      BOOST_LOG_TRIVIAL(debug) << "[graphtyper::constructor] BND variant case 1 @ " << var.pos + 1;
      bnd = std::vector<char>(var.ref);
      std::copy(alt.begin() + 1, find_it, std::back_inserter(bnd)); // Extra insertion

      // Length to extract from the mate locus
      auto const len = EXTRA_SEQUENCE_LENGTH - bnd.size() + 1;
      read_reference_seq(bnd, fasta_index, chrom_idx_2, pos, len); // Read mate locus
      append_sv_tag_to_node(bnd); // Put SV tag
    }
    else
    {
      // Case 2: S [chr:pos[NNNS => Extending reversed sequence left of chr:pos
      BOOST_LOG_TRIVIAL(debug) << "[graphtyper::constructor] BND variant case 2 @ " << var.pos + 1;
      append_sv_tag_to_node(bnd);

      auto find2_it = std::find(find_it + 1, alt.end(), '[');

      if (find2_it == alt.end())
        invalid_bnd_alt();

      auto const len = EXTRA_SEQUENCE_LENGTH - std::distance(find2_it, alt.end()) - 1;
      std::vector<char> seq;
      read_reference_seq(seq, fasta_index, chrom_idx_2, pos, len);

      // Complement the sequence
      std::transform(seq.begin(), seq.end(), seq.begin(), complement);

      // Copy the sequence in reverse
      std::copy(seq.rbegin(), seq.rend(), std::back_inserter(bnd));
      std::copy(find2_it + 1, alt.end(), std::back_inserter(bnd));
    }
  }
  else
  {
    find_it = std::find(alt.begin(), alt.end(), ']');

    if (find_it == alt.end())
      invalid_bnd_alt();

    std::string const chrom = parse_chromosome_name(alt, ']');
    auto const chrom_idx_2 = get_chrom_idx(fasta_index, chrom);
    long pos = parse_position(alt, ']');

    if (find_it == alt.begin())
    {
      // Case 3: S ]chr:pos]NNS => Take sequence from chr:pos and extend it to the left of S
      BOOST_LOG_TRIVIAL(debug) << "[graphtyper::constructor] BND variant case 3 @ " << var.pos + 1;
      append_sv_tag_to_node(bnd);
      auto find2_it = std::find(find_it + 1, alt.end(), ']');

      if (find2_it == alt.end())
        invalid_bnd_alt();

      auto const len = EXTRA_SEQUENCE_LENGTH - std::distance(find2_it, alt.end()) - 1;
      read_reference_seq(bnd, fasta_index, chrom_idx_2, pos - len - 1, len);
      std::copy(find2_it + 1, alt.end(), std::back_inserter(bnd));
    }
    else
    {
      bnd = std::vector<char>(var.ref);
      // Case 4: S SNN]chr:pos] => Take sequence from chr:pos, reverse complement it and extend it
      // to the right of S
      BOOST_LOG_TRIVIAL(debug) << "[graphtyper::constructor] BND variant case 4 @ " << var.pos + 1;
      std::copy(alt.begin() + 1, find_it, std::back_inserter(bnd)); // Extra insertion
      auto const len = EXTRA_SEQUENCE_LENGTH - bnd.size() + 1;
      std::vector<char> seq;
      read_reference_seq(seq, fasta_index, chrom_idx_2, pos - len - 1, len);

      // Copy the reverse complement
      std::transform(seq.begin(), seq.end(), seq.begin(), complement);
      std::copy(seq.rbegin(), seq.rend(), std::back_inserter(bnd));
      append_sv_tag_to_node(bnd);
    }
  }

  // Append SV tag
  var.alts.push_back(std::move(bnd));

  // Add the SV
  graph.SVs.push_back(std::move(sv));
}


void
add_sv_deletion(SV & sv,
                VarRecord & var,
                seqan::FaiIndex const & fasta_index,
                unsigned const chrom_idx,
                uint32_t const EXTRA_SEQUENCE_LENGTH
                )
{
  // Read the first matching reference base
  read_reference_seq(var.ref, fasta_index, chrom_idx, var.pos, 1);

  // Read the alternative allele
  std::vector<char> alt1(var.ref);

  // If there is any inserted sequence, append it
  if (sv.seq.size() > 0 && sv.seq[0] != '.')
    std::copy(sv.seq.begin(), sv.seq.end(), std::back_inserter(alt1));
  else if (sv.ins_seq.size() > 0 && sv.ins_seq[0] != '.')
    std::copy(sv.ins_seq.begin(), sv.ins_seq.end(), std::back_inserter(alt1));

  // Add breakpoint
  if (alt1.size() < EXTRA_SEQUENCE_LENGTH + 1)
  {
    read_reference_seq(alt1,
                       fasta_index,
                       chrom_idx,
                       var.pos + sv.seq.size() + sv.size + 1,
                       static_cast<unsigned>(EXTRA_SEQUENCE_LENGTH + 1 - alt1.size())
                       );
  }

  // Append SV tag
  append_sv_tag_to_node(alt1);
  var.alts.push_back(std::move(alt1));

  // Add the SV
  sv.model = "BREAKPOINT";
  graph.SVs.push_back(std::move(sv));
}


void
add_sv_insertion(SV & sv,
                 VarRecord & var,
                 seqan::VcfRecord const & vcf_record,
                 seqan::FaiIndex const & fasta_index,
                 unsigned const chrom_idx,
                 uint32_t const EXTRA_SEQUENCE_LENGTH
                 )
{
  assert(seqan::length(vcf_record.ref) >= 1);

  if (static_cast<char>(vcf_record.ref[0]) != 'N')
    var.ref = std::vector<char>(begin(vcf_record.ref), end(vcf_record.ref));
  else
    read_reference_seq(var.ref, fasta_index, chrom_idx, var.pos, 1);

  if (!sv.seq.empty())
  {
    BOOST_LOG_TRIVIAL(debug) << "[graphtyper::construtor] Insertion sequence given: "
                             << std::string(sv.seq.begin(), sv.seq.end());

    // Read the first matching reference base
    std::vector<char> alt1;

    // Matching padding base, 1st breakpoint
    read_reference_seq(alt1, fasta_index, chrom_idx, var.pos, 1);
    std::vector<char> alt2(alt1); // Matching padding base, 2nd breakpoint

    // If the given SV sequence is more or equal in length as EXTRA_SEQUENCE_LENGTH we extract
    // only sequence at breakpoints
    if (sv.seq.size() >= EXTRA_SEQUENCE_LENGTH)
    {
      std::copy(sv.seq.begin(),
                sv.seq.begin() + EXTRA_SEQUENCE_LENGTH,
                std::back_inserter(alt1)
                );

      append_sv_tag_to_node(alt1);
      sv.related_sv = static_cast<int>(graph.SVs.size()) + 1; // Next SV is related
      sv.model = "BREAKPOINT1";
      graph.SVs.push_back(sv);
      append_sv_tag_to_node(alt2);

      std::copy(sv.seq.end() - EXTRA_SEQUENCE_LENGTH,
                sv.seq.end(),
                std::back_inserter(alt2)
                );

      sv.related_sv = static_cast<int>(graph.SVs.size()) - 1;  // Previous SV is related
      sv.model = "BREAKPOINT2";
      graph.SVs.push_back(std::move(sv));
    }
    else
    {
      // Padding is the length of padding sequence before and after the SV breakpoint
      long padding_length = EXTRA_SEQUENCE_LENGTH - sv.seq.size();
      assert(padding_length > 0);
      std::copy(sv.seq.begin(), sv.seq.end(), std::back_inserter(alt1));

      read_reference_seq(alt1,
                         fasta_index,
                         chrom_idx,
                         var.pos + 1,
                         padding_length
                         );

      append_sv_tag_to_node(alt1);
      sv.related_sv = static_cast<int>(graph.SVs.size()) + 1; // Next SV is related
      sv.model = "BREAKPOINT1";
      graph.SVs.push_back(sv);
      append_sv_tag_to_node(alt2);

      read_reference_seq(alt2,
                         fasta_index,
                         chrom_idx,
                         var.pos - padding_length,
                         padding_length + 1
                         );

      std::copy(sv.seq.begin(), sv.seq.end(), std::back_inserter(alt2));
      sv.related_sv = static_cast<int>(graph.SVs.size()) - 1;  // Previous SV is related
      sv.model = "BREAKPOINT2";
      graph.SVs.push_back(std::move(sv));
    }

    var.alts.push_back(std::move(alt1));  // Add the alternative path 1
    var.alts.push_back(std::move(alt2));  // Add the alternative path 2
  }
  else if (sv.or_start != -1 && sv.or_end != -1)
  {
    /// No insertion sequence given but or_start and or_end define a sequence
    // Read the first matching reference base
    std::vector<char> alt1;

    // Matching padding base, 1st breakpoint
    read_reference_seq(alt1, fasta_index, chrom_idx, var.pos, 1);
    std::vector<char> alt2;

    std::vector<char> ins = read_reference_genome_ends(fasta_index,
                                                       chrom_idx,
                                                       sv.or_start - 1,
                                                       sv.or_end,
                                                       EXTRA_SEQUENCE_LENGTH
                                                       );

    /// Both breakpoints are entirely in the sequence if it is greater or equal to the
    /// EXTRA_SEQUENCE_LENGTH
    if (ins.size() >= EXTRA_SEQUENCE_LENGTH)
    {
      // Read the beginning of the duplicated sequence
      std::copy(ins.begin(), ins.begin() + EXTRA_SEQUENCE_LENGTH, std::back_inserter(alt1));

      append_sv_tag_to_node(alt1);
      sv.related_sv = static_cast<int>(graph.SVs.size()) + 1; // Next SV
      sv.model = "BREAKPOINT1";
      graph.SVs.push_back(sv);
      append_sv_tag_to_node(alt2);

      // Read the end of the duplicated sequence
      std::copy(ins.end() - EXTRA_SEQUENCE_LENGTH, ins.end(), std::back_inserter(alt2));

      sv.related_sv = static_cast<int>(graph.SVs.size()) - 1; // Previous SV
      sv.model = "BREAKPOINT2";
      graph.SVs.push_back(std::move(sv));
    }
    else
    {
      std::size_t padding_size = EXTRA_SEQUENCE_LENGTH - ins.size();

      /// Read the beginning of the duplicated sequence
      std::copy(ins.begin(), ins.end(), std::back_inserter(alt1));
      read_reference_seq(alt1, fasta_index, chrom_idx, var.pos + 1, padding_size);
      append_sv_tag_to_node(alt1);
      sv.related_sv = static_cast<int>(graph.SVs.size()) + 1; // Next SV
      sv.model = "BREAKPOINT1";
      graph.SVs.push_back(sv);

      /// Read the end of the duplicated sequence
      append_sv_tag_to_node(alt2);
      // Make sure we don't read before the chromosome starts!
      padding_size = std::min(padding_size, static_cast<std::size_t>(var.pos));
      read_reference_seq(alt2, fasta_index, chrom_idx, var.pos - padding_size, padding_size);
      std::move(ins.begin(), ins.end(), std::back_inserter(alt2));

      sv.related_sv = static_cast<int>(graph.SVs.size()) - 1; // Previous SV
      sv.model = "BREAKPOINT2";
      graph.SVs.push_back(std::move(sv));
    }

    var.alts.push_back(std::move(alt1)); // Add the alternative path 1
    var.alts.push_back(std::move(alt2)); // Add the alternative path 2
  }
  else if (sv.ins_seq_left.size() > 0 || sv.ins_seq_right.size() > 0)
  {
    BOOST_LOG_TRIVIAL(debug) << "[graphtyper::construct] Insertion with incomplete sequence.";
    std::vector<char> left;
    std::vector<char> right;

    if (sv.ins_seq_left.size() > 0)
    {
      auto const len = std::min(static_cast<std::size_t>(EXTRA_SEQUENCE_LENGTH),
                                sv.ins_seq_left.size()
                                );

      std::copy(sv.ins_seq_left.begin(),
                sv.ins_seq_left.begin() + len,
                std::back_inserter(left)
                );
    }

    if (sv.ins_seq_right.size() > 0)
    {
      auto const len = std::min(static_cast<std::size_t>(EXTRA_SEQUENCE_LENGTH),
                                sv.ins_seq_right.size()
                                );

      std::copy(sv.ins_seq_right.begin(),
                sv.ins_seq_right.begin() + len,
                std::back_inserter(right)
                );
    }

    if (left.size() > 1 && right.size() > 0)
    {
      BOOST_LOG_TRIVIAL(debug) << "[graphtyper::construct] Both breakpoints defined.";

      {
        // Breakpoint 1
        std::vector<char> alt1(var.ref);
        std::move(begin(left), end(left), std::back_inserter(alt1));
        append_sv_tag_to_node(alt1);
        sv.model = "BREAKPOINT1";
        sv.related_sv = (int)graph.SVs.size() + 1;
        graph.SVs.push_back(sv);
        var.alts.push_back(std::move(alt1));
      }

      {
        // Breakpoint 2
        std::vector<char> alt2;
        append_sv_tag_to_node(alt2);
        std::move(begin(right), end(right), std::back_inserter(alt2));
        sv.model = "BREAKPOINT2";
        sv.related_sv = (int)graph.SVs.size() - 1;
        graph.SVs.push_back(std::move(sv));
        var.alts.push_back(std::move(alt2));
      }
    }
    else if (left.size() > 1)
    {
      BOOST_LOG_TRIVIAL(debug) << "[graphtyper::construct] Only breakpoint 1 defined.";
      std::vector<char> alt1(var.ref);
      std::move(begin(left), end(left), std::back_inserter(alt1));
      append_sv_tag_to_node(alt1);
      sv.model = "BREAKPOINT1";
      graph.SVs.push_back(std::move(sv));
      var.alts.push_back(std::move(alt1));
    }
    else if (right.size() > 0)
    {
      BOOST_LOG_TRIVIAL(debug) << "[graphtyper::construct] Only breakpoint 2 defined.";
      std::vector<char> alt2;
      append_sv_tag_to_node(alt2);
      std::move(begin(right), end(right), std::back_inserter(alt2));
      sv.model = "BREAKPOINT2";
      graph.SVs.push_back(std::move(sv));
      var.alts.push_back(std::move(alt2));
    }
  }
  else
  {
    BOOST_LOG_TRIVIAL(warning) << "[graphtyper::constructor] I do not know how to add an insertion"
                               << " at position " << var.pos;
  }
}


/// Adds a SV duplication variant that will later be added to the graph
void
add_sv_duplication(std::vector<VarRecord> & var_records,
                   SV & sv,
                   VarRecord & var,
                   seqan::FaiIndex const & fasta_index,
                   unsigned const chrom_idx,
                   uint32_t const EXTRA_SEQUENCE_LENGTH
                   )
{
  // Read the first matching reference base
  read_reference_seq(var.ref, fasta_index, chrom_idx, var.pos, 1);

  if (sv.or_end == -1)
  {
    if (sv.or_start == -1)
    {
      /// Case 1: Both breakpoints are known and the duplication is tandem
      /// Duplication starts after: beginPos
      /// Duplicated sequence is [beginPos + 1, beginPos + min(SVLEN, 150)]
      BOOST_LOG_TRIVIAL(debug) << "Case 1: Both breakpoints are known and the duplication is tandem.";

      // Read the duplicated sequence
      std::vector<char> dup = read_reference_genome_ends(fasta_index,
                                                         chrom_idx,
                                                         var.pos + 1,
                                                         var.pos + sv.length + 1,
                                                         EXTRA_SEQUENCE_LENGTH
                                                         );

      BOOST_LOG_TRIVIAL(debug) << "Duplicated sequence: " << std::string(dup.begin(), dup.end());

      VarRecord var2(var); // First breakpoint is at a different location
      // Change the position of the duplication to the end of it (since there is the breakpoint)
      var.pos += sv.length;

      // Re-read the reference base with the new position
      var.ref.clear();
      read_reference_seq(var.ref, fasta_index, chrom_idx, var.pos, 1);
      std::vector<char> dup_begin(var.ref);
      std::copy(sv.ins_seq.begin(), sv.ins_seq.end(), std::back_inserter(dup_begin));
      std::vector<char> dup_end;

      // Both breakpoints are in the duplicated sequence if it is greater or equal to the
      // extra_sequence_length
      if (dup.size() >= EXTRA_SEQUENCE_LENGTH)
      {
        // Read the beginning of the duplicated sequence
        std::copy(dup.begin(), dup.begin() + EXTRA_SEQUENCE_LENGTH, std::back_inserter(dup_begin));

        append_sv_tag_to_node(dup_begin);
        sv.related_sv = static_cast<int>(graph.SVs.size()) + 1; // Next SV
        sv.model = "BREAKPOINT1";
        graph.SVs.push_back(sv);

        /// Read the end of duplication
        append_sv_tag_to_node(dup_end);

        // Read the end of the duplicated sequence
        std::copy(dup.end() - EXTRA_SEQUENCE_LENGTH, dup.end(), std::back_inserter(dup_end));
        std::copy(sv.ins_seq.begin(), sv.ins_seq.end(), std::back_inserter(dup_end));

        sv.related_sv = static_cast<int>(graph.SVs.size()) - 1;  // Previous SV
        sv.model = "BREAKPOINT2";
        graph.SVs.push_back(std::move(sv));
      }
      else
      {
        /// Both breakpoints are known and the duplication is smaller than EXTRA_SEQUENCE_LENGTH
        std::size_t padding_size = EXTRA_SEQUENCE_LENGTH - dup.size();

        /// Read the beginning of the duplicated sequence
        std::copy(dup.begin(), dup.end(), std::back_inserter(dup_begin));
        read_reference_seq(dup_begin, fasta_index, chrom_idx, var.pos + 1, padding_size);
        append_sv_tag_to_node(dup_begin);
        sv.model = "BREAKPOINT1";
        sv.related_sv = static_cast<int>(graph.SVs.size()) + 1; // Next SV
        graph.SVs.push_back(sv);

        /// Read the end of the duplicated sequence
        // Make sure we don't read before the chromosome starts!
        padding_size = std::min(padding_size, static_cast<std::size_t>(var2.pos));
        append_sv_tag_to_node(dup_end);
        read_reference_seq(dup_end,
                           fasta_index,
                           chrom_idx,
                           var2.pos - padding_size + 1,
                           padding_size
                           );

        std::move(dup.begin(), dup.end(), std::back_inserter(dup_end));
        sv.related_sv = static_cast<int>(graph.SVs.size()) - 1;
        sv.model = "BREAKPOINT2";
        graph.SVs.push_back(std::move(sv));
      }

      var.alts.push_back(std::move(dup_begin));
      var2.alts.push_back(std::move(dup_end));
      var_records.push_back(std::move(var2));
    }
    else
    {
      // Case 2: ORSTART given but OREND not. Only one of the breakpoints is known
      // Duplication starts after: beginPos
      // Duplicated sequence is [ORSTART, ORSTART + 150]
      BOOST_LOG_TRIVIAL(debug) << "Case 2: ORSTART="
                               << sv.or_start
                               << " but OREND=N/A. Only one of the breakpoints is known";
      std::vector<char> dup_begin(var.ref);
      std::copy(sv.ins_seq.begin(), sv.ins_seq.end(), std::back_inserter(dup_begin));
      read_reference_seq(dup_begin, fasta_index, chrom_idx, sv.or_start - 1, EXTRA_SEQUENCE_LENGTH);
      append_sv_tag_to_node(dup_begin);
      sv.model = "BREAKPOINT1";
      var.alts.push_back(std::move(dup_begin));
      graph.SVs.push_back(sv);
    }
  }
  else
  {
    assert(sv.or_start == -1); // We should never see duplication with two non-tandem breakpoints

    // Case 3: OREND given but ORSTART not. Only one of the breakpoints is known
    // Duplication starts after: beginPos
    // Duplicated sequence is [OREND - 150, OREND]
    BOOST_LOG_TRIVIAL(debug) << "Case 3: OREND="
                             << sv.or_end
                             << " but ORSTART=N/A. Only one of the breakpoints is known";

    // Do not read before the chromosome starts!
    std::size_t const start_reading_at = std::max(static_cast<std::size_t>(EXTRA_SEQUENCE_LENGTH),
                                                  static_cast<std::size_t>(sv.or_end)
                                                  );

    std::vector<char> dup_begin;
    append_sv_tag_to_node(dup_begin);
    read_reference_seq(dup_begin,
                       fasta_index,
                       chrom_idx,
                       start_reading_at - EXTRA_SEQUENCE_LENGTH,
                       EXTRA_SEQUENCE_LENGTH);

    std::copy(sv.ins_seq.begin(), sv.ins_seq.end(), std::back_inserter(dup_begin));
    var.alts.push_back(std::move(dup_begin));
    sv.model = "BREAKPOINT2";
    graph.SVs.push_back(std::move(sv));
  }
}


/// Adds a SV inversion variant that will later be added to the graph
void
add_sv_inversion(std::vector<VarRecord> & var_records,
                 SV & sv,
                 VarRecord & var,
                 seqan::FaiIndex const & fasta_index,
                 unsigned const chrom_idx,
                 uint32_t const EXTRA_SEQUENCE_LENGTH
                 )
{
  // Read the first matching reference base
  read_reference_seq(var.ref, fasta_index, chrom_idx, var.pos, 1);

  if (sv.inv_type == INV3)
  {
    sv.or_end = sv.end;
  }
  else if (sv.inv_type == INV5)
  {
    sv.or_start = sv.begin;
    sv.begin += sv.size;
    var.pos += sv.size;
    var.ref.clear();
    read_reference_seq(var.ref, fasta_index, chrom_idx, var.pos, 1);
  }

  if (sv.or_end == -1)
  {
    if (sv.or_start == -1)
    {
      /// Case 1: Both breakpoints are known and the duplication is tandem
      /// Duplication starts after: beginPos
      /// Duplicated sequence is [beginPos + 1, beginPos + min(SVLEN, 150)]
      BOOST_LOG_TRIVIAL(debug) << "Case 1: Both breakpoints are known and the duplication is "
                               << "tandem.";

      // Read the duplicated sequence
      std::vector<char> dup = read_reference_genome_ends(fasta_index,
                                                         chrom_idx,
                                                         var.pos + 1,
                                                         var.pos + sv.length + 1,
                                                         EXTRA_SEQUENCE_LENGTH
                                                         );

      std::vector<char> inv(dup.rbegin(), dup.rend());
      std::transform(inv.begin(), inv.end(), inv.begin(), complement);

      BOOST_LOG_TRIVIAL(debug) << "Inverted sequence: " << std::string(inv.begin(), inv.end());

      std::vector<char> inv_begin(var.ref);
      std::copy(sv.ins_seq.begin(), sv.ins_seq.end(), std::back_inserter(inv_begin));

      VarRecord var2(var); // First breakpoint is at a different location
      var2.pos += sv.length;
      var2.ref.clear();
      read_reference_seq(var2.ref, fasta_index, chrom_idx, var2.pos, 1);
      std::vector<char> inv_end;

      // Both breakpoints are in the inverted sequence if it is greater or equal to the
      // extra_sequence_length
      if (inv.size() >= EXTRA_SEQUENCE_LENGTH)
      {
        // Read the beginning of the inverted sequence
        std::copy(inv.begin(), inv.begin() + EXTRA_SEQUENCE_LENGTH, std::back_inserter(inv_begin));

        append_sv_tag_to_node(inv_begin);
        sv.related_sv = static_cast<int>(graph.SVs.size()) + 1;  // Next SV
        sv.model = "BREAKPOINT1";
        graph.SVs.push_back(sv);

        /// Read the end of inversion
        append_sv_tag_to_node(inv_end);

        // Read the end of the inverted sequence
        std::copy(inv.end() - EXTRA_SEQUENCE_LENGTH, inv.end(), std::back_inserter(inv_end));
        std::copy(sv.ins_seq.begin(), sv.ins_seq.end(), std::back_inserter(inv_end));

        sv.related_sv = static_cast<int>(graph.SVs.size()) - 1;  // Previous SV
        sv.model = "BREAKPOINT2";
        graph.SVs.push_back(std::move(sv));
      }
      else
      {
        /// Both breakpoints are known and the duplication is smaller than EXTRA_SEQUENCE_LENGTH
        std::size_t padding_size = EXTRA_SEQUENCE_LENGTH - inv.size();

        /// Read the beginning of the duplicated sequence
        std::copy(inv.begin(), inv.end(), std::back_inserter(inv_begin));
        read_reference_seq(inv_begin, fasta_index, chrom_idx, var.pos + 1, padding_size);
        append_sv_tag_to_node(inv_begin);
        sv.model = "BREAKPOINT1";
        sv.related_sv = static_cast<int>(graph.SVs.size()) + 1; // Next SV
        graph.SVs.push_back(sv);

        /// Read the end of the duplicated sequence
        // Make sure we don't read before the chromosome starts!
        padding_size = std::min(padding_size, static_cast<std::size_t>(var2.pos));
        append_sv_tag_to_node(inv_end);
        read_reference_seq(inv_end,
                           fasta_index,
                           chrom_idx,
                           var2.pos - padding_size + 1,
                           padding_size
                           );

        std::move(inv.begin(), inv.end(), std::back_inserter(inv_end));
        std::copy(sv.ins_seq.begin(), sv.ins_seq.end(), std::back_inserter(inv_end));
        sv.related_sv = static_cast<int>(graph.SVs.size()) - 1;
        sv.model = "BREAKPOINT2";
        graph.SVs.push_back(std::move(sv));
      }

      var.alts.push_back(std::move(inv_begin));
      var2.alts.push_back(std::move(inv_end));
      var_records.push_back(std::move(var2));
    }
    else
    {
      // Case 2: ORSTART given but OREND not. Only one of the breakpoints is known
      // Duplication starts after: beginPos
      // Duplicated sequence is [ORSTART, ORSTART + 150]
      BOOST_LOG_TRIVIAL(debug) << "Case 2: ORSTART="
                               << sv.or_start
                               << " but OREND=N/A. Only one of the breakpoints is known";
      std::vector<char> dup;
      read_reference_seq(dup, fasta_index, chrom_idx, sv.or_start - 1, EXTRA_SEQUENCE_LENGTH);
      std::transform(dup.begin(), dup.end(), dup.begin(), complement);

      std::vector<char> inv;
      append_sv_tag_to_node(inv);
      std::copy(dup.rbegin(), dup.rend(), std::back_inserter(inv));
      std::copy(sv.ins_seq.begin(), sv.ins_seq.end(), std::back_inserter(inv));
      sv.model = "BREAKPOINT2";
      var.alts.push_back(std::move(inv));
      graph.SVs.push_back(sv);
    }
  }
  else
  {
    assert(sv.or_start == -1); // We should never see inversion with two non-tandem breakpoints

    // Case 3: OREND given but ORSTART not. Only one of the breakpoints is known
    // Duplication starts after: beginPos
    // Duplicated sequence is [OREND - 150, OREND]
    BOOST_LOG_TRIVIAL(debug) << "Case 3: OREND="
                             << sv.or_end
                             << " but ORSTART=N/A. Only one of the breakpoints is known";

    // Do not read before the chromosome starts!
    std::size_t const start_reading_at = std::max(static_cast<std::size_t>(EXTRA_SEQUENCE_LENGTH),
                                                  static_cast<std::size_t>(sv.or_end)
                                                  );

    std::vector<char> dup;
    read_reference_seq(dup,
                       fasta_index,
                       chrom_idx,
                       start_reading_at - EXTRA_SEQUENCE_LENGTH,
                       EXTRA_SEQUENCE_LENGTH);

    std::transform(dup.begin(), dup.end(), dup.begin(), complement);

    std::vector<char> inv(var.ref);
    std::copy(sv.ins_seq.begin(), sv.ins_seq.end(), std::back_inserter(inv));
    std::copy(dup.rbegin(), dup.rend(), std::back_inserter(inv));
    append_sv_tag_to_node(inv);

    var.alts.push_back(std::move(inv));
    sv.model = "BREAKPOINT1";
    graph.SVs.push_back(std::move(sv));
  }
}


std::vector<seqan::VcfRecord>
split_multi_allelic(seqan::VcfRecord && vcf_record)
{
  std::vector<seqan::VcfRecord> vcf_records;
  long const NUM_ALT_ALLELES = 1 + std::count(begin(vcf_record.alt), end(vcf_record.alt), ',');
  vcf_records.reserve(NUM_ALT_ALLELES);

  if (NUM_ALT_ALLELES == 1)
  {
    vcf_records.push_back(std::move(vcf_record));
    return vcf_records;
  }

  // Parse the INFO field
  seqan::StringSet<seqan::CharString> alts;
  seqan::strSplit(alts, vcf_record.alt, seqan::EqualsChar<','>());
  assert(NUM_ALT_ALLELES == (long)length(alts));

  for (long i = 0; i < NUM_ALT_ALLELES; ++i)
  {
    seqan::VcfRecord new_rec(vcf_record);
    new_rec.alt = alts[i];
    vcf_records.push_back(std::move(new_rec));
  }

  return vcf_records;
}


void
transform_sv_records(seqan::VcfRecord & vcf_record,
                     seqan::FaiIndex const & fasta_index,
                     GenomicRegion const & genomic_region
                     )
{
  if (vcf_record.beginPos == 0)
    return;

  // Make sure it is a biallelic variant
  assert(seqan::length(vcf_record.alt) > 0ul);
  assert(std::count(seqan::begin(vcf_record.alt), seqan::end(vcf_record.alt), ',') == 0);

  if (std::find(begin(vcf_record.alt), end(vcf_record.alt), '<') != end(vcf_record.alt) ||
      std::find(begin(vcf_record.alt), end(vcf_record.alt), ']') != end(vcf_record.alt) ||
      std::find(begin(vcf_record.alt), end(vcf_record.alt), '[') != end(vcf_record.alt)
      )
  {
    return; // Already transformed
  }


  // Change record to a SV if the difference between the allele sizes is large enough
  int const size_diff = static_cast<int>(seqan::length(vcf_record.alt)) -
                        static_cast<int>(seqan::length(vcf_record.ref));

  if (size_diff <= -50) // DEL
  {
    BOOST_LOG_TRIVIAL(debug) << "[graphtyper::constructor] Transformed an SV deletion @ "
                             << vcf_record.beginPos + 1
                             << " with size diff " << size_diff;

    std::string seq = "";

    if (static_cast<char>(vcf_record.ref[0]) != static_cast<char>(vcf_record.alt[0]))
    {
      unsigned const chrom_idx = get_chrom_idx(fasta_index, genomic_region.chr);
      --vcf_record.beginPos;

      std::vector<char> ref_before;
      read_reference_seq(ref_before, fasta_index, chrom_idx, vcf_record.beginPos, 1);
      assert(ref_before.size() == 1);

      std::string ref_before_str(ref_before.begin(), ref_before.end());
      vcf_record.ref = ref_before_str.c_str();
      seqan::CharString new_alt(vcf_record.ref);
      seqan::append(new_alt, vcf_record.alt);
      vcf_record.alt = new_alt;
    }

    if (seqan::length(vcf_record.alt) > 1)
      seq = std::string(begin(vcf_record.alt) + 1, end(vcf_record.alt));

    std::ostringstream ss;

    if (seqan::length(vcf_record.info) > 0)
      ss << ";";

    ss << "SVTYPE=DEL;SVLEN=" << (-size_diff)
       << ";SVSIZE=" << (-size_diff)
       << ";END=" << (seq.size() + vcf_record.beginPos + 1 - size_diff);

    if (seq.size() > 0)
      ss << ";SEQ=" << seq;

    vcf_record.ref = vcf_record.ref[0]; // Only first matching base
    vcf_record.alt = "<DEL>";

    seqan::CharString new_info = ss.str();
    seqan::append(vcf_record.info, new_info);
  }
  else if (size_diff >= 50)
  {
    BOOST_LOG_TRIVIAL(debug) << "[graphtyper::constructor] Transformed an SV insertion @ "
                             << vcf_record.beginPos + 1
                             << " with size diff " << size_diff;


    std::string deleted_seq;
    std::string seq;

    // Check if first base does not match
    if (static_cast<char>(vcf_record.ref[0]) != static_cast<char>(vcf_record.alt[0]))
    {
      unsigned const chrom_idx = get_chrom_idx(fasta_index, genomic_region.chr);
      std::vector<char> ref_before;
      --vcf_record.beginPos;

      read_reference_seq(ref_before, fasta_index, chrom_idx, vcf_record.beginPos, 1);
      assert(ref_before.size() == 1);

      {
        seqan::CharString new_ref;
        appendValue(new_ref, ref_before[0]);
        append(new_ref, vcf_record.ref);
        vcf_record.ref = new_ref;
      }

      seq = std::string(begin(vcf_record.alt), end(vcf_record.alt));
    }
    else
    {
      seq = std::string(begin(vcf_record.alt) + 1, end(vcf_record.alt));
    }

    std::ostringstream ss;

    if (seqan::length(vcf_record.info) > 0 && seqan::back(vcf_record.info) != ';')
      ss << ";";

    ss << "SVTYPE=INS;SVLEN=" << size_diff
       << ";SVSIZE=" << size_diff
       << ";SEQ=" << seq;

    vcf_record.alt = "<INS>";
    seqan::CharString new_info = ss.str();
    seqan::append(vcf_record.info, new_info);
  }
}


void
add_var_record(std::vector<VarRecord> & var_records,
               seqan::VcfRecord const & vcf_record,
               seqan::FaiIndex const & fasta_index,
               GenomicRegion genomic_region,
               bool is_sv_graph
               )
{
  assert(vcf_record.beginPos >= 0);
  VarRecord var(static_cast<uint32_t>(vcf_record.beginPos));

  var.is_sv = false;
  auto const & v_alt = vcf_record.alt;

  /// Set var.is_sv by checking if the alt. allele is an SV
  if (seqan::length(v_alt) >= 5)
  {
    for (long i = 0; i < static_cast<long>(seqan::length(v_alt)); ++i)
    {
      if (v_alt[i] == '<' || v_alt[i] == '[' || v_alt[i] == ']')
      {
        var.is_sv = true;

        if (not is_sv_graph)
        {
          BOOST_LOG_TRIVIAL(error) << "[graphtyper::constructor] Found an SV in a non-SV graph at "
                                   << genomic_region.chr << ":" << (vcf_record.beginPos + 1);
          std::exit(1234);
        }

        break;
      }
    }
  }
  /// Done checking if the alt. allele is an SV

  if (var.is_sv)
  {
    // Read the FASTA index of the chromosome (since it won't change)
    unsigned chrom_idx = get_chrom_idx(fasta_index, genomic_region.chr);

    // Replace Ns with correct base from the reference
    if (var.ref.size() == 1 && var.ref[0] == 'N')
    {
      var.ref.clear();
      read_reference_seq(var.ref, fasta_index, chrom_idx, var.pos, 1);
    }

    SV sv;
    sv.begin = var.pos + 1;
    sv.chrom = genomic_region.chr;

    if (seqan::length(vcf_record.id) > 0)
      sv.old_variant_id = seqan::toCString(vcf_record.id);

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

      if (parse_info_sv_type("SVTYPE", key, val, sv.type))
      {}
      else if (parse_info_int("END", key, val, sv.end))
      {}
      else if (parse_info_int("SVSIZE", key, val, sv.size))
      {}
      else if (parse_info_int("SVLEN", key, val, sv.length))
      {}
      else if (parse_info_int("NCLUSTERS", key, val, sv.n_clusters))
      {}
      else if (parse_info_int("ORSTART", key, val, sv.or_start))
      {}
      else if (parse_info_int("OREND", key, val, sv.or_end))
      {}
      else if (parse_info_int("NUM_MERGED_SVS", key, val, sv.num_merged_svs))
      {}
      else if (parse_info_str("SEQ", key, val, sv.seq))
      {}
      else if (parse_info_str("SVINSSEQ", key, val, sv.ins_seq))
      {}
      else if (parse_info_str("LEFT_SVINSSEQ", key, val, sv.ins_seq_left))
      {}
      else if (parse_info_str("RIGHT_SVINSSEQ", key, val, sv.ins_seq_right))
      {}
      else if (parse_info_inv_type(key, sv.inv_type))
      {}
    }

    if (sv.type == NOT_SV)
    {
      BOOST_LOG_TRIVIAL(error) << "[graphtyper::constructor] Improper VCF input. Detected SV with "
                               << "no SVTYPE defined at pos. "
                               << var.pos;
      std::exit(1);
    }

    if (sv.length < 0)
      sv.length = -sv.length; // Make SVLEN positive

    // Make sure sv.length is set (for non-BND)
    if (sv.type != BND && sv.length == 0)
    {
      sv.length = sv.size;

      if (sv.length == 0)
      {
        sv.length = (int)sv.seq.size();

        if (sv.length == 0)
          sv.length = (int)sv.ins_seq.size();
      }

      // If length is less than 50 then it is not an SV
      if (sv.length < 50 && ((int)sv.ins_seq_left.size() + (int)sv.ins_seq_right.size()) < 20)
      {
        BOOST_LOG_TRIVIAL(warning) << "[graphtyper::constructor] WARNING: Ignored SV at " << var.pos
                                   << " because it was less than 50 bp.";
        return;
      }
    }

    // Make sure sv.size is set
    if (sv.size == 0)
      sv.size = sv.length;

    // Make sure sv.end is set
    if (sv.end == 0)
      sv.end = sv.begin + sv.size;

    if (sv.type == INS && sv.seq.size() > 0)
    {
      auto is_similar = [](std::vector<char> const & seq1, std::vector<char> const & seq2) -> bool
                        {
                          seqan::Dna5String seqan_seq1;
                          seqan::Dna5String seqan_seq2;

                          std::size_t constexpr MAX_SIZE = 1000;

                          // At most align MAX_SIZE bp
                          if (seq1.size() > MAX_SIZE && seq2.size() > MAX_SIZE)
                          {
                            seqan_seq1 = std::string(begin(seq1), begin(seq1) + MAX_SIZE).c_str();
                            seqan_seq2 = std::string(begin(seq2), begin(seq2) + MAX_SIZE).c_str();
                          }
                          else
                          {
                            seqan_seq1 = std::string(begin(seq1), end(seq1)).c_str();
                            seqan_seq2 = std::string(begin(seq2), end(seq2)).c_str();
                          }

                          auto const larger_size = std::max(length(seqan_seq1), length(seqan_seq2));
                          seqan::Align<seqan::Dna5String> align;
                          seqan::resize(seqan::rows(align), 2);
                          seqan::assignSource(seqan::row(align, 0), seqan_seq1);
                          seqan::assignSource(seqan::row(align, 1), seqan_seq2);

                          seqan::Score<int, seqan::Simple> scoringScheme(1, -1, -1);
                          seqan::AlignConfig<> alignConfig;

                          int score = seqan::globalAlignment(align, scoringScheme, alignConfig);
                          return static_cast<double>(score) / static_cast<double>(larger_size) >= 0.8;
                        };

      if (static_cast<long>(var.pos - 1l - sv.seq.size()) >= 0l)
      {
        std::vector<char> ref_before;

        read_reference_seq(ref_before,
                           fasta_index,
                           chrom_idx,
                           (unsigned)std::max(0l, static_cast<long>(var.pos - 1l - sv.seq.size())),
                           (unsigned)sv.seq.size()
                           );

        assert(ref_before.size() == sv.seq.size());

        if (is_similar(ref_before, sv.seq))
        {
          // Change to a duplication
          BOOST_LOG_TRIVIAL(debug) << "[graphtyper::constructor] Changed insertion at position "
                                   << var.pos
                                   << " to a duplication and moved it.";
          var.pos -= sv.seq.size();
          sv.type = DUP;
        }
      }

      if (sv.type == INS)
      {
        std::vector<char> ref_after;
        read_reference_seq(ref_after, fasta_index, chrom_idx, var.pos + 1, (uint32_t)sv.seq.size());

        if (is_similar(ref_after, sv.seq))
        {
          // Change to a duplication
          BOOST_LOG_TRIVIAL(debug) << "[graphtyper::constructor] Changed insertion at position "
                                   << var.pos
                                   << " to a duplication.";
          sv.type = DUP;
        }
      }
    }

    uint32_t constexpr EXTRA_SEQUENCE_LENGTH = 152;
    static_assert(EXTRA_SEQUENCE_LENGTH % 2 == 0, "Extra sequence length must be an even number");

    // Handle structural variants
    switch (sv.type)
    {
    case BND:
    {
      BOOST_LOG_TRIVIAL(debug) << "A breakend " << std::string(begin(v_alt), end(v_alt)) << " @ "
                               << (vcf_record.beginPos + 1);
      add_sv_breakend(sv, var, vcf_record, fasta_index, chrom_idx, EXTRA_SEQUENCE_LENGTH);
      break;
    }


    case DEL:
    case DEL_ALU:
    {
      // Handle deletions
      BOOST_LOG_TRIVIAL(debug) << "A deletion of size " << sv.size << " @ "
                               << (vcf_record.beginPos + 1);
      add_sv_deletion(sv, var, fasta_index, chrom_idx, EXTRA_SEQUENCE_LENGTH);
      break;
    }


    case DUP:
    {
      BOOST_LOG_TRIVIAL(debug) << "A duplication of size " << sv.size << " @ "
                               << (vcf_record.beginPos + 1);
      add_sv_duplication(var_records, sv, var, fasta_index, chrom_idx, EXTRA_SEQUENCE_LENGTH);
      break;
    }


    case INS:
    {
      BOOST_LOG_TRIVIAL(debug) << "An insertion of size "
                               << sv.size
                               << " @ " << (vcf_record.beginPos + 1);

      add_sv_insertion(sv, var, vcf_record, fasta_index, chrom_idx, EXTRA_SEQUENCE_LENGTH);
      break;
    }


    case INV:
    {
      BOOST_LOG_TRIVIAL(debug) << "An inversion of size "
                               << sv.size
                               << " @ "
                               << (vcf_record.beginPos + 1);

      add_sv_inversion(var_records, sv, var, fasta_index, chrom_idx, EXTRA_SEQUENCE_LENGTH);
      break;
    }

    default: {
      return;
    }                 // Skip adding this variant
    }
  }
  else
  {
    // Handle typical small variants (SNPs, indels <50 bps)
    var.ref = std::vector<char>(seqan::begin(vcf_record.ref), seqan::end(vcf_record.ref));
    var.alts.push_back(std::vector<char>(seqan::begin(vcf_record.alt), seqan::end(vcf_record.alt)));
    //seqan::StringSet<seqan::CharString> seqan_alts;
    //seqan::strSplit(seqan_alts, vcf_record.alt, seqan::EqualsChar<','>());
    //std::cerr << "t " << seqan::length(seqan_alts) << "\n";
    //var.alts.reserve(seqan::length(seqan_alts));
    //
    //for (unsigned i = 0; i < seqan::length(seqan_alts); ++i)
    //  var.alts.push_back(std::vector<char>(seqan::begin(seqan_alts[i]), seqan::end(seqan_alts[i])));

    // Remove alleles with non-ACGT
    var.alts.erase(
      std::remove_if(
        var.alts.begin(),
        var.alts.end(),
        [](std::vector<char> const & alt)
      {
        if (alt.size() == 0)
          return true;

        for (char const c : alt)
        {
          if (c != 'A' && c != 'C' && c != 'G' && c != 'T')
            return true;
        }

        return false;
      }
        ), var.alts.end());

  }

  // Only add var if there are some alternative alleles
  if (var.alts.size() > 0)
    var_records.push_back(std::move(var));
}


void
construct_graph(std::string const & reference_filename,
                std::string const & vcf_filename,
                std::string const & region,
                bool const is_sv_graph,
                bool const use_absolute_positions,
                bool const check_index)
{
  graph = Graph(use_absolute_positions);
  graph.is_sv_graph = is_sv_graph;

  if (region.substr(0, 3) != "chr" && region.substr(0, 1) != ".")
    graph.use_prefix_chr = false;

  BOOST_LOG_TRIVIAL(debug) << "[graphtyper::constructor] Constructing graph for region "
                           << region;

  GenomicRegion genomic_region(region);

  BOOST_LOG_TRIVIAL(debug) << "[graphtyper::constructor] Reading FASTA file located at "
                           << reference_filename;

  // Load the reference genome
  seqan::FaiIndex fasta_index;
  open_reference_genome(fasta_index, reference_filename);
  absolute_pos.calculate_offsets(graph);

  // Read the reference sequence
  std::vector<char> reference_sequence;
  read_reference_genome(reference_sequence, fasta_index, genomic_region);

  // Read variant records
  std::vector<VarRecord> var_records;

  if (vcf_filename.size() > 0)
  {
    BOOST_LOG_TRIVIAL(debug) << "[graphtyper::constructor] INFO: Reading VCF file located at "
                             << vcf_filename;

    if (check_index)
    {
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
            (vcf_record.beginPos + seqan::length(vcf_record.ref)) <=
            static_cast<int64_t>(genomic_region.end)
            )
        {
          std::vector<seqan::VcfRecord> records = split_multi_allelic(std::move(vcf_record));

          for (auto & rec : records)
          {
            if (is_sv_graph)
              transform_sv_records(rec, fasta_index, genomic_region);

            add_var_record(var_records, rec, fasta_index, genomic_region, is_sv_graph);
          }
        }

        is_read_record = seqan::readRegion(vcf_record, tabix_file);
      }
    }
    else
    {
      // test
      igzstream igz(vcf_filename.c_str()); // input vcf.gz file

      if (!igz.rdbuf()->is_open())
      {
        BOOST_LOG_TRIVIAL(error) << "Could not open VCF " << vcf_filename;
        std::exit(1);
      }

      std::string line;

      while (std::getline(igz, line))
      {
        // Skip header
        if (line.size() > 0 && line[0] == '#')
          continue;

        seqan::VcfRecord vcf_record;
        _insertDataToVcfRecord(vcf_record, line.c_str(), 0);

        if (vcf_record.beginPos >= static_cast<int64_t>(genomic_region.begin) &&
            (vcf_record.beginPos + seqan::length(vcf_record.ref)) <=
            static_cast<int64_t>(genomic_region.end)
            )
        {
          std::vector<seqan::VcfRecord> records = split_multi_allelic(std::move(vcf_record));

          for (auto & rec : records)
          {
            if (is_sv_graph)
              transform_sv_records(rec, fasta_index, genomic_region);

            add_var_record(var_records, rec, fasta_index, genomic_region, is_sv_graph);
          }
        }
      }
    }

    // Remove duplicate alternative alleles
    for (auto & var_record : var_records)
    {
      std::sort(var_record.alts.begin(), var_record.alts.end());
      var_record.alts.erase(std::unique(var_record.alts.begin(), var_record.alts.end()),
                            var_record.alts.end()
                            );
    }

#ifndef NDEBUG
    genomic_region.check_if_var_records_match_reference_genome(var_records, reference_sequence);
#endif // NDEBUG

    for (auto & var_record : var_records)
    {
      genomic_region.add_reference_to_record_if_they_have_a_matching_prefix(var_record,
                                                                            reference_sequence
                                                                            );
    }

#ifndef NDEBUG
    genomic_region.check_if_var_records_match_reference_genome(var_records, reference_sequence);
#endif // NDEBUG
  }

  seqan::clear(fasta_index); // Close reference genome FASTA

  // Sort var_records by position in increasing order
  std::sort(var_records.begin(),
            var_records.end(),
            [](VarRecord const & a, VarRecord const & b){
      return a.pos < b.pos;
    }
            );

  graph.add_genomic_region(std::move(reference_sequence),
                           std::move(var_records),
                           std::move(genomic_region)
                           );

#ifndef NDEBUG
  if (!graph.check())
  {
    BOOST_LOG_TRIVIAL(error) << "[graphtyper::graph] Problem creating graph. Printing graph:";
    gyper::graph.print();
    std::exit(1);
  }
#endif // NDEBUG

  BOOST_LOG_TRIVIAL(debug) << "[graphtyper::constructor] Graph was successfully constructed.";

  // Create all specials positions
  graph.create_special_positions();
}


} // namespace gyper
