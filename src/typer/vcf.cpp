#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>

#include <boost/algorithm/string/predicate.hpp> // boost::algorithm::ends_with
#include <boost/log/trivial.hpp>

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/absolute_position.hpp>
#include <graphtyper/graph/genomic_region.hpp>
#include <graphtyper/graph/graph.hpp>
#include <graphtyper/graph/haplotype_calls.hpp>
#include <graphtyper/graph/reference_depth.hpp>
#include <graphtyper/graph/var_record.hpp>
#include <graphtyper/typer/vcf.hpp>
#include <graphtyper/utilities/graph_help_functions.hpp>
#include <graphtyper/utilities/options.hpp> // gyper::options::instance()
#include <graphtyper/utilities/type_conversions.hpp>

#include <paw/align.hpp>

#include "tbx.h"

namespace
{

std::string
current_date()
{
  time_t  now = time(0);
  struct tm time_structure;
  char buf[16];
  time_structure = *localtime(&now);
  strftime(buf, sizeof(buf), "%Y%m%d", &time_structure);
  std::string const date(buf);
  return date;
}


std::vector<uint8_t>
get_haplotype_phred(gyper::HapSample const & sample, uint32_t const cnum)
{
  using namespace gyper;
  assert(sample.log_score.size() > 0);
  std::vector<uint8_t> hap_phred;

  // Check if all phred scores are zero by finding the first non-zero phred score
  auto find_it = std::find_if(sample.log_score.begin(),
                              sample.log_score.end(),
                              [](uint16_t const val){
      return val != 0;
    });

  if (find_it == sample.log_score.end())
  {
    // If no non-zero phred score is found, the phred scores of the haplotype should be all zero
    hap_phred = std::vector<uint8_t>(cnum * (cnum + 1) / 2, 0u);
  }
  else
  {
    hap_phred = std::vector<uint8_t>(cnum * (cnum + 1) / 2, 255u);

    // First find out what the maximum log score is
    uint16_t const max_log_score = *std::max_element(sample.log_score.begin(),
                                                     sample.log_score.end()
                                                     );

    std::size_t const num_alleles = cnum * (cnum + 1) / 2;
    assert(sample.log_score.size() == num_alleles);

    for (std::size_t i = 0; i < num_alleles; ++i)
    {
      double const LOG10_HALF_times_10 = 3.01029995663981195213738894724493026768189881462108541;

      auto const phred_score =
        std::llround((max_log_score - sample.log_score[i]) * LOG10_HALF_times_10);

      if (phred_score < 255u)
        hap_phred[i] = static_cast<uint8_t>(phred_score);
    }
  }

  return hap_phred;
}


std::vector<std::vector<uint8_t> >
get_genotype_phred(gyper::HapSample const & sample,
                   std::vector<gyper::Genotype> const & gts
                   )
{
  using namespace gyper;

  uint32_t cnum = 1;

  for (auto const & gt : gts)
    cnum *= gt.num;

  std::vector<uint8_t> hap_phred = get_haplotype_phred(sample, cnum);
  std::vector<std::vector<uint8_t> > phred(gts.size());

  auto find_it = std::find_if(hap_phred.begin(), hap_phred.end(), [](uint8_t const val){
      return val != 0;
    });

  if (find_it == hap_phred.end())
  {
    // Every PHRED score is 0
    for (long i = 0; i < static_cast<long>(gts.size()); ++i)
    {
      auto const n = static_cast<long>(gts[i].num) * static_cast<long>(gts[i].num + 1) / 2;
      phred[i] = std::vector<uint8_t>(n, 0u);
    }

    return phred;
  }

  // Some PHRED scores are not 0
  for (long i = 0; i < static_cast<long>(gts.size()); ++i)
  {
    auto const n = static_cast<long>(gts[i].num) * static_cast<long>(gts[i].num + 1) / 2;
    phred[i] = std::vector<uint8_t>(n, 255u);
  }

  for (long i = 0; i < static_cast<long>(hap_phred.size()); ++i)
  {
    // No need to consider this case
    if (hap_phred[i] == 255u)
      continue;

    std::pair<uint32_t, uint32_t> call = to_pair(uint32_t(i));

    uint32_t nnum = cnum;
    uint32_t rem1 = call.first;
    uint32_t rem2 = call.second;

    for (unsigned j = 0; j < gts.size(); ++j)
    {
      nnum /= gts[j].num;
      assert(nnum != 0);
      uint32_t const call1 = rem1 / nnum;
      uint32_t const call2 = rem2 / nnum;
      uint32_t const index = call1 > call2 ?
                             static_cast<uint32_t>(to_index(call2, call1)) :
                             static_cast<uint32_t>(to_index(call1, call2));

      phred[j][index] = std::min(phred[j][index], hap_phred[i]);
      rem1 %= nnum;
      rem2 %= nnum;
    }
  }

  return phred;
}


} // anon namespace


namespace gyper
{

Vcf::Vcf()
{}

Vcf::Vcf(VCF_FILE_MODE const _filemode, std::string const & _filename)
{
  open(_filemode, _filename);
}


void
Vcf::open(VCF_FILE_MODE const _filemode, std::string const & _filename)
{
  filename = _filename;
  set_filemode(_filemode);
}


void
Vcf::set_filemode(VCF_FILE_MODE const _filemode)
{
  if (_filemode == READ_MODE)
  {
    if (boost::algorithm::ends_with(filename, ".vcf.gz"))
      filemode = READ_BGZF_MODE;
    else
      filemode = READ_UNCOMPRESSED_MODE;
  }
  else if (_filemode == WRITE_MODE)
  {
    if (boost::algorithm::ends_with(filename, ".vcf.gz"))
      filemode = WRITE_BGZF_MODE;
    else
      filemode = WRITE_UNCOMPRESSED_MODE;
  }
  else
  {
    filemode = _filemode;
  }
}


/******************
 * CLASS MODIFERS *
 ******************/

// I/O
void
Vcf::open_vcf_file_for_reading()
{
  switch (filemode)
  {
  case READ_UNCOMPRESSED_MODE:
    BOOST_LOG_TRIVIAL(error) << "Cannot read uncompressed VCF, gzip it please.";
    std::exit(1);
    break;

  case READ_BGZF_MODE:
    bgzf_in.open(filename.c_str());

    if (!bgzf_in.rdbuf()->is_open())
    {
      BOOST_LOG_TRIVIAL(error) << "Could not open " << filename;
      std::exit(1);
    }

    break;

  default:
    BOOST_LOG_TRIVIAL(error) << "Trying to read in writing mode.";
    std::exit(1);
  }
}


bool
Vcf::read_record(bool const SITES_ONLY)
{
  std::string line;

  if (!std::getline(bgzf_in, line))
    return false;

  // Get all positions of TABs in the line
  std::vector<std::size_t> const tabs = get_all_pos(line, '\t');
  // We ignore the following fields: qual (5), filter (6), info (7)

  std::string const chrom = get_string_at_tab_index(line, tabs, 0);
  uint32_t const pos = static_cast<uint32_t>(std::stoul(get_string_at_tab_index(line, tabs, 1)));
  std::string const id = get_string_at_tab_index(line, tabs, 2);
  std::string const ref = get_string_at_tab_index(line, tabs, 3);
  std::string const alts = get_string_at_tab_index(line, tabs, 4);
  std::vector<std::size_t> const alt_commas = get_all_pos(alts, ',');

  Variant new_var; // Create a new variant for this position
  new_var.abs_pos = absolute_pos.get_absolute_position(chrom, pos); // Parse positions

  // Check for graphtyper variant ID suffix
  {
    auto start_it = std::find(id.begin(), id.end(), '[');

    if (start_it != id.end())
    {
      auto end_it = std::find(start_it + 1, id.end(), ']');

      if (end_it != id.end())
      {
        new_var.suffix_id = std::string(start_it + 1, end_it);
      }
    }
  }

  // Parse sequences
  new_var.seqs.push_back(gyper::to_vec(std::string(ref)));
  assert(alt_commas.size() >= 2);

  for (int a = 0; a < static_cast<int>(alt_commas.size()) - 1; ++a)
    new_var.seqs.push_back(gyper::to_vec(get_string_at_tab_index(alts, alt_commas, a)));

  // Parse infos
  {
    std::string const info = get_string_at_tab_index(line, tabs, 7);

    // Don't parse anything if the INFO field is empty
    if (info.size() > 0)
    {
      std::vector<std::size_t> const info_semicolons = get_all_pos(info, ';');

      std::unordered_set<std::string> const keys_to_parse(
        {
//          "AC",
          "CR",
          "END",
          "HOMSEQ"
//          "INV3", "INV5",
          "LEFT_SVINSSEQ",
          "NCLUSTERS", "NUM_MERGED_SVS",
          "MQsquared",
          "OLD_VARIANT_ID", "OREND", "ORSTART",
          "PS",
          "RELATED_SV_ID", "RIGHT_SVINSSEQ",
          "SBF", "SBF1", "SBF2",
          "SBR", "SBR1", "SBR2",
          "SVLEN", "SVTYPE", "SVSIZE", "SVMODEL", "SEQ", "SVINSSEQ", "SV_ID"
        }
        );

      for (int s = 0; s < static_cast<int>(info_semicolons.size()) - 1; ++s)
      {
        std::string const info_key_value = get_string_at_tab_index(info, info_semicolons, s);
        auto eq_it = std::find(info_key_value.begin(), info_key_value.end(), '=');

        if (eq_it == info_key_value.end())
        {
          if (keys_to_parse.count(info_key_value) == 1)
            new_var.infos[info_key_value] = "";

          continue;
        }

        std::string const key(info_key_value.begin(), eq_it);

        if (keys_to_parse.count(key) == 1)
          new_var.infos[key] = std::string(eq_it + 1, info_key_value.end());
      }
    }
  }

  // Parse samples, if any
  if (!SITES_ONLY && sample_names.size() > 0)
  {
    std::string const format = get_string_at_tab_index(line, tabs, 8);
    std::vector<std::size_t> all_format_colon = get_all_pos(format, ':');
    int ad_field = -1;
    int gt_field = -1;
    int pl_field = -1;
    int md_field = -1;
    int ra_field = -1;
    int pp_field = -1;
    int ft_field = -1;

    for (int32_t f = 0; f < static_cast<int>(all_format_colon.size() - 1); ++f)
    {
      std::string const field = get_string_at_tab_index(format, all_format_colon, f);

      if (field == "AD")
        ad_field = f;
      else if (field == "GT")
        gt_field = f;
      else if (field == "PL")
        pl_field = f;
      else if (field == "MD")
        md_field = f;
      else if (field == "RA")
        ra_field = f;
      else if (field == "PP")
        pp_field = f;
      else if (field == "FT")
        ft_field = f;
    }

    assert(ad_field != -1);
    assert(gt_field != -1);
    assert(pl_field != -1);
    int constexpr FIELD_OFFSET = 9;

    for (int i = FIELD_OFFSET; i < static_cast<int>(sample_names.size()) + FIELD_OFFSET; ++i)
    {
      // Create a new sample call
      SampleCall new_call;

      // Parse string of sample i
      std::string sample_string = get_string_at_tab_index(line, tabs, i);
      std::vector<std::size_t> sample_string_colon = get_all_pos(sample_string, ':');

      // Parse GT (for phase)
      std::string const gt_str = get_string_at_tab_index(sample_string, sample_string_colon, gt_field);

      // Check if the genotypes are phased
      /*
      if (std::find(gt_str.begin(), gt_str.end(), '|') != gt_str.end())
      {
        // Check if gt is missing
        if (gt_str[0] == '.')
        {
          new_var.phase.push_back(0);
        }
        else
        {
          std::vector<std::size_t> gt_str_pos = get_all_pos(gt_str, '|');
          assert(gt_str_pos.size() == 3);
          uint8_t call1 = static_cast<uint8_t>(std::stoul(get_string_at_tab_index(gt_str, gt_str_pos, 0)));
          uint8_t call2 = static_cast<uint8_t>(std::stoul(get_string_at_tab_index(gt_str, gt_str_pos, 1)));

          if (call1 <= call2)
            new_var.phase.push_back(0);
          else
            new_var.phase.push_back(1);
        }
      }
      */

      // Parse AD
      std::string const ad_str = get_string_at_tab_index(sample_string, sample_string_colon, ad_field);
      std::vector<std::size_t> ad_str_comma = get_all_pos(ad_str, ',');

      for (int j = 0; j < static_cast<int>(ad_str_comma.size()) - 1; ++j)
        new_call.coverage.push_back(static_cast<uint16_t>(std::stoul(get_string_at_tab_index(ad_str, ad_str_comma,
                                                                                             j))));

      // Parse MD
      if (md_field != -1)
      {
        std::string const md_str = get_string_at_tab_index(sample_string, sample_string_colon, md_field);
        unsigned long const md = std::stoul(md_str);
        assert(md <= 0xFFu);
        new_call.ambiguous_depth = static_cast<uint8_t>(md);
      }

      // Parse RA
      if (ra_field != -1)
      {
        std::string const ra_str = get_string_at_tab_index(sample_string, sample_string_colon, ra_field);
        std::vector<std::size_t> ra_str_comma = get_all_pos(ra_str, ',');

        assert(ra_str_comma.size() == 3);
        new_call.ref_total_depth = static_cast<uint16_t>(std::stoul(get_string_at_tab_index(ra_str, ra_str_comma, 0)));
        new_call.alt_total_depth = static_cast<uint16_t>(std::stoul(get_string_at_tab_index(ra_str, ra_str_comma, 1)));
      }

      // Parse PP
      if (pp_field != -1)
      {
        std::string const pp_str = get_string_at_tab_index(sample_string, sample_string_colon, pp_field);
        new_call.alt_proper_pair_depth = static_cast<uint8_t>(std::stoul(std::move(pp_str)));
      }

      // Parse FT
      if (ft_field != -1)
      {
        std::string const ft_str = get_string_at_tab_index(sample_string, sample_string_colon, ft_field);

        if (ft_str == "PASS")
        {
          new_call.filter = 0;
        }
        else
        {
          new_call.filter = static_cast<int8_t>(std::atoi(ft_str.substr(4).c_str()));
        }
      }

      // Parse PL
      std::string const pl_str = get_string_at_tab_index(sample_string, sample_string_colon, pl_field);
      std::vector<std::size_t> pl_str_comma = get_all_pos(pl_str, ',');

      for (int j = 0; j < static_cast<int>(pl_str_comma.size()) - 1; ++j)
        new_call.phred.push_back(static_cast<uint8_t>(std::stoul(get_string_at_tab_index(pl_str, pl_str_comma, j))));

      assert(new_call.coverage.size() * (new_call.coverage.size() + 1) / 2 == new_call.phred.size());
      new_var.calls.push_back(std::move(new_call));
    }
  }

  variants.push_back(std::move(new_var));
  return true;
}


void
Vcf::read_samples()
{
  if (!bgzf_in.rdbuf()->is_open())
    return;

  bool const is_checking_contigs = gyper::graph.contigs.size() == 0ull;

  while (true)
  {
    std::string line;

    if (!std::getline(bgzf_in, line))
    {
      BOOST_LOG_TRIVIAL(error) << "[vcf] Could not find any line with samples in '"
                               << filename << "'.";
      std::exit(1);
    }

    assert(line.size() > 2);

    // Read contigs
    if (is_checking_contigs && line[0] == '#' && line[1] == '#' &&
        line.substr(2, 11) == std::string("contig=<ID=")
        )
    {
      std::size_t comma_pos = line.find(',', 13);

      if (comma_pos == std::string::npos)
        continue;

      std::size_t closed_pos = line.find('>', comma_pos);
      assert(closed_pos != std::string::npos);
      Contig contig;
      contig.name = line.substr(13, comma_pos - 13);
      std::string length = line.substr(comma_pos + 8, closed_pos - comma_pos - 8);
      contig.length = static_cast<uint32_t>(std::stoul(length));
      gyper::graph.contigs.push_back(std::move(contig));
    }
    else if (line[0] == '#' && line[1] != '#')
    {
      // Gather list of all samples using this line
      assert(sample_names.size() == 0);
      std::vector<std::size_t> const all_line_tabs = get_all_pos(line, '\t');

      // Add all samples
      for (int i = 9; i < static_cast<int>(all_line_tabs.size()) - 1; ++i)
        sample_names.push_back(get_string_at_tab_index(line, all_line_tabs, i));

      break;
    }
  }

  // Recalculate contig offsets since they may have been changed
  if (is_checking_contigs)
    absolute_pos.calculate_offsets(gyper::graph.contigs);
}


void
Vcf::read(bool const SITES_ONLY)
{
  open_vcf_file_for_reading();
  read_samples();

  // Read all records
  while (read_record(SITES_ONLY))
  {}
  // The while loop will stop when all records have been read

  close_vcf_file();
}


void
Vcf::open_for_writing()
{
  switch (filemode)
  {
  case WRITE_UNCOMPRESSED_MODE:
    bgzf_stream.open("-", "w");
    break;

  case WRITE_BGZF_MODE:
    bgzf_stream.open(filename, "wb");
    break;

  default:
    BOOST_LOG_TRIVIAL(error) << "Trying to write in reading mode.";
    std::exit(1);
  }
}


void
Vcf::write_header()
{
  // Basic info
  bgzf_stream << "##fileformat=VCFv4.2\n"
              << "##fileDate=" << current_date() << "\n"
              << "##source=Graphtyper\n"
              << "##graphtyperVersion=" << graphtyper_VERSION_MAJOR << "." << graphtyper_VERSION_MINOR;

  if (std::string(GIT_NUM_DIRTY_LINES) != std::string("0"))
    bgzf_stream << "-dirty";

  bgzf_stream << "\n"
              << "##graphtyperGitBranch=" << GIT_BRANCH << '\n'
              << "##graphtyperSHA1=" << GIT_COMMIT_LONG_HASH << '\n';

  // Definitions of contigs
  for (auto const & contig : graph.contigs)
    bgzf_stream << "##contig=<ID=" << contig.name << ",length=" << contig.length << ">\n";

  // INFO definitions
  {
    bgzf_stream
      << "##INFO=<ID=ABHet,Number=1,Type=Float,Description=\"Allele Balance for heterozygous"
       "calls (read count of call2/(call1+call2)) where the called genotype is call1/call2. "
       "-1 if no heterozygous calls.\">\n"
      << "##INFO=<ID=ABHom,Number=1,Type=Float,Description=\"Allele Balance for homozygous calls"
       "(read count of A/(A+O)) where A is the called allele and O is anything else. -1 if no "
       "homozygous calls.\">\n"
      << "##INFO=<ID=ABHetMulti,Number=A,Type=Float,Description=\"List of Allele Balance values for"
       "multiallelic heterozygous calls (alt/(ref+alt)). Each value corresponds to an alt in the "
       "same order as they appear. -1 if not available.\">\n"
      << "##INFO=<ID=ABHomMulti,Number=R,Type=Float,Description=\"List of Allele Balance values for"
       "multiallelic homozygous calls (A/(A+0)) where A is the called allele and O is anything "
       "else. Each value corresponds to a ref or alt in the same order as they appear. -1 if not "
       "available.\">\n"
      << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Number of alternate alleles in called genotypes.\">\n"
      << "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Number of alleles in called genotypes.\">\n"
      << "##INFO=<ID=CR,Number=1,Type=Integer,Description=\"Number of clipped reads by Graphtyper.\">\n"
      << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of an SV.\">\n"
      << "##INFO=<ID=HOMSEQ,Number=.,Type=String,Description=\"Sequence of base pair identical "
       "homology at event breakpoints.\">\n"
      << "##INFO=<ID=INV3,Number=0,Type=Flag,Description=\"Inversion breakends open 3' of reported location\">\n"
      << "##INFO=<ID=INV5,Number=0,Type=Flag,Description=\"Inversion breakends open 5' of reported location\">\n"
      << "##INFO=<ID=LEFT_SVINSSEQ,Number=.,Type=String,Description=\"Known left side of insertion "
       "for an insertion of unknown length.\">\n"
      << "##INFO=<ID=MaxAAS,Number=A,Type=Integer,Description=\"Maximum alternative allele support "
       "per alt. allele.\">\n"
      << "##INFO=<ID=MaxAASR,Number=A,Type=Float,Description=\"Maximum alternative allele support "
       "ratio per alt. allele.\">\n"
      << "##INFO=<ID=MaxAltPP,Number=1,Type=Integer,Description=\"Maximum number of proper pairs "
       "support the alternative allele.\">\n"
      << "##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"Root-mean-square mapping quality.\">\n"
      << "##INFO=<ID=MQsquared,Number=1,Type=Integer,Description=\"Sum of squared mapping qualities. "
       "Used to calculate MQ.\">\n"
      << "##INFO=<ID=NCLUSTERS,Number=1,Type=Integer,Description=\"Number of SV candidates in "
       "cluster, as reported by popSVar.\">\n"
      << "##INFO=<ID=NGT,Number=3,Type=Integer,Description=\"Number of REF/REF, REF/ALT and ALT/ALT"
       "genotypes, respectively.\">\n"
      << "##INFO=<ID=NHet,Number=1,Type=Integer,Description=\"Number of heterozygous genotype "
       "calls.\">\n"
      << "##INFO=<ID=NHomRef,Number=1,Type=Integer,Description=\"Number of homozygous reference genotype calls.\">\n"
      << "##INFO=<ID=NHomAlt,Number=1,Type=Integer,Description=\"Number of homozygous alternative genotype calls.\">\n"
      << "##INFO=<ID=NRP,Number=0,Type=Flag,Description=\"Set if no Non-Reference has PASS.\">\n"
      << "##INFO=<ID=NUM_MERGED_SVS,Number=1,Type=Integer,Description=\"Number of SVs merged.\">\n"
      << "##INFO=<ID=OLD_VARIANT_ID,Number=1,Type=String,Description=\"Variant ID from a VCF (SVs only).\">\n"
      << "##INFO=<ID=ORSTART,Number=1,Type=Integer,Description=\"Start coordinate of sequence origin.\">\n"
      << "##INFO=<ID=OREND,Number=1,Type=Integer,Description=\"End coordinate of sequence origin.\">\n"
      << "##INFO=<ID=QD,Number=1,Type=Float,Description=\"QUAL divided by NonReferenceSeqDepth.\">\n"
      << "##INFO=<ID=PASS_AC,Number=A,Type=Integer,Description=\"Number of alternate alleles in "
       "called genotyped that have FT = PASS.\">\n"
      << "##INFO=<ID=PASS_AN,Number=1,Type=Integer,Description=\"Number of genotype calls that have"
       "FT = PASS.\">\n"
      << "##INFO=<ID=PASS_ratio,Number=1,Type=Float,Description=\"Ratio of genotype calls that have"
       "FT = PASS.\">\n"
      << "##INFO=<ID=PS,Number=1,Type=Integer,Description=\"Unique ID of the phase set this variant"
       "is a member of. If the calls are unphased, it is used to represent which haplotype set "
       "the variant is a member of.\">\n"
//      << "##INFO=<ID=RACount,Number=R,Type=Integer,Description=\"Number of realigned reads which "
//       "were uniquely aligned to this allele.\">\n"
//      << "##INFO=<ID=RADist,Number=R,Type=Integer,Description=\"Total realignment distance of all "
//       "reads per allele. Realignment distance = abs((Original aligned begin position) - "
//       "(New aligned begin position)).\">\n"
//      << "##INFO=<ID=RAMeanDist,Number=R,Type=Integer,Description=\"Mean realignment distance per "
//       "allele.\">\n"
      << "##INFO=<ID=RefLen,Number=1,Type=Integer,Description=\"Length of the reference "
       "allele.\">\n"
      << "##INFO=<ID=RELATED_SV_ID,Number=1,Type=Integer,Description=\"GraphTyper ID of a related "
       "SV.\">\n"
      << "##INFO=<ID=RIGHT_SVINSSEQ,Number=.,Type=String,Description=\"Known right side of "
       "insertion for an insertion of unknown length.\">\n"
      << "##INFO=<ID=SB,Number=1,Type=Float,Description=\"Strand bias (F/(F+R)) where F and R are "
       "forward and reverse strands, respectively. -1 if not available.\">\n"
      << "##INFO=<ID=SBAlt,Number=1,Type=Float,Description=\"Strand bias of alternative alleles only. "
       "-1 if not available.\">\n"
      << "##INFO=<ID=SBF,Number=R,Type=Integer,Description=\"Number of forward stranded reads per "
       "allele.\">\n"
      << "##INFO=<ID=SBF1,Number=R,Type=Integer,Description=\"Number of first forward stranded "
       "reads per allele.\">\n"
      << "##INFO=<ID=SBF2,Number=R,Type=Integer,Description=\"Number of second forward stranded "
       "reads per allele.\">\n"
      << "##INFO=<ID=SBR,Number=R,Type=Integer,Description=\"Number of reverse stranded reads per "
       "allele.\">\n"
      << "##INFO=<ID=SBR1,Number=R,Type=Integer,Description=\"Number of first reverse stranded "
       "reads per allele.\">\n"
      << "##INFO=<ID=SBR2,Number=R,Type=Integer,Description=\"Number of second reverse stranded "
       "reads per allele.\">\n"
      << "##INFO=<ID=SEQ,Number=1,Type=String,Description=\"Inserted sequence at variant site.\">\n"
      << "##INFO=<ID=SeqDepth,Number=1,Type=Integer,Description=\"Total accumulated sequencing "
       "depth over all the samples.\">\n"
      << "##INFO=<ID=SV_ID,Number=1,Type=Integer,Description=\"GraphTyper's ID on SV\">\n"
      << "##INFO=<ID=SVINSSEQ,Number=.,Type=String,Description=\"Sequence of insertion\">\n"
      << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of structural variant in bp."
       " Negative lengths indicate a deletion.\">\n"
      << "##INFO=<ID=SVMODEL,Number=1,Type=String,Description=\"Model used for SV genotyping.\">\n"
      << "##INFO=<ID=SVSIZE,Number=1,Type=Integer,Description=\"Size of structural variant in bp."
       " Always 50 or more.\">\n"
      << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant.\">\n"
      << "##INFO=<ID=VarType,Number=1,Type=String,Description=\"First letter is program identifier,"
       "the second letter is variant type.\">\n";
  }

  // FORMAT definitions
  {
    bgzf_stream
      << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GenoType call. ./. is called if there is no "
       "coverage at the variant site.\">\n"
      << "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Filter. PASS or FAILN where N is a number.\">\n"
      << "##FORMAT=<ID=AD,Number=R,Type=Integer,Description="
       "\"Allelic depths for the ref and alt alleles in the order listed.\">\n"
      << "##FORMAT=<ID=MD,Number=1,Type=Integer,Description=\"Read depth of multiple alleles.\">\n"
      << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth.\">\n"
      << "##FORMAT=<ID=RA,Number=2,Type=Integer,Description="
       "\"Total read depth of the reference allele and all alternative alleles, "
       "including reads that support more than one allele.\">\n"
      << "##FORMAT=<ID=PP,Number=1,Type=Integer,Description=\"Number of reads that "
       "support non-reference haplotype that are proper pairs.\">\n"
      << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality.\">\n"
      << "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"PHRED-scaled genotype likelihoods.\">\n";
  }

  // FILTER definitions
  {
    bgzf_stream << "##FILTER=<ID=LowABHet,Description=\"Allele balance of heterozygous carriers is below 17.5%.\">\n"
                << "##FILTER=<ID=LowABHom,Description=\"Allele balance of homozygous carriers is below 90%.\">\n"
                << "##FILTER=<ID=LowQD,Description=\"QD (quality by depth) is below 9.0.\">\n"
                << "##FILTER=<ID=LowQUAL,Description=\"QUAL score is less than 10.\">\n"
                << "##FILTER=<ID=LowPratio,Description=\"Ratio of PASSed calls was too low.\">\n";

  }

  // Column names
  bgzf_stream << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";

  if (sample_names.size() > 0)
  {
    // Only a "format" column if there are any samples
    bgzf_stream << "\tFORMAT";

    for (auto const & sample_name : sample_names)
      bgzf_stream << "\t" << sample_name;
  }

  bgzf_stream << "\n";
}


void
Vcf::write_record(Variant const & var, std::string const & suffix, bool const FILTER_ZERO_QUAL)
{
  // Parse the position
  auto contig_pos = absolute_pos.get_contig_position(var.abs_pos, gyper::graph.contigs);

  if (!Options::instance()->output_all_variants && var.calls.size() > 0 && var.seqs.size() > 200)
  {
    BOOST_LOG_TRIVIAL(info) << "Skipped outputting variant at position "
                            << contig_pos.first << ":" << contig_pos.second
                            << " because there are " << var.seqs.size() << " alleles.";
    return;
  }

  if (!Options::instance()->output_all_variants)
  {
    std::size_t total_allele_length = 0;

    for (auto const & seq : var.seqs)
    {
      total_allele_length += seq.size();

      if (total_allele_length > 16000)
      {
        BOOST_LOG_TRIVIAL(warning) << "Skipped outputting variant at position "
                                   << contig_pos.first << ":" << contig_pos.second
                                   << " because the total length of alleles is too high.";
        return;
      }
    }
  }

  // Calculate qual
  const uint64_t variant_qual = var.get_qual();

  if (FILTER_ZERO_QUAL && (sample_names.size() > 0 && variant_qual == 0))
  {
    //BOOST_LOG_TRIVIAL(warning) << "Zero qual variant at pos " << contig_pos.first << ":" << contig_pos.second;
    BOOST_LOG_TRIVIAL(debug) << "Skipped outputting variant at position "
                             << contig_pos.first << ":" << contig_pos.second
                             << " because the variant had zero quality.";
    return;
  }

  bgzf_stream << contig_pos.first << '\t';
  bgzf_stream << contig_pos.second << '\t';

  // Write the ID field
  bgzf_stream << contig_pos.first; // Keep the 'chr'

  bgzf_stream << ':' << contig_pos.second << ':' << var.determine_variant_type();

  if (var.suffix_id.size() > 0)
    bgzf_stream << "[" << var.suffix_id << "]";

  bgzf_stream << suffix;

  // Parse the sequences
  assert(var.seqs.size() >= 2);
  bgzf_stream << '\t' << std::string(var.seqs[0].begin(), var.seqs[0].end()) << '\t'
              << std::string(var.seqs[1].begin(), var.seqs[1].end());

  // Print other allele sequences if it is multi-allelic marker
  for (long a = 2; a < static_cast<long>(var.seqs.size()); ++a)
    bgzf_stream << ',' << std::string(var.seqs[a].begin(), var.seqs[a].end());

  // Parse qual
  bgzf_stream << "\t" << std::to_string(variant_qual) << "\t";

  // Parse filter
  if (sample_names.size() == 0 || Options::const_instance()->ploidy > 2)
  {
    bgzf_stream << ".\t";
  }
  else
  {
    bool is_pass = true;

    if (var.infos.count("ABHet") == 1 &&
        var.infos.at("ABHet") != std::string("-1") &&
        std::stod(var.infos.at("ABHet")) < 0.175)
    {
      if (!is_pass)
        bgzf_stream << ";";

      bgzf_stream << "LowABHet";
      is_pass = false;
    }

    if (var.infos.count("ABHom") == 1 &&
        var.infos.at("ABHom") != std::string("-1") &&
        std::stod(var.infos.at("ABHom")) < 0.85)
    {
      if (!is_pass)
        bgzf_stream << ";";

      bgzf_stream << "LowABHom";
      is_pass = false;
    }

    if (var.infos.count("AN") == 1 && std::stoi(var.infos.at("AN")) >= 6 &&
        var.infos.count("QD") == 1 && std::stod(var.infos.at("QD")) < 9.0)
    {
      if (!is_pass)
        bgzf_stream << ";";

      bgzf_stream << "LowQD";
      is_pass = false;
    }

    if (variant_qual < 10)
    {
      if (!is_pass)
        bgzf_stream << ";";

      bgzf_stream << "LowQUAL";
      is_pass = false;
    }

    // Only filter on PASS_ratio if we have sufficient amount of samples
    if (var.infos.count("AN") == 1 && std::stoi(var.infos.at("AN")) >= 500 &&
        var.infos.count("PASS_ratio") == 1 && std::stod(var.infos.at("PASS_ratio")) < 0.05
        )
    {
      if (!is_pass)
        bgzf_stream << ";";

      bgzf_stream << "LowPratio";
      is_pass = false;
    }

    if (is_pass)
      bgzf_stream << "PASS";

    bgzf_stream << "\t";
  }


  // Parse info
  if (var.infos.empty())
  {
    bgzf_stream << ".";
  }
  else
  {
    auto write_info = [&](std::map<std::string, std::string>::const_iterator it)
                      {
                        bgzf_stream << it->first;

                        if (it->second.size() > 0)
                          bgzf_stream << '=' << it->second;
                      };

    write_info(var.infos.cbegin());

    for (auto map_it = std::next(var.infos.cbegin(), 1); map_it != var.infos.cend(); ++map_it)
    {
      bgzf_stream << ';';
      write_info(map_it);
    }

  }

  assert(sample_names.size() == var.calls.size());
  bool const is_sv = var.is_sv();

  // Parse FORMAT
  if (var.calls.size() > 0)
  {
    if (is_sv)
      bgzf_stream << "\tGT:FT:AD:MD:DP:RA:PP:GQ:PL";
    else
      bgzf_stream << "\tGT:AD:MD:DP:GQ:PL";
  }

  for (long i = 0; i < static_cast<long>(var.calls.size()); ++i)
  {
    auto const & call = var.calls[i];

    // Write GT
    auto pl_non_zero = [](uint8_t const pl) -> bool
                       {
                         return pl != 0;
                       };

    if (std::find_if(call.phred.begin(), call.phred.end(), pl_non_zero) == call.phred.end())
    {
      bgzf_stream << "\t./.";
    }
    else
    {
      std::pair<uint16_t, uint16_t> const gt_call = call.get_gt_call();
      bgzf_stream << "\t" << gt_call.first << "/" << gt_call.second;
    }

    // Write FT
    long const gq = call.get_gq();

    if (is_sv)
    {
      bgzf_stream << ":";
      long filter = call.check_filter(gq);

      if (filter == 0)
      {
        bgzf_stream << "PASS";
      }
      else
      {
        assert(filter > 0);
        bgzf_stream << "FAIL" << filter;
      }
    }

    // Write AD
    assert(call.coverage.size() > 0);
    bgzf_stream << ":" << call.coverage[0];

    for (auto ad_it = call.coverage.begin() + 1; ad_it != call.coverage.end(); ++ad_it)
      bgzf_stream << "," << *ad_it;

    // Write MD (Multi-depth)
    bgzf_stream << ":" << static_cast<uint64_t>(call.ambiguous_depth);

    // Write DP
    bgzf_stream << ":" << call.get_depth();

    // Write RA
    if (is_sv)
    {
      bgzf_stream << ":" << call.ref_total_depth << "," << call.alt_total_depth;
    }

    // Write PP
    if (is_sv)
    {
      bgzf_stream << ':' << static_cast<std::size_t>(call.alt_proper_pair_depth);
    }

    // Write GQ
    bgzf_stream << ':' << std::min(99l, gq);

    // Write PL
    bgzf_stream << ':' << static_cast<uint16_t>(call.phred[0]);
    // This cast to uint16_t is needed! Otherwise uint8_t is represented as a char

    for (std::size_t p = 1; p < call.phred.size(); ++p)
      bgzf_stream << ',' << static_cast<uint16_t>(call.phred[p]);
  }

  // Fin.
  bgzf_stream << '\n';
}


void
Vcf::write(std::string const & region)
{
  this->open_for_writing();
  this->write_header();
  this->write_records(region);
  this->close_vcf_file();
}


void
Vcf::write_records(uint32_t const region_begin,
                   uint32_t const region_end,
                   bool const FILTER_ZERO_QUAL,
                   std::vector<Variant> const & vars
                   )
{
  if (vars.size() == 0)
    return;

  // Sort the variants
  auto compare_variants = [&](long const i, long const j) -> bool
                          {
                            assert(i < static_cast<long>(vars.size()));
                            assert(j < static_cast<long>(vars.size()));
                            gyper::Variant const & a = vars[i];
                            gyper::Variant const & b = vars[j];
                            assert(a.seqs.size() >= 2);
                            assert(b.seqs.size() >= 2);
                            return a.abs_pos < b.abs_pos ||
                                   (a.abs_pos == b.abs_pos &&
                                    a.determine_variant_type() < b.determine_variant_type()
                                   ) ||
                                   (a.abs_pos == b.abs_pos &&
                                    a.determine_variant_type() == b.determine_variant_type() &&
                                    a.seqs < b.seqs
                                   );
                          };

  auto same_variant_id = [](Variant const & a, Variant const & b) -> bool
                         {
                           return a.abs_pos == b.abs_pos && a.determine_variant_type() == b.determine_variant_type();
                         };

  auto inside_region = [region_begin, region_end](uint32_t const pos) -> bool
                       {
                         return pos >= region_begin && pos <= region_end;
                       };

  std::vector<long> indexes;

  for (long i = 0; i < static_cast<long>(vars.size()); ++i)
    indexes.push_back(i);

  std::sort(indexes.begin(), indexes.end(), compare_variants);
  assert(indexes.size() > 0);
  assert(indexes[0] < static_cast<long>(vars.size()));

  {
    auto const & first_var = vars[indexes[0]];

    if (inside_region(first_var.abs_pos))
      write_record(first_var, "" /*suffix*/, FILTER_ZERO_QUAL);
  }

  long dup = -1; // -1 means no duplication

  // Write the variants in the correct order
  for (long i = 1; i < static_cast<long>(indexes.size()); ++i)
  {
    assert(i - 1 >= 0);
    assert(indexes[i - 1] >= 0);
    assert(i < static_cast<long>(vars.size()));
    assert(indexes[i] < static_cast<long>(vars.size()));
    auto const & prev_var = vars[indexes[i - 1]];
    auto const & curr_var = vars[indexes[i]];

    // stop if the we have gone beyond the region
    if (curr_var.abs_pos > region_end)
      break;

    // skip variants before the region
    if (curr_var.abs_pos < region_begin)
      continue;

    // skip duplicate variants
    if (curr_var == prev_var)
      continue;

    // make sure the variant ID is unique
    if (!same_variant_id(curr_var, prev_var))
    {
      write_record(curr_var, "" /*suffix*/, FILTER_ZERO_QUAL);
      dup = -1;
    }
    else
    {
      ++dup;
      assert(dup >= 0);
      write_record(curr_var, std::string(".") + std::to_string(dup), FILTER_ZERO_QUAL);
    }
  }
}


void
Vcf::write_records(std::string const & region, bool const FILTER_ZERO_QUAL)
{
  uint32_t region_begin = 0;
  uint32_t region_end = 0xFFFFFFFFull;

  // Restrict to a region if it is given
  if (region != ".")
  {
    GenomicRegion genomic_region(region);

    if (absolute_pos.is_contig_available(genomic_region.chr))
    {
      region_begin = 1 + absolute_pos.get_absolute_position(genomic_region.chr,
                                                            genomic_region.begin
                                                            );

      region_end = absolute_pos.get_absolute_position(genomic_region.chr,
                                                      genomic_region.end
                                                      );
    }
  }

  this->write_records(region_begin, region_end, FILTER_ZERO_QUAL, variants);
}


void
Vcf::write_segments()
{
  BOOST_LOG_TRIVIAL(debug) << "[graphtyper::vcf] Writing "
                           << segments.size()
                           << " segments to "
                           << filename;

  for (auto const & segment : segments)
  {
    assert(sample_names.size() == segment.segment_calls.size());
    auto contig_pos = absolute_pos.get_contig_position(segment.id, gyper::graph.contigs);

    // Write CHROM and POS
    bgzf_stream << contig_pos.first << "\t" << contig_pos.second;

    // Write the ID
    bgzf_stream << "\t" << contig_pos.first << ":" << contig_pos.second << ":" << segment.var_type;

    if (segment.segment_name.size() > 0)
      bgzf_stream << ":" << segment.segment_name;
    else if (segment.extra_id != -1)
      bgzf_stream << ":" << std::to_string(segment.extra_id);

    // Write the sequences
    bgzf_stream << "\t" << segment.get_ref_string() << "\t" << segment.get_alt_string();

    // Parse QUAL
    bgzf_stream << "\t" << ".";

    // Parse FILTER
    bgzf_stream << "\t" << ".";

    // Parse INFO
    bgzf_stream << "\t" << "RefLen=" << std::to_string(segment.ref_size);

    // Parse FORMAT
    bgzf_stream << "\t" << "GT:PL";

    for (auto const & segment_call : segment.segment_calls)
    {
      bgzf_stream << "\t" << segment_call.call.first << "/" << segment_call.call.second << ":";

      assert(segment_call.phred.size() > 1);
      bgzf_stream << std::to_string(static_cast<uint16_t>(segment_call.phred[0]));

      for (unsigned p = 1; p < segment_call.phred.size(); ++p)
        bgzf_stream << "," << std::to_string(static_cast<uint16_t>(segment_call.phred[p]));
    }

    bgzf_stream << "\n";
  }
}


void
Vcf::close_vcf_file()
{
  bgzf_stream.close();

  if (bgzf_in)
  {
    bgzf_in.close();
  }
}


void
Vcf::write_tbi_index() const
{
  int ret = tbx_index_build(filename.c_str(), 0, &tbx_conf_vcf);

  if (ret < 0)
  {
    BOOST_LOG_TRIVIAL(error) << "Could not build VCF index";
    std::exit(ret);
  }
}


void
Vcf::clear()
{
  sample_names.clear();
  variants.clear();
}


void
Vcf::add_segment(Segment && segment)
{
  assert(segment.allele_names.size() > 1);

  if (segment.allele_names.size() == 2 /*check if it is already biallelic*/)
  {
    this->segments.push_back(std::move(segment));
    return;
  }

  // Split the segments into biallelic VCF records
  std::vector<Segment> new_segments = segment.get_biallelic_segments();
  segment.clear();
  std::move(new_segments.begin(), new_segments.end(), std::back_inserter(this->segments));
}


void
Vcf::add_haplotype(Haplotype & haplotype, bool const clear_haplotypes, uint32_t const phase_set)
{
  assert(haplotype.gts.size() > 0);

  // Add each genotype
  std::vector<Variant> new_vars;
  new_vars.reserve(haplotype.gts.size());

  for (auto const & gt : haplotype.gts)
    new_vars.push_back(Variant(gt));

  assert(new_vars.size() == haplotype.gts.size());
  assert(new_vars.size() == haplotype.var_stats.size());

  // Add variant stats
  for (long i = 0; i < static_cast<long>(new_vars.size()); ++i)
  {
    auto const & stat = haplotype.var_stats[i];
    auto & new_var = new_vars[i];

    new_var.infos["CR"] = std::to_string(stat.clipped_reads);
    new_var.infos["PS"] = std::to_string(phase_set);
    new_var.infos["MQsquared"] = std::to_string(stat.mapq_squared);
    new_var.infos["SBF"] = stat.get_forward_strand_bias();
    new_var.infos["SBR"] = stat.get_reverse_strand_bias();
    new_var.infos["SBF1"] = stat.get_r1_forward_strand_bias();
    new_var.infos["SBF2"] = stat.get_r2_forward_strand_bias();
    new_var.infos["SBR1"] = stat.get_r1_reverse_strand_bias();
    new_var.infos["SBR2"] = stat.get_r2_reverse_strand_bias();
  }

  //long const ploidy = Options::const_instance()->ploidy;

  //if (ploidy <= 2)
  {
    // Normal genotyping
    for (auto & hap_sample : haplotype.hap_samples)
    {
      std::vector<std::vector<uint8_t> > gt_phred =
        get_genotype_phred(hap_sample, haplotype.gts);

      assert(gt_phred.size() == new_vars.size());
      assert(hap_sample.gt_coverage.size() == new_vars.size());

      for (long i = 0; i < static_cast<long>(new_vars.size()); ++i)
      {
        assert(i < static_cast<long>(gt_phred.size()));
        assert(i < static_cast<long>(hap_sample.gt_coverage.size()));
        assert(new_vars[i].seqs.size() == hap_sample.gt_coverage[i].size());

        SampleCall new_sample_call(std::move(gt_phred[i]),
                                   std::move(hap_sample.gt_coverage[i]),
                                   hap_sample.get_ambiguous_depth(),
                                   hap_sample.get_ambiguous_depth_alt(),
                                   hap_sample.get_alt_proper_pair_depth());

        new_vars[i].calls.push_back(std::move(new_sample_call)); // Add new call
      }
    }
  }
  /*
  else
  {
    // Camou genotyping
    for (auto & hap_sample : haplotype.hap_samples)
    {
      for (long i = 0; i < static_cast<long>(new_vars.size()); ++i)
      {
        assert(i < static_cast<long>(hap_sample.gt_coverage.size()));

        auto const & cov = hap_sample.gt_coverage[i];
        assert(cov.size() > 0);
        long const alt_coverage = std::accumulate(cov.begin() + 1, cov.end(), 0l);
        long const total_coverage = alt_coverage + cov[0];
        long const cnum = cov.size();
        std::vector<uint8_t> phred;

        if (total_coverage == 0)
        {
          // Not enough data
          phred.resize(cnum * (cnum + 1) / 2, 0);
        }
        else
        {
          phred.resize(cnum * (cnum + 1) / 2, 255);
          phred[0] = 0;
          std::vector<long> normalized_cov;
          normalized_cov.reserve(cnum);

          for (long k{0}; k < cnum; ++k)
            normalized_cov.push_back(cov[k] * ploidy / 2l);

          for (long y{1}; y < cnum; ++y)
          {
            long constexpr ERROR = 4;
            auto const norm_cov = normalized_cov[y];
            long phred00 = norm_cov * ERROR;
            long phred01_or_11 = cov[0];
            long const min_phred = std::min(phred00, phred01_or_11);
            phred00 -= min_phred;
            phred01_or_11 -= min_phred;

            phred00 = std::min(99l, phred00 * 3l);
            phred01_or_11 = std::min(99l, phred01_or_11 * 3l);

            if (phred00 > phred[0])
              phred[0] = phred00;

            for (long x{0}; x < cnum; ++x)
            {
              auto const index = to_index(x, y);

              if (phred01_or_11 < phred[index])
                phred[index] = phred01_or_11;
            }
          }
        }

        assert(new_vars[i].seqs.size() == hap_sample.gt_coverage[i].size());

        SampleCall new_sample_call(std::move(phred),
                                   std::move(hap_sample.gt_coverage[i]),
                                   hap_sample.get_ambiguous_depth(),
                                   hap_sample.get_ambiguous_depth_alt(),
                                   0); // propar pair

        new_vars[i].calls.push_back(std::move(new_sample_call)); // Add new call
      }
    }
  }
  //*/

  // Set variant suffix ID
  if (Options::instance()->variant_suffix_id.size() > 0)
  {
    for (auto & new_var : new_vars)
      new_var.suffix_id = Options::instance()->variant_suffix_id;
  }

  // Add the variants
  std::move(new_vars.begin(), new_vars.end(), std::back_inserter(variants));

  if (clear_haplotypes)
    haplotype.clear();
}


void
Vcf::add_haplotypes_for_extraction(std::vector<HaplotypeCall> const & hap_calls, bool const is_splitting_vars)
{
  assert(graph.size() > 0);
  bool const is_sv_graph = graph.is_sv_graph;

  for (long i = 0; i < static_cast<long>(hap_calls.size()); ++i)
  {
    auto & hap_call = hap_calls[i];
    auto const & gts = hap_call.gts;

    assert(hap_call.calls.size() >= 1);
    assert(hap_call.calls[0] == 0);

    // Only add variants if there is something else than the reference called
    if (hap_call.calls.size() >= 2)
    {
      Variant var(gts, hap_call.calls);

      if (is_sv_graph)
      {
        // Reformat SVs when we have an SV graph
        for (int a = 1; a < static_cast<int>(var.seqs.size()); ++a)
        {
          auto & seq = var.seqs[a];
          auto find_it = std::find(seq.cbegin(), seq.cend(), '<');

          if (std::distance(find_it, seq.cend()) > 11)
          {
            std::istringstream ss {
              std::string(find_it + 4, find_it + 11)
            };
            long sv_id;
            ss >> sv_id;

            // If we can't parse correctly the SV ID so we ignore it
            assert(ss.eof());
            assert(sv_id < static_cast<long>(graph.SVs.size()));

            auto const & sv = graph.SVs[sv_id];
            std::string const sv_allele = std::string(1, var.seqs[0][0]) + sv.get_allele();
            var.seqs[a] = std::vector<char>(sv_allele.cbegin(), sv_allele.cend());
            var.infos["SVTYPE"] = sv.get_type();
            var.infos["SVLEN"] = sv.type == DEL ?
                                 std::to_string(-sv.length) :
                                 std::to_string(sv.length);

            var.infos["SVSIZE"] = std::to_string(sv.length);

            if (sv.old_variant_id.size() > 0)
              var.infos["OLD_VARIANT_ID"] = sv.old_variant_id;
          }
        }

        this->variants.push_back(std::move(var));
      }
      else
      {
        // Split variants further down if they are large and the graph is not an SV graph
        long constexpr SPLIT_THRESHOLD = 25;

        if (is_splitting_vars && var.seqs.size() >= 3 && std::all_of(var.seqs.begin(),
                                                                     var.seqs.end(),
                                                                     [](std::vector<char> const & seq)
          {
            return static_cast<long>(seq.size()) > SPLIT_THRESHOLD + 25;
          })
            )
        {
          //this->variants.push_back(std::move(var));
          //*
          //std::cout << var.print() << "\n";
          paw::Skyr skyr(var.seqs);

          bool constexpr is_normalize {
            false
          };                                  // The events must not be normalized
          bool constexpr use_asterisks {
            false
          };                                   // We should not use asterisks
          skyr.find_all_edits(is_normalize);
          skyr.find_variants_from_edits();
          skyr.populate_variants_with_calls(use_asterisks);
          std::vector<paw::Variant> spl_vars = skyr.split_variants(SPLIT_THRESHOLD);

          if (spl_vars.size() == 1)
          {
            this->variants.push_back(std::move(var));
          }
          else
          {
            BOOST_LOG_TRIVIAL(debug) << "Variant to be splitted: " << var.print();

            for (auto & spl_var : spl_vars)
            {
              Variant new_var;
              new_var.abs_pos = var.abs_pos + spl_var.pos;
              bool is_any_zero_size = false;

              for (auto const & variant_seq : spl_var.seqs)
              {
                if (variant_seq.size() == 0)
                  is_any_zero_size = true;

                new_var.seqs.push_back(std::vector<char>(variant_seq.begin(), variant_seq.end()));
              }

              if (is_any_zero_size)
                new_var.add_base_in_front(true); // Add N is true

              BOOST_LOG_TRIVIAL(debug) << "...into " << new_var.print();
              this->variants.push_back(std::move(new_var));
            }
          }//*/
        }
        else
        {
          this->variants.push_back(std::move(var));
        }
      }
    }
  }
}


/*
void
Vcf::post_process_camou_variants(long const ploidy)
{
  // The variant is not camou
  if (ploidy <= 2)
    return;

  for (auto & var : variants)
  {

  }
}



void
Vcf::post_process_variants(bool const NORMALIZE, bool const TRIM_SEQUENCES)
{
  if (NORMALIZE)
  {
    for (auto & var : variants)
      var.normalize();
  }
  else if (TRIM_SEQUENCES)
  {
    for (auto & var : variants)
      var.trim_sequences(false); // Don't keep one match
  }
}*/


/** Non-member functions */
std::vector<std::size_t>
get_all_pos(std::string const & line, char const delim)
{
  std::vector<std::size_t> tabs(1, 0ul);

  while (true)
  {
    std::size_t pos = line.find(delim, tabs.back() + 1);

    if (pos == std::string::npos)
    {
      tabs.push_back(pos);
      return tabs;
    }
    else
    {
      ++pos;
      tabs.push_back(pos);
    }
  }
}


std::string
get_string_at_tab_index(std::string const & line,
                        std::vector<std::size_t> const & tabs,
                        int const index
                        )
{
  assert(index >= 0);
  assert(index + 1 < static_cast<int>(tabs.size()));
  assert(tabs[index + 1] > tabs[index]);

  if (tabs[index + 1] == std::string::npos)
    return line.substr(tabs[index], std::string::npos);
  else
    return line.substr(tabs[index], tabs[index + 1] - tabs[index] - 1);
}


} // namespace gyper
