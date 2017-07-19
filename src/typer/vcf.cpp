#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>

#include "InputFile.h" // Included in the StatGen library

#include <boost/algorithm/string/predicate.hpp> // boost::algorithm::ends_with
#include <boost/log/trivial.hpp>

#include <graphtyper/graph/absolute_position.hpp>
#include <graphtyper/graph/genomic_region.hpp>
#include <graphtyper/graph/var_record.hpp>
#include <graphtyper/typer/vcf.hpp>
#include <graphtyper/utilities/graph_help_functions.hpp>
#include <graphtyper/utilities/options.hpp> // gyper::options::instance()
#include <graphtyper/utilities/type_conversions.hpp>


namespace
{

bool
same_prefix(std::string const & s1, std::string const & s2)
{
  std::size_t const SMALLER_STRING_SIZE = std::min(s1.size(), s2.size());

  for (std::size_t i = 0; i < SMALLER_STRING_SIZE; ++i)
  {
    if (s1[i] != s2[i])
      return false;
  }

  return true;
}


std::string
current_date()
{
  time_t  now = time(0);
  struct tm tstruct;
  char buf[16];
  tstruct = *localtime(&now);
  strftime(buf, sizeof(buf), "%Y%m%d", &tstruct);
  return buf;
}


std::vector<uint8_t>
get_haplotype_phred(gyper::HapSample const & sample, uint32_t const cnum)
{
  using namespace gyper;
  assert(sample.log_score.size() > 0);
  std::vector<uint8_t> hap_phred;

  // Check if all phred scores are zero by finding the first non-zero phred score
  auto find_it = std::find_if(sample.log_score.begin(), sample.log_score.end(), [](uint16_t const val){
      return val != 0;
    });

  if (find_it == sample.log_score.end())
  {
    // If no non-zero phred score is found, the phred scores of the haplotype should be all zero as well
    hap_phred = std::vector<uint8_t>(cnum * (cnum + 1) / 2, 0u);
  }
  else
  {
    hap_phred = std::vector<uint8_t>(cnum * (cnum + 1) / 2, 255u);

    // First find out what the maximum log score is
    uint16_t const max_log_score = *std::max_element(sample.log_score.begin(), sample.log_score.end());
    std::size_t const num_alleles = cnum * (cnum + 1) / 2;
    assert(sample.log_score.size() == num_alleles);

    for (std::size_t i = 0; i < num_alleles; ++i)
    {
      double const LOG10_HALF_times_10 = 3.0102999566398119521373889472449302676818988146210854131;
      uint64_t const phred_score = std::llround((max_log_score - sample.log_score[i]) * LOG10_HALF_times_10);

      if (phred_score < 255u)
        hap_phred[i] = static_cast<uint8_t>(phred_score);
    }
  }

  return hap_phred;
}


std::vector<std::vector<uint8_t> >
get_genotype_phred(gyper::HapSample const & sample, std::vector<uint8_t> & phase, std::vector<gyper::Genotype> const & gts)
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

  phase.resize(gts.size(), 0);

  if (find_it == hap_phred.end())
  {
    for (unsigned i = 0; i < gts.size(); ++i)
      phred[i] = std::vector<uint8_t>(gts[i].num * (gts[i].num + 1) / 2, 0u);

    return phred;
  }
  else
  {
    for (unsigned i = 0; i < gts.size(); ++i)
      phred[i] = std::vector<uint8_t>(gts[i].num * (gts[i].num + 1) / 2, 255u);
  }

  for (uint32_t i = 0; i < hap_phred.size(); ++i)
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
      uint32_t const index = call1 > call2 ? to_index(call2, call1) : to_index(call1, call2);

      // Check if we need to flip the phasing compared to the first variant in the phase set
      if (call1 > call2 && hap_phred[i] == 0)
        phase[j] = 1;

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

/***************
 * CLASS ACESS *
 ***************/

std::size_t
Vcf::get_sample_index(std::string const & sample_name) const
{
  auto find_it = std::lower_bound(sample_names.begin(), sample_names.end(), sample_name);

  // If sample already exists, don't add it again
  if (find_it == sample_names.end() || *find_it != sample_name)
  {
    BOOST_LOG_TRIVIAL(error) << "[vcf] Could not find sample " << sample_name << ".";
    std::exit(1);
  }

  return std::distance(sample_names.begin(), find_it);
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
    vcf_file = ifopen(filename.c_str(), "r", InputFile::UNCOMPRESSED);
    break;

  case READ_BGZF_MODE:
    vcf_file = ifopen(filename.c_str(), "rb", InputFile::BGZF);
    break;

  default:
    std::cerr << "[vcf] ERROR: Trying to read in writing mode." << std::endl;
    std::exit(1);
  }
}


std::string
Vcf::read_line()
{
  std::string line;
  vcf_file >> line;
  return line;
}


bool
Vcf::read_record()
{
  std::string const line = read_line();

  if (line.size() == 0)
    return false;

  // Get all positions of TABs in the line
  std::vector<std::size_t> const tabs = get_all_pos(line, '\t');
  // We ignore the following fields: id (2), qual (5), filter (6), info (7)

  std::string const chrom = get_string_at_tab_index(line, tabs, 0);
  uint32_t const pos = std::stoul(get_string_at_tab_index(line, tabs, 1));
  std::string ref = get_string_at_tab_index(line, tabs, 3);
  std::string const alts = get_string_at_tab_index(line, tabs, 4);
  std::vector<std::size_t> const alt_commas = get_all_pos(alts, ',');

  Variant new_var; // Create a new variant for this position
  new_var.abs_pos = absolute_pos.get_absolute_position(chrom, pos); // Parse positions

  // Parse sequences
  new_var.seqs.push_back(gyper::to_vec(std::move(ref)));
  assert(alt_commas.size() >= 2);

  for (std::size_t a = 0; a < alt_commas.size() - 1; ++a)
    new_var.seqs.push_back(gyper::to_vec(get_string_at_tab_index(alts, alt_commas, a)));

  // Parse infos
  {
    std::string const info = get_string_at_tab_index(line, tabs, 7);

    // Don't parse anything if the INFO field is empty
    if (info.size() > 0)
    {
      std::vector<std::size_t> const info_semicolons = get_all_pos(info, ';');

      for (std::size_t s = 0; s < info_semicolons.size() - 1; ++s)
      {
        std::string const info_key_value = get_string_at_tab_index(info, info_semicolons, s);

        if (same_prefix(info_key_value, std::string("AC=")))
          new_var.infos["AC"] = std::string(info_key_value.begin() + 3, info_key_value.end()); // Get AC so vcf_break_down can first remove uncalled alleles
        else if (same_prefix(info_key_value, std::string("CR=")))
          new_var.infos["CR"] = std::string(info_key_value.begin() + 3, info_key_value.end());
        else if (same_prefix(info_key_value, std::string("GX=")))
          new_var.infos["GX"] = std::string(info_key_value.begin() + 3, info_key_value.end());
        else if (same_prefix(info_key_value, std::string("PS=")))
          new_var.infos["PS"] = std::string(info_key_value.begin() + 3, info_key_value.end());
        else if (same_prefix(info_key_value, std::string("MQ=")))
          new_var.infos["MQ"] = std::string(info_key_value.begin() + 3, info_key_value.end());
        else if (same_prefix(info_key_value, std::string("MQ0=")))
          new_var.infos["MQ0"] = std::string(info_key_value.begin() + 4, info_key_value.end());
        else if (same_prefix(info_key_value, std::string("MQperAllele=")))
          new_var.infos["MQperAllele"] = std::string(info_key_value.begin() + 12, info_key_value.end());
        else if (same_prefix(info_key_value, std::string("RACount=")))
          new_var.infos["RACount"] = std::string(info_key_value.begin() + 8, info_key_value.end());
        else if (same_prefix(info_key_value, std::string("RADist=")))
          new_var.infos["RADist"] = std::string(info_key_value.begin() + 7, info_key_value.end());
        else if (same_prefix(info_key_value, std::string("SBF=")))
          new_var.infos["SBF"] = std::string(info_key_value.begin() + 4, info_key_value.end());
        else if (same_prefix(info_key_value, std::string("SBF1=")))
          new_var.infos["SBF1"] = std::string(info_key_value.begin() + 5, info_key_value.end());
        else if (same_prefix(info_key_value, std::string("SBF2=")))
          new_var.infos["SBF2"] = std::string(info_key_value.begin() + 5, info_key_value.end());
        else if (same_prefix(info_key_value, std::string("SBR1=")))
          new_var.infos["SBR1"] = std::string(info_key_value.begin() + 5, info_key_value.end());
        else if (same_prefix(info_key_value, std::string("SBR2=")))
          new_var.infos["SBR2"] = std::string(info_key_value.begin() + 5, info_key_value.end());
        else if (same_prefix(info_key_value, std::string("SBR=")))
          new_var.infos["SBR"] = std::string(info_key_value.begin() + 4, info_key_value.end());
        else if (same_prefix(info_key_value, std::string("Unaligned=")))
          new_var.infos["Unaligned"] = std::string(info_key_value.begin() + 10, info_key_value.end());
      }
    }
  }

  // Parse samples, if any
  if (sample_names.size() > 0)
  {
    std::string const format = get_string_at_tab_index(line, tabs, 8);
    std::vector<std::size_t> all_format_colon = get_all_pos(format, ':');
    int32_t ad_field = -1;
    int32_t gt_field = -1;
    int32_t pl_field = -1;

    for (int32_t f = 0; f < static_cast<int32_t>(all_format_colon.size() - 1); ++f)
    {
      std::string const field = get_string_at_tab_index(format, all_format_colon, f);

      if (field == "AD")
        ad_field = f;
      else if (field == "GT")
        gt_field = f;
      else if (field == "PL")
        pl_field = f;
    }

    assert(ad_field != -1);
    assert(gt_field != -1);
    assert(pl_field != -1);
    std::size_t const FIELD_OFFSET = 9;

    for (std::size_t i = FIELD_OFFSET; i < sample_names.size() + FIELD_OFFSET; ++i)
    {
      // Create a new sample call
      SampleCall new_call;

      // Parse string of sample i
      std::string sample_string = get_string_at_tab_index(line, tabs, i);
      std::vector<std::size_t> sample_string_colon = get_all_pos(sample_string, ':');

      // Parse GT (for phase)
      std::string const gt_str = get_string_at_tab_index(sample_string, sample_string_colon, gt_field);

      // Check if the genotypes are phased
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

      // Parse AD
      std::string const ad_str = get_string_at_tab_index(sample_string, sample_string_colon, ad_field);
      std::vector<std::size_t> ad_str_comma = get_all_pos(ad_str, ',');

      for (std::size_t j = 0; j < ad_str_comma.size() - 1; ++j)
        new_call.coverage.push_back(static_cast<uint16_t>(std::stoul(get_string_at_tab_index(ad_str, ad_str_comma, j))));

      // Parse PL
      std::string const pl_str = get_string_at_tab_index(sample_string, sample_string_colon, pl_field);
      std::vector<std::size_t> pl_str_comma = get_all_pos(pl_str, ',');

      for (std::size_t j = 0; j < pl_str_comma.size() - 1; ++j)
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
  while (true)
  {
    std::string line;
    vcf_file >> line;

    if (line.size() == 0)
    {
      BOOST_LOG_TRIVIAL(error) << "[vcf] Could not find any line with samples in '" << filename << "'.";
      std::exit(1);
    }

    assert(line.size() > 2);

    if (line[0] == '#' && line[1] != '#')
    {
      // Gather list of all samples using this line
      assert(sample_names.size() == 0);
      std::vector<std::size_t> const all_line_tabs = get_all_pos(line, '\t');

      // Add all samples
      for (std::size_t i = 9; i < all_line_tabs.size() - 1; ++i)
        sample_names.push_back(get_string_at_tab_index(line, all_line_tabs, i));

      return;
    }
  }
}


void
Vcf::read()
{
  open_vcf_file_for_reading();
  read_samples();

  // Read all records
  while (read_record())
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
    vcf_file = ifopen(filename.c_str(), "w", InputFile::UNCOMPRESSED);
    break;

  case WRITE_BGZF_MODE:
    vcf_file = ifopen(filename.c_str(), "wb", InputFile::BGZF);
    break;

  default:
    std::cerr << "[vcf] ERROR: Trying to write in reading mode." << std::endl;
    std::exit(1);
  }
}


void
Vcf::write_header()
{
  assert(vcf_file);
  // Basic info
  *vcf_file << "##fileformat=VCFv4.2\n";
  *vcf_file << "##fileDate=" << current_date() << "\n";
  *vcf_file << "##source=Graphtyper\n";
  //*vcf_file << "##reference=\n"; // TODO: Add file location of the reference genome used.

  // Contigs
  {
    *vcf_file << "##contig=<ID=chr1,length=248956422>\n";
    *vcf_file << "##contig=<ID=chr2,length=242193529>\n";
    *vcf_file << "##contig=<ID=chr3,length=198295559>\n";
    *vcf_file << "##contig=<ID=chr4,length=190214555>\n";
    *vcf_file << "##contig=<ID=chr5,length=181538259>\n";
    *vcf_file << "##contig=<ID=chr6,length=170805979>\n";
    *vcf_file << "##contig=<ID=chr7,length=159345973>\n";
    *vcf_file << "##contig=<ID=chr8,length=145138636>\n";
    *vcf_file << "##contig=<ID=chr9,length=138394717>\n";
    *vcf_file << "##contig=<ID=chr10,length=133797422>\n";
    *vcf_file << "##contig=<ID=chr11,length=135086622>\n";
    *vcf_file << "##contig=<ID=chr12,length=133275309>\n";
    *vcf_file << "##contig=<ID=chr13,length=114364328>\n";
    *vcf_file << "##contig=<ID=chr14,length=107043718>\n";
    *vcf_file << "##contig=<ID=chr15,length=101991189>\n";
    *vcf_file << "##contig=<ID=chr16,length=90338345>\n";
    *vcf_file << "##contig=<ID=chr17,length=83257441>\n";
    *vcf_file << "##contig=<ID=chr18,length=80373285>\n";
    *vcf_file << "##contig=<ID=chr19,length=58617616>\n";
    *vcf_file << "##contig=<ID=chr20,length=64444167>\n";
    *vcf_file << "##contig=<ID=chr21,length=46709983>\n";
    *vcf_file << "##contig=<ID=chr22,length=50818468>\n";
    *vcf_file << "##contig=<ID=chrX,length=156040895>\n";
    *vcf_file << "##contig=<ID=chrY,length=57227415>\n";
    *vcf_file << "##contig=<ID=chrM>\n";
    *vcf_file << "##contig=<ID=chrUn>\n";
  }

  // INFO
  if (sample_names.size() > 0 && segments.size() == 0)
  {
    *vcf_file << "##INFO=<ID=ABHet,Number=1,Type=Float,Description=\"Allele Balance for heterozygous calls (call2/(call1+call2)) where are called genotype is call1/call2. -1 if not available.\">\n";
    *vcf_file << "##INFO=<ID=ABHom,Number=1,Type=Float,Description=\"Allele Balance for homozygous calls (A/(A+O)) where A is the called allele and O is anything else. -1 if not available.\">\n";
    *vcf_file << "##INFO=<ID=ABHetMulti,Number=A,Type=Float,Description=\"List of Allele Balance values for multiallelic heterozygous calls (alt/(ref+alt)). Each value corresponds to an alt in the same order as they appear. -1 if not available.\">\n";
    *vcf_file << "##INFO=<ID=ABHomMulti,Number=R,Type=Float,Description=\"List of Allele Balance values for multiallelic homozygous calls (A/(A+0)) where A is the called allele and O is anything else. Each value corresponds to a ref or alt in the same order as they appear. -1 if not available.\">\n";
    *vcf_file << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Number of alternate alleles in called genotypes.\">\n";
    *vcf_file << "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Number of alleles in called genotypes.\">\n";
    *vcf_file << "##INFO=<ID=CR,Number=1,Type=Integer,Description=\"Number of clipped reads.\">\n";
    *vcf_file << "##INFO=<ID=GX,Number=1,Type=Integer,Description=\"Graph complexity, 10*log10(#paths around the variant).\">\n";
    *vcf_file << "##INFO=<ID=MaxAAS,Number=A,Type=Integer,Description=\"Maximum alternative allele support per alt. allele.\">\n";
    *vcf_file << "##INFO=<ID=MaxAASR,Number=A,Type=Float,Description=\"Maximum alternative allele support ratio per alt. allele.\">\n";
    *vcf_file << "##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"Root-mean-square mapping quality.\">\n";
    *vcf_file << "##INFO=<ID=MQ0,Number=1,Type=Integer,Description=\"Number of reads with MQ=0.\">\n";
    *vcf_file << "##INFO=<ID=MQperAllele,Number=R,Type=Integer,Description=\"Mapping quality of reads aligned to each allele.\">\n";
    *vcf_file << "##INFO=<ID=QD,Number=1,Type=Float,Description=\"QUAL divided by NonReferenceSeqDepth.\">\n";
    *vcf_file << "##INFO=<ID=PS,Number=1,Type=Integer,Description=\"Unique ID of the phase set this variant is a member of. "
              << "If the calls are unphased, it is used to represent which haplotype set the variant is a member of.\">\n";
    *vcf_file << "##INFO=<ID=RACount,Number=R,Type=Integer,Description=\"Number of realigned reads which were uniquely aligned to this allele.\">\n";
    *vcf_file << "##INFO=<ID=RADist,Number=R,Type=Integer,Description=\"Total realignment distance of all reads per allele."
              << " Realignment distance = abs((Original aligned begin position) - (New aligned begin position)).\">\n";
    *vcf_file << "##INFO=<ID=RAMeanDist,Number=R,Type=Integer,Description=\"Mean realignment distance per allele.\">\n";
    *vcf_file << "##INFO=<ID=RefLen,Number=1,Type=Integer,Description=\"Length of the reference allele.\">\n";
    *vcf_file << "##INFO=<ID=SB,Number=1,Type=Float,Description=\"Strand bias (F/(F+R)) where F and R are forward and reverse strands, respectively. -1 if not available.\">\n";
    *vcf_file << "##INFO=<ID=SBF,Number=R,Type=Integer,Description=\"Number of forward stranded reads per allele.\">\n";
    *vcf_file << "##INFO=<ID=SBF1,Number=R,Type=Integer,Description=\"Number of first forward stranded reads per allele.\">\n";
    *vcf_file << "##INFO=<ID=SBF2,Number=R,Type=Integer,Description=\"Number of second forward stranded reads per allele.\">\n";
    *vcf_file << "##INFO=<ID=SBR,Number=R,Type=Integer,Description=\"Number of reverse stranded reads per allele.\">\n";
    *vcf_file << "##INFO=<ID=SBR1,Number=R,Type=Integer,Description=\"Number of first reverse stranded reads per allele.\">\n";
    *vcf_file << "##INFO=<ID=SBR2,Number=R,Type=Integer,Description=\"Number of second reverse stranded reads per allele.\">\n";
    *vcf_file << "##INFO=<ID=SeqDepth,Number=1,Type=Integer,Description=\"Total accumulated sequencing depth over all the samples.\">\n";
    *vcf_file << "##INFO=<ID=Unaligned,Number=1,Type=Integer,Description=\"Number of previously unaligned reads that got re-aligned.\">\n";
    *vcf_file << "##INFO=<ID=VarType,Number=1,Type=String,Description=\"First letter is program identifier, the second letter is variant type.\">\n";
  }
  else
  {
    *vcf_file << "##INFO=<ID=RefLen,Number=1,Type=Integer,Description=\"Length of the reference allele.\">\n";
  }

  // FORMAT
  if (sample_names.size() > 0 && segments.size() == 0)
  {
    *vcf_file << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    *vcf_file << "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n";
    *vcf_file << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth\">\n";
    *vcf_file << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality.\">\n";
    *vcf_file << "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"PHRED-scaled genotype likelihoods\">\n";
  }
  else if (segments.size() > 0)
  {
    *vcf_file << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    *vcf_file << "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"PHRED-scaled genotype likelihoods\">\n";
  }

  *vcf_file << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";

  if (sample_names.size() > 0)
    *vcf_file << "\tFORMAT";

  for (auto const & sample_name : sample_names)
    *vcf_file << "\t" << sample_name;

  *vcf_file << "\n";
}


void
Vcf::write_record(Variant const & var, std::string const & suffix, const bool FILTER_ZERO_QUAL)
{
  // Parse the position
  auto contig_pos = absolute_pos.get_contig_position(var.abs_pos);

  if (!Options::instance()->output_all_variants && var.calls.size() > 0 && var.seqs.size() > 85)
  {
    BOOST_LOG_TRIVIAL(info) << "[vcf] Skipped outputting variant at position " << contig_pos.first << ":" << contig_pos.second
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
        BOOST_LOG_TRIVIAL(info) << "[vcf] Skipped outputting variant at position " << contig_pos.first << ":" << contig_pos.second
                                << " because the total length of alleles is too high.";
        return;
      }
    }
  }

  // Calculate qual
  const uint64_t variant_qual = var.get_qual();

  if (FILTER_ZERO_QUAL && variant_qual == 0)
    return;

  *vcf_file << contig_pos.first << "\t" << contig_pos.second;

  // Parse the ID
  *vcf_file << "\t" << contig_pos.first << ":" << contig_pos.second << ":" << var.determine_variant_type() << suffix;

  // Parse the sequences
  assert(var.seqs.size() >= 2);
  *vcf_file << "\t" << std::string(var.seqs[0].begin(), var.seqs[0].end()) << "\t" << std::string(var.seqs[1].begin(), var.seqs[1].end());

  for (auto it = var.seqs.begin() + 2; it != var.seqs.end(); ++it)
    *vcf_file << "," << std::string(it->begin(), it->end());

  // Parse qual
  *vcf_file << "\t" << std::to_string(variant_qual);

  // Parse filter
  *vcf_file << "\t.";

  // Parse info
  if (var.infos.empty())
  {
    *vcf_file << "\t.";
  }
  else if (var.calls.size() > 0)
  {
    *vcf_file << "\t" << var.infos.begin()->first << "=" << var.infos.begin()->second;

    for (auto map_it = std::next(var.infos.begin(), 1); map_it != var.infos.end(); ++map_it)
      *vcf_file << ";" << map_it->first << "=" << map_it->second;
  }
  else
  {
    if (var.infos.count("RefLen") == 1)
      *vcf_file << "\tRefLen=" << var.infos.at("RefLen");
    else
      *vcf_file << "\t.";
  }

  if (sample_names.size() != var.calls.size())
  {
    std::cerr << "[vcf] ERROR: Sample name and variant calls sizes don't match: " << sample_names.size() << " and " << var.calls.size() << std::endl;
    std::exit(1);
  }

  // Parse FORMAT
  if (sample_names.size() > 0)
    *vcf_file << "\tGT:AD:DP:GQ:PL";

  for (std::size_t i = 0; i < var.calls.size(); ++i)
  {
    auto const & call = var.calls[i];

    // Write GT
    {
      // If all PHRED scores are zero, print ./. or .|.
      if (std::find_if(call.phred.begin(), call.phred.end(), [](uint8_t const pl){return pl != 0;}) == call.phred.end())
      {
        if (var.phase.size() > 0)
          *vcf_file << "\t.|.";
        else
          *vcf_file << "\t./.";
      }
      else
      {
        std::pair<uint16_t, uint16_t> const gt_call = call.get_gt_call();

        if (var.phase.size() > 0)
        {
          assert(i < var.phase.size());

          if (var.phase[i] == 0)
            *vcf_file << "\t" << gt_call.first << "|" << gt_call.second;
          else
            *vcf_file << "\t" << gt_call.second << "|" << gt_call.first;
        }
        else
        {
          *vcf_file << "\t" << gt_call.first << "/" << gt_call.second;
        }
      }
    }

    // Write AD:DP
    {
      assert(call.coverage.size() > 0);
      *vcf_file << ":" << static_cast<uint16_t>(call.coverage[0]);

      for (auto ad_it = call.coverage.begin() + 1; ad_it != call.coverage.end(); ++ad_it)
        *vcf_file << "," << static_cast<uint16_t>(*ad_it);
    }

    // Write DP
    {
      *vcf_file << ":" << std::accumulate(call.coverage.begin(), call.coverage.end(), static_cast<uint32_t>(0u));
    }

    // Write GQ
    {
      *vcf_file << ":" << static_cast<uint16_t>(call.get_gq());
    }

    // Write PL
    {
      *vcf_file << ":" << static_cast<uint16_t>(call.phred[0]); // This cast to uint16_t is needed! Otherwise it is converted to char!

      for (unsigned i = 1; i < call.phred.size(); ++i)
        *vcf_file << "," << static_cast<uint16_t>(call.phred[i]); // ..it is also needed here.
    }
  }

  // Fin.
  *vcf_file << "\n";
}


void
Vcf::write(std::string const & region)
{
  open_for_writing();
  write_header();
  write_records(region);
  close_vcf_file();
}


void
Vcf::write_records(std::string const & region, const bool FILTER_ZERO_QUAL)
{
  if (variants.size() == 0)
    return;

  // Sort the variants
  auto compare_variants = [&](std::size_t const i, std::size_t const j) -> bool
                          {
                            assert(i < variants.size());
                            assert(j < variants.size());
                            gyper::Variant const & a = variants[i];
                            gyper::Variant const & b = variants[j];
                            assert(a.seqs.size() >= 2);
                            assert(b.seqs.size() >= 2);
                            return a.abs_pos < b.abs_pos ||
                              (a.abs_pos == b.abs_pos && a.determine_variant_type() < b.determine_variant_type());
                          };

  auto same_variant_id = [](Variant const & a, Variant const & b) -> bool
  {
    return a.abs_pos == b.abs_pos && a.determine_variant_type() == b.determine_variant_type();
  };

  std::vector<std::size_t> indexes;

  for (std::size_t i = 0; i < variants.size(); ++i)
    indexes.push_back(i);

  GenomicRegion genomic_region(region);
  uint32_t const region_begin = 1 + absolute_pos.get_absolute_position(genomic_region.chr, genomic_region.begin);
  uint32_t const region_end = absolute_pos.get_absolute_position(genomic_region.chr, genomic_region.end);
  std::sort(indexes.begin(), indexes.end(), compare_variants);
  assert(indexes.size() > 0);

  auto inside_region = [&](uint32_t const pos) -> bool
  {
    return pos >= region_begin && pos <= region_end;
  };

  std::array<std::string, 4u> latin_suffix = {".bis", ".ter", ".quater", ".quinquies"};
  int32_t dup = -1; // -1 means no duplication

  if (inside_region(variants[indexes[0]].abs_pos))
    write_record(variants[indexes[0]], "" /*suffix*/, FILTER_ZERO_QUAL);

  // Write the variants in the correct order
  for (std::size_t i = 1; i < indexes.size(); ++i)
  {
    // Make sure the variant is unique
    if (inside_region(variants[indexes[i]].abs_pos))
    {
      if (!same_variant_id(variants[indexes[i]], variants[indexes[i - 1]]))
      {
        write_record(variants[indexes[i]], "" /*suffix*/, FILTER_ZERO_QUAL);
        dup = -1;
      }
      else
      {
        ++dup;
        assert(dup >= 0);

        if (dup < static_cast<int8_t>(latin_suffix.size()))
          write_record(variants[indexes[i]], latin_suffix[dup], FILTER_ZERO_QUAL);
        else
          write_record(variants[indexes[i]], std::string(".") + std::to_string(dup), FILTER_ZERO_QUAL);
      }
    }
  }
}


void
Vcf::write_segments()
{
  BOOST_LOG_TRIVIAL(info) << "[graphtyper::vcf] Writing " << segments.size() << " segments to " << filename;

  for (auto const & segment : segments)
  {
    assert(sample_names.size() == segment.segment_calls.size());
    auto contig_pos = absolute_pos.get_contig_position(segment.id);

    // Parse CHROM and POS
    *vcf_file << contig_pos.first << "\t" << contig_pos.second;

    // Parse the ID
    *vcf_file << "\t" << contig_pos.first << ":" << contig_pos.second << ":" << segment.var_type;

    if (segment.extra_id != -1)
      *vcf_file << ":" << std::to_string(segment.extra_id);

    // Parse the sequences
    *vcf_file << "\t" << segment.get_ref_string() << "\t" << segment.get_alt_string();

    // Parse QUAL
    *vcf_file << "\t" << ".";

    // Parse FILTER
    *vcf_file << "\t" << ".";

    // Parse INFO
    *vcf_file << "\t" << "RefLen=" << std::to_string(segment.ref_size);

    // Parse FORMAT
    *vcf_file << "\t" << "GT:PL";

    for (auto const & segment_call : segment.segment_calls)
    {
      *vcf_file << "\t" << segment_call.call.first << "/" << segment_call.call.second << ":";

      assert(segment_call.phred.size() > 1);
      *vcf_file << std::to_string(static_cast<uint16_t>(segment_call.phred[0]));

      for (unsigned p = 1; p < segment_call.phred.size(); ++p)
        *vcf_file << "," << std::to_string(static_cast<uint16_t>(segment_call.phred[p]));
    }

    *vcf_file << "\n";
  }
}


void
Vcf::close_vcf_file()
{
  ifclose(vcf_file);
  delete vcf_file;
  vcf_file = nullptr;
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
  bool constexpr FORCE_BIALLELIC = true;

  if (!FORCE_BIALLELIC || segment.allele_names.size() == 2 /*check if it is already biallelic*/)
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
  // Add each genotype individually
  std::vector<Variant> new_vars;
  assert(haplotype.gts.size() > 0);

  for (auto const & gt : haplotype.gts)
    new_vars.push_back(Variant(gt));

  assert(new_vars.size() == haplotype.gts.size());
  assert(new_vars.size() == haplotype.var_stats.size());

  // Add variant stats
  for (std::size_t i = 0; i < new_vars.size(); ++i)
  {
    new_vars[i].infos["CR"] = std::to_string(haplotype.var_stats[i].clipped_reads);
    new_vars[i].infos["GX"] = std::to_string(haplotype.var_stats[i].graph_complexity);
    new_vars[i].infos["MQ"] = std::to_string(haplotype.var_stats[i].get_rms_mapq());
    new_vars[i].infos["MQ0"] = std::to_string(haplotype.var_stats[i].mapq_zero_count);
    new_vars[i].infos["MQperAllele"] = haplotype.var_stats[i].get_rms_mapq_per_allele();
    new_vars[i].infos["PS"] = std::to_string(phase_set);
    new_vars[i].infos["RACount"] = haplotype.var_stats[i].get_realignment_count();
    new_vars[i].infos["RADist"] = haplotype.var_stats[i].get_realignment_distance();
    new_vars[i].infos["SBF"] = haplotype.var_stats[i].get_forward_strand_bias();
    new_vars[i].infos["SBR"] = haplotype.var_stats[i].get_reverse_strand_bias();
    new_vars[i].infos["SBF1"] = haplotype.var_stats[i].get_r1_forward_strand_bias();
    new_vars[i].infos["SBF2"] = haplotype.var_stats[i].get_r2_forward_strand_bias();
    new_vars[i].infos["SBR1"] = haplotype.var_stats[i].get_r1_reverse_strand_bias();
    new_vars[i].infos["SBR2"] = haplotype.var_stats[i].get_r2_reverse_strand_bias();
    new_vars[i].infos["Unaligned"] = haplotype.var_stats[i].get_unaligned_count();
  }

  for (auto const & sample : haplotype.hap_samples)
  {
    // Calculate the haplotypic phred scores
    std::vector<uint8_t> phase; // 0 = same as unphased, 1 other way around
    std::vector<std::vector<uint8_t> > gt_phred = get_genotype_phred(sample, phase, haplotype.gts);
    assert(phase.size() == new_vars.size());

    for (std::size_t i = 0; i < new_vars.size(); ++i)
    {
      assert(i < gt_phred.size());
      assert(i < sample.gt_coverage.size());
      new_vars[i].calls.push_back(SampleCall(gt_phred[i], sample.gt_coverage[i])); // Add new call

      if (Options::instance()->phased_output)
      {
        new_vars[i].phase.push_back(phase[i]);
        assert(new_vars[i].calls.size() == new_vars[i].phase.size());
      }

      assert(new_vars[i].seqs.size() == sample.gt_coverage[i].size());
    }
  }

  // Add the variants
  std::move(new_vars.begin(), new_vars.end(), std::back_inserter(variants));

  if (clear_haplotypes)
    haplotype.clear();
}


void
Vcf::add_haplotypes_for_extraction(std::vector<std::vector<Genotype> > const & gts, std::vector<std::vector<uint32_t> > const & hap_calls)
{
  assert(hap_calls.size() == gts.size());
  std::vector<Variant> new_vars;

  for (std::size_t i = 0; i < gts.size(); ++i)
  {
    // Only add variants if there is something else than the reference called
    if (hap_calls[i].size() > 1)
    {
      assert(hap_calls[i][0] == 0);
      new_vars.push_back(Variant(gts[i], hap_calls[i]));
    }
    else
    {
      assert(hap_calls[i].size() == 1);
      assert(hap_calls[i][0] == 0);
    }
  }

  bool const BREAK_DOWN_EXTRACTED_HAPLOTYPES = !Options::instance()->skip_breaking_down_extracted_haplotypes;

  for (auto && var : new_vars)
  {
    // Only try to split multiallelic variants
    if (BREAK_DOWN_EXTRACTED_HAPLOTYPES)
    {
      std::vector<Variant> new_broken_down_vars = break_down_variant(std::move(var), SPLIT_VAR_THRESHOLD);
      std::move(new_broken_down_vars.begin(), new_broken_down_vars.end(), std::back_inserter(this->variants));
    }
    else
    {
      this->variants.push_back(std::move(var));
    }
  }
}


void
Vcf::post_process_variants(bool const NORMALIZE)
{
  if (NORMALIZE)
  {
    for (auto & var : variants)
      var.normalize();
  }
  else
  {
    for (auto & var : variants)
      var.trim_sequences(false); // Don't keep one match
  }

  // Generate the INFO field if there are any samples
  if (sample_names.size() > 0)
  {
    for (auto & var : variants)
      var.generate_infos();
  }
}


void
Vcf::add_read_group_to_pn(std::unordered_map<std::string, std::string> const & new_rg2pn)
{
  for (auto map_it = new_rg2pn.begin(); map_it != new_rg2pn.end(); ++map_it)
  {
    assert(rg2pn.count(map_it->first) == 0);
    rg2pn[map_it->first] = map_it->second;
  }
}


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
get_string_at_tab_index(std::string const & line, std::vector<std::size_t> const & tabs, std::size_t const index)
{
  assert(index + 1 < tabs.size());
  assert(tabs[index + 1] > tabs[index]);

  if (tabs[index + 1] == std::string::npos)
    return line.substr(tabs[index], std::string::npos);
  else
    return line.substr(tabs[index], tabs[index + 1] - tabs[index] - 1);
}


} // namespace gyper
