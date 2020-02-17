#include <iostream>
#include <sstream>
#include <fstream>

#include <boost/algorithm/string.hpp>
#include <boost/log/trivial.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/bam_io.h>

#include <graphtyper/constants.hpp>
#include <graphtyper/utilities/io.hpp>


namespace gyper
{

void
get_sample_name_from_bam_header(std::string const & hts_filename,
                                std::vector<std::string> & samples,
                                std::unordered_map<std::string, int> & rg2sample_i)
{
  seqan::HtsFileIn hts_file;

  if (!open(hts_file, hts_filename.c_str()))
  {
    BOOST_LOG_TRIVIAL(error) << "[graphtyper::utilities::io] Could not open " << hts_filename << " for reading";
    std::exit(1);
  }

  std::string const header_text(hts_file.hdr->text, hts_file.hdr->l_text);
  std::vector<std::string> header_lines;

  // Split the header text into lines
  boost::split(header_lines, header_text, boost::is_any_of("\n"));

  for (auto & line_it : header_lines)
  {
    if (boost::starts_with(line_it, "@RG"))
    {
      std::size_t const pos_id = line_it.find("\tID:");
      std::size_t const pos_samp = line_it.rfind("\tSM:");

      if (pos_samp == std::string::npos || pos_id == std::string::npos)
      {
        std::cerr << "[graphtyper::utilities::io] ERROR: Could not parse RG and sample from header line:"
                  << line_it << std::endl;
        std::exit(1);
      }

      std::size_t pos_id_ends = line_it.find("\t", pos_id + 1);

      // Check if this is the last field
      if (pos_id_ends == std::string::npos)
        pos_id_ends = line_it.size();

      std::size_t pos_samp_ends = line_it.find("\t", pos_samp + 1);

      // Check if this is the last field
      if (pos_samp_ends == std::string::npos)
        pos_samp_ends = line_it.size();

      std::string new_id = line_it.substr(pos_id + 4, pos_id_ends - pos_id - 4);
      std::string new_sample = line_it.substr(pos_samp + 4, pos_samp_ends - pos_samp - 4);

#ifndef NDEBUG
      BOOST_LOG_TRIVIAL(debug) << "[graphtyper::utilities::io] Added RG: '"
                               << new_id << "' => '" << new_sample << "'";
#endif // NDEBUG

      auto find_it = std::find(samples.begin(), samples.end(), new_sample);

      // check if this is a new sample
      if (find_it == samples.end())
      {
        rg2sample_i[new_id] = samples.size();
        samples.push_back(new_sample);
      }
      else
      {
        rg2sample_i[new_id] = std::distance(samples.begin(), find_it);
      }
    }
  }
}


std::unordered_map<std::string, long>
get_contig_to_lengths(std::string const & fai_filename)
{
  std::ifstream fai(fai_filename.c_str());

  if (!fai.is_open())
  {
    BOOST_LOG_TRIVIAL(error) << "Could not open fasta.fai at: " << fai_filename;
    std::exit(1);
  }

  std::unordered_map<std::string, long> contig_to_lengths;
  std::string line;

  while (std::getline(fai, line))
  {
    if (line.size() == 0)
      continue;

    auto find_it = std::find(line.begin(), line.end(), '\t');

    if (find_it == line.end())
    {
      BOOST_LOG_TRIVIAL(warning) << "[graphtyper::io] Could not parse contig line: " << line;
      continue;
    }

    std::string contig(line.begin(), find_it);
    auto find2_it = std::find(find_it + 1, line.end(), '\t');

    if (find2_it == line.end())
    {
      BOOST_LOG_TRIVIAL(warning) << "[graphtyper::io] Could not parse contig line: " << line;
      continue;
    }

    std::string length(find_it + 1, find2_it);

    if (contig_to_lengths.count(contig) == 1)
      BOOST_LOG_TRIVIAL(warning) << "[graphtyper::io] Duplicated contig in FAI: " << contig;

    contig_to_lengths[contig] = std::stoull(length);
  }

  return contig_to_lengths;
}


/*
std::vector<std::pair<seqan::CharString, seqan::Dna5String> >
read_fasta_sequences(std::string const & fasta_filename)
{
  std::vector<std::pair<seqan::CharString, seqan::Dna5String> > fasta_sequences;

  seqan::StringSet<seqan::CharString> ids;
  seqan::StringSet<seqan::Dna5String> seqs;

  seqan::SeqFileIn fasta_file(fasta_filename.c_str());
  seqan::readRecords(ids, seqs, fasta_file);

  fasta_sequences.reserve(seqan::length(ids));

  auto id_it = seqan::begin(ids);
  auto seq_it = seqan::begin(seqs);

  while (seq_it != seqan::end(seqs))
  {
    fasta_sequences.push_back({*id_it, *seq_it});

    ++seq_it;
    ++id_it;
  }

  return fasta_sequences;
}


std::map<std::string, std::vector<seqan::Dna5String> >
read_haplotypes_from_fasta(std::string const & fasta_filename)
{
  std::map<std::string, std::vector<seqan::Dna5String> > haplotypes;
  seqan::SeqFileIn fasta_file(fasta_filename.c_str());
  seqan::StringSet<seqan::CharString> ids;
  seqan::StringSet<seqan::Dna5String> seqs;
  seqan::readRecords(ids, seqs, fasta_file);

  for (unsigned i = 0; i < seqan::length(ids); ++i)
  {
    seqan::CharString allele;
    bool found_star = false;

    for (unsigned j = 0; j < seqan::length(ids[i]); ++j)
    {
      char const & c = ids[i][j];

      if (c == '*')
      {
        found_star = true;
        seqan::appendValue(allele, c);
      }
      else if (c == '_')
      {
        break;
      }
      else if (c == ' ')
      {
        if (found_star)
          break;
        else
          seqan::clear(allele);
      }
      else
      {
        seqan::appendValue(allele, c);
      }
    }

    assert(seqan::length(allele) > 1);
    std::string allele_str(seqan::toCString(allele));

    // Add "HLA-" in front if it is missing
    if (allele_str[1] == '*')
    {
      allele_str = std::string("HLA-").append(allele_str);
    }

    auto found_it = haplotypes.find(allele_str);

    if (found_it == haplotypes.end())
    {
      haplotypes[allele_str] = {seqs[i]};
    }
    else
    {
      found_it->second.push_back(seqs[i]);
    }
  }

  assert(haplotypes.size() > 0);
  return haplotypes;
}
*/


void
append_to_file(std::string && data, std::string const & file_name)
{
  if (file_name == "-")
  {
    std::cout << data;
  }
  else
  {
    std::ofstream myfile;
    myfile.open(file_name.c_str(), std::ios_base::app);

    if (!myfile.is_open())
    {
      BOOST_LOG_TRIVIAL(error) << "[graphtyper::io] Cannot write to " << file_name;
      std::exit(1);
    }

    myfile << data;
    myfile.flush();
    myfile.close();
  }
}


void
write_to_file(std::string && data, std::string const & file_name)
{
  if (file_name == "-")
  {
    std::cout << data;
  }
  else
  {
    std::ofstream myfile;
    myfile.open(file_name.c_str(), std::ios_base::trunc);

    if (!myfile.is_open())
    {
      BOOST_LOG_TRIVIAL(error) << "[graphtyper::io] Cannot write to " << file_name;
      std::exit(1);
    }

    myfile << data;
    myfile.flush();
    myfile.close();
  }
}


void
write_gzipped_to_file(std::stringstream & ss, std::string const & file_name, bool const append)
{
  std::ofstream compressed(file_name.c_str(),
                           append ? (std::ofstream::binary | std::ofstream::app) : std::ofstream::binary
                           );

  if (!compressed.is_open())
  {
    BOOST_LOG_TRIVIAL(error) << "[graphtyper::io] Could not open file '" << file_name << "'";
    std::exit(3);
  }

  boost::iostreams::filtering_streambuf<boost::iostreams::input> out;
  out.push(boost::iostreams::gzip_compressor());
  out.push(ss);
  boost::iostreams::copy(out, compressed);
}


} // namespace gyper
