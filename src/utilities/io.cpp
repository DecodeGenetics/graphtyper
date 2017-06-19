#include <iostream>
#include <sstream>
#include <fstream>

#include <boost/algorithm/string.hpp>
#include <boost/log/trivial.hpp>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/bam_io.h>

#include <graphtyper/constants.hpp>
#include <graphtyper/utilities/io.hpp>


namespace gyper
{

std::vector<std::string>
get_sample_names_from_bam_header(std::string const & hts_filename, std::unordered_map<std::string, std::string> & rg2sample)
{
  std::vector<std::string> samples;
  seqan::HtsFileIn hts_file;

  if (!open(hts_file, hts_filename.c_str()))
  {
    BOOST_LOG_TRIVIAL(error) << "Could not open " << hts_filename << " for reading";
    return samples;
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
      std::size_t const pos_samp = line_it.rfind("\tSM:"); // use rfind for superoptimization!!!!1one

      if (pos_samp != std::string::npos and pos_id != std::string::npos)
      {
        std::size_t pos_id_ends = line_it.find("\t", pos_id + 1);

        // Check if this is the last field
        if (pos_id_ends == std::string::npos)
          pos_id_ends = line_it.size();

        std::size_t pos_samp_ends = line_it.find("\t", pos_samp + 1);

        // Check if this is the last field
        if (pos_samp_ends == std::string::npos)
          pos_samp_ends = line_it.size();

        std::string new_id = line_it.substr(pos_id + 4, pos_id_ends - pos_id - 4);
        assert(pos_samp + 11 < line_it.size());

        std::string new_sample = line_it.substr(pos_samp + 4, pos_samp_ends - pos_samp - 4);

        // Make sure this read group has not been seen
        assert(rg2sample.count(new_id) == 0);
        rg2sample[new_id] = new_sample;
        BOOST_LOG_TRIVIAL(info) << "[graphtyper::io] Added RG: '" << new_id << "' => '" << new_sample << "'";

        auto find_it = std::find(samples.begin(), samples.end(), new_sample);

        if (find_it == samples.end())
          samples.push_back(new_sample);
      }
    }
  }

  std::sort(samples.begin(), samples.end());
  return samples;
}


std::vector<std::pair<seqan::CharString, seqan::Dna5String> >
read_fasta_sequences(std::string fasta_filename)
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
read_haplotypes_from_fasta(std::string fasta_filename)
{
  std::map<std::string, std::vector<seqan::Dna5String> > haplotypes;
  seqan::SeqFileIn fasta_file(fasta_filename.c_str());
  seqan::StringSet<seqan::CharString> ids;
  seqan::StringSet<seqan::Dna5String> seqs;
  seqan::readRecords(ids, seqs, fasta_file);

  for (unsigned i = 0; i < seqan::length(ids); ++i)
  {
    // std::cout << "Id, seq = " << ids[i] << '\n' << seqan::length(seqs[i]) << std::endl;
    // if (seqan::length(seqs[i]) < n * K - (n - 1))
    // {
    //   std::size_t const bases_needed = n * K / 2 - (n - 1) / 2 - seqan::length(seqs[i]) / 2;
    //
    //   // std::cout << "Sequence is short: " << seqs[i] << std::endl;
    //   seqan::Dna5String prev;
    //   seqan::Dna5String next;
    //
    //   if (i > 0 && seqan::back(ids[i]) != 'P') // Dont check prev at the 5P feature
    //   {
    //     prev = seqs[i - 1];
    //     std::cout << "prev is " << prev << std::endl;
    //   }
    //
    //   if (i < seqan::length(ids) - 1 && seqan::back(ids[i]) != '3') // Dont check next at the P3 feature
    //   {
    //     next = seqs[i + 1];
    //     std::cout << "next is " << next << std::endl;
    //   }
    //
    //   if (bases_needed < seqan::length(prev))
    //     seqan::erase(prev, 0, seqan::length(prev) - bases_needed);
    //
    //   if (bases_needed < seqan::length(next))
    //     seqan::resize(next, bases_needed);
    //
    //   seqan::append(prev, seqs[i]);
    //   seqan::append(prev, next);
    //
    //   if (seqan::length(prev) > (n * K - (n - 1)))
    //     seqan::resize(prev, (n * K - (n - 1)));
    //
    //   seqs[i] = prev;
    //   // std::cout << "Sequence is now longer " << seqs[i] << std::endl;
    // }

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

  assert (haplotypes.size() > 0);

  for (unsigned i = 0; i < haplotypes.begin()->second.size(); ++i)
  {
    auto find_it = std::find_if(haplotypes.begin(),
                                haplotypes.end(),
                                [i]
                                (std::pair<std::string, std::vector<seqan::Dna5String> > const & p)
                                {
                                  return seqan::length(p.second[i]) < 50;
                                });

    if (find_it != haplotypes.end())
    {
      std::cout << "I found a short sequence" << std::endl;
      std::cout << "!! hap first i = " << i << " " << seqan::length(find_it->second[i]) << ", " << find_it->first << std::endl;
    }
  }


  return haplotypes;
}


void
append_to_file(std::string && data, std::string file_name)
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
      std::cerr << "[io] ERROR: Cannot write to " << file_name << std::endl;
      myfile.close();
      std::exit(1);
    }

    myfile << data;
    myfile.flush();
    myfile.close();
  }
}


void
write_to_file(std::string && data, std::string file_name)
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
      std::cerr << "[io] ERROR: Cannot write to " << file_name << std::endl;
      myfile.close();
      std::exit(1);
    }

    myfile << data;
    myfile.flush();
    myfile.close();
  }
}


} // namespace gyper
