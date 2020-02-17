#pragma once

#include <string>
#include <sstream>
#include <unordered_map>
#include <vector>

//#include <seqan/basic.h>
//#include <seqan/sequence.h>


namespace gyper
{

void
get_sample_name_from_bam_header(std::string const & hts_filename,
                                std::vector<std::string> & samples,
                                std::unordered_map<std::string, int> & rg2sample_i);

//std::vector<std::pair<seqan::CharString, seqan::Dna5String> > read_fasta_sequences(std::string const & fasta_filename);
//std::map<std::string, std::vector<seqan::Dna5String> > read_haplotypes_from_fasta(std::string const & fasta_filename);
void append_to_file(std::string && data, std::string const & file_name);
void write_to_file(std::string && data, std::string const & file_name);
void write_gzipped_to_file(std::stringstream & ss, std::string const & file_name, bool const append = false);
std::unordered_map<std::string, long> get_contig_to_lengths(std::string const & fai_filename);

} // namespace gyper
