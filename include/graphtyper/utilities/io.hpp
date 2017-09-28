#pragma once

#include <string>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include <seqan/basic.h>
#include <seqan/sequence.h>


namespace gyper
{

std::vector<std::string>
get_sample_names_from_bam_header(std::string const & hts_filename,
                                 std::unordered_map<std::string, std::string> & rg2pn
                                 );

std::vector<std::pair<seqan::CharString, seqan::Dna5String> > read_fasta_sequences(std::string const & fasta_filename);
std::map<std::string, std::vector<seqan::Dna5String> > read_haplotypes_from_fasta(std::string const & fasta_filename);
void append_to_file(std::string && data, std::string const & file_name);
void write_to_file(std::string && data, std::string const & file_name);
void write_gzipped_to_file(std::stringstream & ss, std::string const & file_name);

} // namespace gyper
