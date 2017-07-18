#pragma once

#include <mutex> // std::mutex
#include <string> // std::string
#include <unordered_map> // std::unordered_map
#include <vector> // std::vector

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/genotype.hpp>
#include <graphtyper/graph/haplotype.hpp>
#include <graphtyper/typer/segment.hpp>
#include <graphtyper/typer/variant.hpp>

class InputFile; // Implemented in libStatGen

namespace gyper
{

enum VCF_FILE_MODE
{
  READ_UNCOMPRESSED_MODE,
  READ_BGZF_MODE,
  WRITE_UNCOMPRESSED_MODE,
  WRITE_BGZF_MODE,
  READ_MODE, // Selects whether to use uncompressed or bgzf based on file extension
  WRITE_MODE
};


class Vcf
{
public:
  Vcf(VCF_FILE_MODE const _filemode= READ_UNCOMPRESSED_MODE, std::string const & filename = "hi_i_am_a.vcf");
  void open(VCF_FILE_MODE const _filemode, std::string const & _filename);
  void set_filemode(VCF_FILE_MODE const _filemode);

  /***************
   * CLASS ACESS *
   ***************/
  std::size_t get_sample_index(std::string const & sample_name) const;

  /******************
   * CLASS MODIFERS *
   ******************/
  /** I/O member functions */
  InputFile * vcf_file = nullptr;
  void open_vcf_file_for_reading();
  void read_samples();
  std::string read_line();
  bool read_record();
  void read(); /** \brief Reads the VCF file. */
  void open_for_writing();
  void write_header();
  void write_record(Variant const & var, std::string const & suffix = "", const bool FILTER_ZERO_QUAL = false);
  void write_segments();
  void write(std::string const & region = "."); /** \brief Writes the VCF file. */
  void write_records(std::string const & region = ".", const bool FILTER_ZERO_QUAL = false);
  void close_vcf_file();

  /** Common */
  void clear();

  // Adding data
  void add_segment(Segment && segment);
  void add_haplotype(Haplotype & haplotype, bool const clear_haplotypes, uint32_t const phase_set = 0);
  void add_haplotypes_for_extraction(std::vector<std::vector<Genotype> > const & gts, std::vector<std::vector<uint32_t> > const & hap_calls);
  void add_read_group_to_pn(std::unordered_map<std::string, std::string> const & new_rg2pn);

  // Modify data
  void post_process_variants(bool const NORMALIZE = true, bool const TRIM_SEQUENCES = true);

  VCF_FILE_MODE filemode;
  std::string filename;
  std::vector<std::string> sample_names;
  std::vector<Variant> variants;
  std::vector<Segment> segments;
  std::unordered_map<std::string, std::string> rg2pn;


  /************************
   * CLASS PRIVATE ACCESS *
   ************************/
private:
  std::mutex mutable vcf_mutex;
  std::size_t get_or_add_sample_index(std::string const & sample_name);
};

std::vector<std::size_t> get_all_pos(std::string const & line, char const delim = '\t');
std::string get_string_at_tab_index(std::string const & line, std::vector<std::size_t> const & tabs, std::size_t const index);

} // namespace gyper
