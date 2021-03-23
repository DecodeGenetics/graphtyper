#pragma once

#include <memory>
#include <string> // std::string
#include <unordered_map> // std::unordered_map
#include <unordered_set> // std::unordered_map
#include <vector> // std::vector

#include <boost/serialization/access.hpp>

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/genotype.hpp>
#include <graphtyper/graph/haplotype.hpp>
#include <graphtyper/typer/segment.hpp>
#include <graphtyper/typer/variant.hpp>
#include <graphtyper/utilities/bgzf_stream.hpp>
#include <graphtyper/utilities/gzstream.hpp>


namespace gyper
{

class HaplotypeCall;

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
  friend class boost::serialization::access;

public:
  Vcf();
  explicit Vcf(VCF_FILE_MODE _filemode, std::string const & filename);
  Vcf(Vcf const &) = delete;
  Vcf(Vcf &&) = delete;

  void open(VCF_FILE_MODE _filemode, std::string const & _filename);
  void set_filemode(VCF_FILE_MODE _filemode);

  /******************
   * CLASS MODIFERS *
   ******************/
  /** I/O member functions */
  igzstream bgzf_in;
  BGZF_stream bgzf_stream;

  void open_vcf_file_for_reading();
  void read_samples();
  bool read_record(bool SITES_ONLY = false);
  void read(bool SITES_ONLY = false); /** \brief Reads the VCF file. */
  void open_for_writing(long const n_threads = 1);
  void write_header(bool const is_dropping_genotypes = false);

  void write_record(Variant const & var,
                    std::string const & suffix = "",
                    bool const FILTER_ZERO_QUAL = false,
                    bool const is_dropping_genotypes = false);

  void write_tbi_index() const;
  void write_segments();
  void write(std::string const & region = ".", long const n_threads = 1); /** \brief Writes the VCF file. */

  void write_records(uint32_t region_begin,
                     uint32_t region_end,
                     bool FILTER_ZERO_QUAL,
                     bool is_dropping_genotypes,
                     std::vector<Variant> const & vars);

  void write_records(std::string const & region = ".",
                     bool FILTER_ZERO_QUAL = false,
                     bool const is_dropping_genotypes = false);

  void close_vcf_file();

  /** Common */
  void clear();

  // Adding data
  void add_hla_haplotypes(std::vector<Haplotype> & haplotypes,
                          std::vector<std::unordered_map<uint32_t, uint32_t> > const & allele_hap_gts_ptr);
  void add_haplotype(Haplotype & haplotype, int32_t phase_set = -1);

  void add_haplotypes_for_extraction(std::vector<HaplotypeCall> const & hap_calls, bool const is_splitting_vars);

  VCF_FILE_MODE filemode;
  std::string filename;

  std::vector<std::string> sample_names;
  std::vector<Variant> variants;

private:
  template <class Archive>
  void serialize(Archive & ar, unsigned int version);
};

std::vector<std::size_t> get_all_pos(std::string const & line, char const delim = '\t');
std::string get_string_at_tab_index(std::string const & line,
                                    std::vector<std::size_t> const & tabs,
                                    int index);

void save_vcf(Vcf const & vcf, std::string const & filename);
void load_vcf(Vcf & vcf, std::string const & filename, long n_batch);
bool append_vcf(Vcf & vcf, std::string const & filename, long n_batch);

} // namespace gyper
