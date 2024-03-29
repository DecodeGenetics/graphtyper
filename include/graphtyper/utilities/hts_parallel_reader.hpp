#pragma once

#include <cassert>  // assert
#include <fstream>  // std::ifstream
#include <iostream> // std::cout, std::endl
#include <queue>    // std::priority_queue
#include <string>   // std::string
#include <vector>   // std::vector

#include <graphtyper/index/ph_index.hpp> // PHIndex
#include <graphtyper/typer/vcf_writer.hpp>
#include <graphtyper/utilities/hts_reader.hpp>
#include <graphtyper/utilities/hts_store.hpp>

namespace gyper
{
class Primers;
class VariantMap;

class HtsParallelReader
{
private:
  std::vector<HtsRecord> heap;
  HtsStore store;
  std::vector<std::string> samples; // List of sample names
  long num_rg{0};                   // Number of read groups

public:
  HtsParallelReader() = default;
  HtsParallelReader(HtsParallelReader const &) = delete;
  HtsParallelReader(HtsParallelReader &&) = delete;
  HtsParallelReader & operator=(HtsParallelReader const &) = delete;
  HtsParallelReader & operator=(HtsParallelReader &&) = delete;
  ~HtsParallelReader();

  // all hts files the parallel reader has
  std::vector<HtsReader> hts_files;

  // Opens multiple hts paths
  void open(std::vector<std::string> const & hts_file_paths,
            std::string const & reference = "",
            std::string const & region = ".");

  // Closes all hts files
  void close();

  // read the next hts record from the heap
  bool read_record(HtsRecord & hts_record);

  // move a record from 'from' to 'to'
  void move_record(HtsRecord & to, HtsRecord & from);

  // get a list of all sample names
  std::vector<std::string> const & get_samples() const;

  // get the sample and read group index for this record
  void get_sample_and_rg_index(long & sample_i, long & rg_i, HtsRecord const & hts_rec) const;

  // get the number of read groups
  long get_num_rg() const;

  // get the number of samples
  long get_num_samples() const;

  // get pointer to a header for the opened hts_files
  bam_hdr_t * get_header() const;
};

void parallel_reader_genotype_only(
  long const thread_id,
  std::string * out_path,
  std::vector<std::string> const * hts_paths_ptr,
  std::vector<double> const * avg_cov_ptr,
  std::string const * output_dir_ptr,
  std::string const * reference_fn_ptr,
  std::string const * region_ptr,
  PHIndex const * ph_index_ptr,
  Primers const * primers,
  std::vector<std::map<std::pair<uint16_t, uint16_t>, std::map<std::pair<uint16_t, uint16_t>, int8_t>>> * ph_ptr,
  bool const is_writing_calls_vcf,
  bool const is_writing_hap,
  std::vector<std::unordered_map<uint32_t, uint32_t>> * allele_hap_gts_ptr = nullptr);

void parallel_reader_with_discovery(std::string * out_path,
                                    std::vector<std::string> const * hts_paths_ptr,
                                    std::string const * output_dir_ptr,
                                    std::string const * reference_fn_ptr,
                                    std::string const * region_ptr,
                                    PHIndex const * ph_index_ptr,
                                    Primers const * primers,
                                    long const minimum_variant_support,
                                    double const minimum_variant_support_ratio,
                                    bool const is_writing_calls_vcf,
                                    bool const is_writing_hap);

void sam_merge(std::string const & output_sam, std::vector<std::string> const & input_sams);

} // namespace gyper
