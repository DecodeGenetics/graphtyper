#pragma once

#include <string> // std::string
#include <vector> // std::vector


namespace gyper
{

class PHIndex;
class Vcf;
class Primers;

// returns the prefix to the output files
std::vector<std::string>
call(std::vector<std::string> const & hts_path,
     std::vector<double> const & avg_cov_by_readlen,
     std::string const & graph_path,
     PHIndex const & ph_index,
     std::string const & output_dir,
     std::string const & reference,
     std::string const & region,
     Primers const * primers,
     std::map<std::pair<uint16_t, uint16_t>, std::map<std::pair<uint16_t, uint16_t>, int8_t> > & ph,
     long const minimum_variant_support,
     double const minimum_variant_support_ratio,
     bool const is_writing_calls_vcf,
     bool const is_discovery,
     bool const is_writing_hap);


// returns the written variant maps
std::vector<std::string>
discover_directly_from_bam(std::string const & graph_path,
                           std::vector<std::string> const & hts_paths,
                           std::string const & region_str,
                           std::string const & output_dir,
                           long minimum_variant_support,
                           double minimum_variant_support_ratio);

void
streamlined_discovery(std::vector<std::string> const & hts_paths,
                      std::string const & reference_fn,
                      std::string const & region_str,
                      gyper::Vcf & vcf);

} // namespace gyper
