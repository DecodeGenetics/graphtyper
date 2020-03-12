#pragma once

#include <string> // std::string
#include <vector> // std::vector


namespace gyper
{


// returns the prefix to the output files
std::vector<std::string>
call(std::vector<std::string> const & hts_path,
     std::string const & graph_path,
     std::string const & index_path,
     std::string const & output_dir,
     std::string const & reference,
     std::string const & region,
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

} // namespace gyper
