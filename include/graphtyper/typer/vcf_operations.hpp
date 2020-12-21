#pragma once

#include <string>
#include <vector>


namespace gyper
{

void vcf_merge(std::vector<std::string> & vcfs, std::string const & output);

void vcf_concatenate(std::vector<std::string> const & vcfs,
                     std::string const & output,
                     bool const SKIP_SORT,
                     bool const SITES_ONLY,
                     bool const WRITE_TBI,
                     std::string const & region);

void vcf_break_down(std::string const & vcf, std::string const & output, std::string const & region);

void
vcf_merge_and_filter(std::vector<std::string> const & vcfs,
                     std::string const & output,
                     std::map<std::pair<uint16_t, uint16_t>,
                              std::map<std::pair<uint16_t, uint16_t>, int8_t> > const & ph);

void vcf_merge_and_break(std::vector<std::string> const & vcfs,
                         std::string const & output,
                         std::string const & region,
                         bool const FILTER_ZERO_QUAL,
                         bool const force_no_variant_overlap,
                         bool const force_no_break_down);

void vcf_update_info(std::string const & vcf, std::string const & output);

} // namespace gyper
