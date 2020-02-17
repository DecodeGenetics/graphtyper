#pragma once

#include <string>
#include <vector>

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/graph.hpp>

namespace gyper
{

void construct_graph(std::string const & reference_filename,
                     std::string const & vcf_filename,
                     std::string const & region,
                     bool is_sv_graph = false,
                     bool use_absolute_positions = true,
                     bool check_index = true);

} // namespace gyper
