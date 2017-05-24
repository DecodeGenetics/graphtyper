#pragma once

#include <string> // std::string
#include <vector> // std::vector


namespace gyper
{

void
call(std::vector<std::string> hts_path,
     std::string graph_path,
     std::string index_path,
     std::vector<std::string> const & regions,
     std::vector<std::string> const & segment_fasta_files,
     std::string const & output_dir
     );

} // namespace gyper
