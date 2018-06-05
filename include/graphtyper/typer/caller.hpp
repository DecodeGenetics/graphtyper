#pragma once

#include <string> // std::string
#include <vector> // std::vector


namespace gyper
{

void
call(std::vector<std::string> const & hts_path,
     std::string const & graph_path,
     std::string const & index_path,
     std::vector<std::string> const & regions,
     std::vector<std::string> const & segment_fasta_files,
     std::string const & output_dir
     );

void
discover_directly_from_bam(std::string const & graph_path,
                           std::vector<std::string> const & sams,
                           std::vector<std::string> const & regions,
                           std::string const & output_dir
  );

} // namespace gyper
