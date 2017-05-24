#pragma once

#include <string> // std::string
#include <vector> // std::vector

#include <graphtyper/typer/vcf_writer.hpp>


namespace gyper
{

void
segment_calling(std::vector<std::string> const & segment_fasta_files,
                VcfWriter & writer,
                std::string const & segment_path,
                std::vector<std::string> const & samples
                );

}
