#pragma once

#include <string>
#include <vector>

#include <graphtyper/graph/genomic_region.hpp>


namespace gyper
{

void
create_dir(std::string const & dir, unsigned mode = 0755);

std::string
create_temp_dir(GenomicRegion const & region);

void
remove_file_tree(std::string const & path);

bool
is_file(std::string const & filename);

void
check_file_exists(std::string const & filename);

void
check_file_exists_or_empty(std::string const & filename);

bool
is_directory(std::string const & filename);

bool
is_defined_in_env(std::string const & var);

} // namespace gyper
