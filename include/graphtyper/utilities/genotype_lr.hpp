#include <string>


namespace gyper
{

void
genotype_lr_regions(std::string ref_path,
                    std::vector<std::string> const & sams,
                    std::vector<GenomicRegion> const & regions,
                    std::string const & output_path,
                    bool const is_copy_reference);

} // namespace gyper
