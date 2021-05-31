#include <fstream>
#include <string>
#include <unordered_set>

#include <graphtyper/utilities/logging.hpp>

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/graph.hpp>
#include <graphtyper/typer/path.hpp>
#include <graphtyper/typer/genotype_paths.hpp>
#include <graphtyper/typer/primers.hpp>
#include <graphtyper/utilities/string.hpp>
#include <graphtyper/utilities/system.hpp>


namespace gyper
{

Primers::Primers(std::string const & primer_bedpe)
{
  read(primer_bedpe);
}


void
Primers::read(std::string const & primer_bedpe)
{
  if (!gyper::is_file(primer_bedpe))
  {
    print_log(log_severity::error, __HERE__, "] Could not find file ", primer_bedpe);
    std::exit(1);
  }

  // Open file
  std::ifstream fs(primer_bedpe.c_str());

  if (!fs.is_open())
  {
    print_log(log_severity::error, __HERE__, "] Could not open file ", primer_bedpe);
    std::exit(1);
  }

  std::string line;

  while (std::getline(fs, line))
  {
    std::vector<std::string_view> fields = split_on_delim(line, '\t');

    if (fields.size() < 6)
    {
      print_log(log_severity::error, __HERE__, " Expected at least 6 fields in Primer BEDPE, got ", fields.size());
      std::exit(1);
    }

    GenomicRegion left_region(fields[0], stoi64(fields[1]), stoi64(fields[2]));
    GenomicRegion right_region(fields[3], stoi64(fields[4]), stoi64(fields[5]));

    print_log(log_severity::debug, __HERE__, " Got primer regions ", left_region.to_string(), " "
                            , right_region.to_string());

    left.push_back(std::move(left_region));
    right.push_back(std::move(right_region));
  }
}


void
Primers::check(GenotypePaths & genos) const
{
  bool const is_reversed = (genos.flags & IS_SEQ_REVERSED) != 0;

  if (is_reversed)
    check_right(genos);
  else
    check_left(genos);
}


void
Primers::check_left(GenotypePaths & genos) const
{
  for (auto & path : genos.paths)
  {
    if (path.var_order.size() == 0)
      continue; // Nothing to do

    // Check if the path starts inside the primer region
    std::vector<Location> s_locs = graph.get_locations_of_a_position(path.start, path);

    // Loop over all left primer regions
    for (auto const & l : left)
    {
      // Get the absolute paths of the left primer region
      long constexpr PADDING = 5;
      auto const abs_begin = std::max(static_cast<long>(l.get_absolute_begin_position()) - PADDING, 1l);
      auto const abs_end = l.get_absolute_end_position();

      for (auto const & loc : s_locs)
      {
        auto pos = loc.node_order + loc.offset;

        // Detected that the path starts in the primer
        if (pos >= abs_begin && pos <= abs_end)
        {
          std::vector<uint32_t> var_orders = graph.get_var_orders(abs_begin, abs_end);

          for (long i = path.var_order.size() - 1; i >= 0l; --i)
          {
            if (std::find(var_orders.begin(), var_orders.end(), path.var_order[i]) != var_orders.end())
            {
              print_log(log_severity::debug, __HERE__, " LEFT Removed var_order=", path.var_order[i]);

              // Found a var_order inside the region
              path.erase_ref_support(i);
            }
          }
        }
      }
    }
  }
}


void
Primers::check_right(GenotypePaths & genos) const
{
  for (auto & path : genos.paths)
  {
    if (path.var_order.size() == 0)
      continue; // Nothing to do

    // Check if the path starts inside the primer region
    std::vector<Location> e_locs = graph.get_locations_of_a_position(path.end, path);

    // Loop over all right primer regions
    for (auto const & r : right)
    {
      // Get the absolute paths of the right primer region with some padding
      long constexpr PADDING = 5;
      long const abs_begin = r.get_absolute_begin_position();
      long const abs_end = r.get_absolute_end_position() + PADDING;

      for (auto const & loc : e_locs)
      {
        auto pos = loc.node_order + loc.offset;

        // Detected that the path ends in the primer
        if (pos >= abs_begin && pos <= abs_end)
        {
          std::vector<uint32_t> var_orders = graph.get_var_orders(abs_begin, abs_end);

          for (long i = path.var_order.size() - 1; i >= 0l; --i)
          {
            if (std::find(var_orders.begin(), var_orders.end(), path.var_order[i]) != var_orders.end())
            {
              print_log(log_severity::debug, __HERE__, " RIGHT Removed var_order=", path.var_order[i]);

              // Found a var_order inside the region
              path.erase_ref_support(i);
            }
          }
        }
      }
    }
  }
}


} // namespace gyper
