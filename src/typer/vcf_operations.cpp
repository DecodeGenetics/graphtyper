#include <cmath> // sqrt
#include <string> // std::string
#include <sstream> // std::ostringstream
#include <vector> // std::vector

#include <boost/log/trivial.hpp> // BOOST_LOG_TRIVIAL

#include <graphtyper/graph/absolute_position.hpp>
#include <graphtyper/graph/genomic_region.hpp>
#include <graphtyper/typer/variant.hpp> // gyper::break_down_variant
#include <graphtyper/typer/vcf_operations.hpp>
#include <graphtyper/typer/vcf.hpp> // gyper::Vcf
#include <graphtyper/typer/var_stats.hpp> // gyper::join_strand_bias, gyper::split_bias_to_strings
#include <graphtyper/utilities/options.hpp> // gyper::Options


namespace gyper
{

void
vcf_merge(std::vector<std::string> & vcfs, std::string const & output)
{
  // Skip if the filename contains '*'
  vcfs.erase(std::remove_if(vcfs.begin(), vcfs.end(), [](std::string const & vcf){
      return std::count(vcf.begin(), vcf.end(), '*') > 0;
    }), vcfs.end());

  if (vcfs.size() == 0)
    return;

  gyper::Vcf vcf;
  vcf.open(READ_MODE, vcfs.at(0));
  vcf.read(); // Read the entire file
  vcf.open(WRITE_MODE, output); // Change to write mode
  vcf.open_for_writing();

  std::vector<gyper::Vcf> next_vcfs(vcfs.size() - 1);

  // For checking if we have duplicated IDs
  long dup = -1l;
  uint32_t old_abs_pos = static_cast<uint32_t>(-1);
  std::string old_variant_type = "";

  // Open all VCFs and add sample names
  for (std::size_t i = 1; i < vcfs.size(); ++i)
  {
    if (std::count(vcfs[i].begin(), vcfs[i].end(), '*') > 0)
      continue;

    gyper::Vcf & next_vcf = next_vcfs[i - 1];
    next_vcf.open(READ_MODE, vcfs[i]);

    // Open the VCF file
    next_vcf.open_vcf_file_for_reading();

    if (next_vcf.bgzf_in)
    {
      // Read the sample names and add them
      next_vcf.read_samples();

      // Add samples names (we chose to copy them for assertions when reading the records)
      vcf.sample_names.insert(vcf.sample_names.end(),
                              next_vcf.sample_names.begin(),
                              next_vcf.sample_names.end());
    }
  }

  BOOST_LOG_TRIVIAL(debug) << __HERE__ << " Total number of samples read is "
                           << vcf.sample_names.size();

  vcf.write_header(); // Now that we know all the sample names we can write the header

  // For each variant in the first VCF, add the calls from the other VCFs
  for (auto & var : vcf.variants)
  {
    // Only output variant if it is in the genotyping region
    // CR
    uint32_t number_of_clipped_reads = 0;

    {
      auto find_it = var.infos.find("CR");

      if (find_it != var.infos.end())
        number_of_clipped_reads = std::strtoull(find_it->second.c_str(), NULL, 10);
    }

    // MQ. Keep track of the total mapping quality rooted/squared
    uint64_t mapq_squared = 0;

    {
      auto find_it = var.infos.find("MQsquared");

      if (find_it != var.infos.end())
        mapq_squared = std::strtoull(find_it->second.c_str(), NULL, 10);
    }

    // SBF and SBR
    std::vector<uint32_t> strand_forward;
    std::vector<uint32_t> strand_reverse;

    std::vector<uint32_t> r1_strand_forward;
    std::vector<uint32_t> r1_strand_reverse;
    std::vector<uint32_t> r2_strand_forward;
    std::vector<uint32_t> r2_strand_reverse;

//    std::vector<uint32_t> realignment_distance;
//    std::vector<uint32_t> realignment_count;

    {
      auto find_it = var.infos.find("SBF");

      if (find_it != var.infos.end())
        strand_forward = split_bias_to_numbers(find_it->second);

      find_it = var.infos.find("SBR");

      if (find_it != var.infos.end())
        strand_reverse = split_bias_to_numbers(find_it->second);

      // Read specific strand bias
      find_it = var.infos.find("SBF1");

      if (find_it != var.infos.end())
        r1_strand_forward = split_bias_to_numbers(find_it->second);

      find_it = var.infos.find("SBF2");

      if (find_it != var.infos.end())
        r2_strand_forward = split_bias_to_numbers(find_it->second);

      find_it = var.infos.find("SBR1");

      if (find_it != var.infos.end())
        r1_strand_reverse = split_bias_to_numbers(find_it->second);

      find_it = var.infos.find("SBR2");

      if (find_it != var.infos.end())
        r2_strand_reverse = split_bias_to_numbers(find_it->second);

//      // Realignment count and distance
//      find_it = var.infos.find("RACount");
//
//      if (find_it != var.infos.end())
//        realignment_count = split_bias_to_numbers(find_it->second);
//
//      find_it = var.infos.find("RADist");
//
//      if (find_it != var.infos.end())
//        realignment_distance = split_bias_to_numbers(find_it->second);
    }

    for (auto & next_vcf : next_vcfs)
    {
      if (!next_vcf.bgzf_in)
        continue;

      assert(next_vcf.variants.size() == 0);
      bool const SUCCESS = next_vcf.read_record();

      if (!SUCCESS)
      {
        BOOST_LOG_TRIVIAL(error) << __HERE__ << " There was a problem reading "
                                 << next_vcf.filename << ".";
        std::exit(1);
      }

      assert(next_vcf.variants.size() == 1);

      // Get CR
      {
        auto find_it = next_vcf.variants[0].infos.find("CR");

        if (find_it != next_vcf.variants[0].infos.end())
          number_of_clipped_reads += std::strtoull(find_it->second.c_str(), NULL, 10);
      }

      // Get MQ
      {
        auto find_it = next_vcf.variants[0].infos.find("MQsquared");

        if (find_it != next_vcf.variants[0].infos.end())
          mapq_squared += std::strtoull(find_it->second.c_str(), NULL, 10);
      }


      // Get SBF and SBR
      {
        auto add_to_bias_lambda = [&](std::string id, std::vector<uint32_t> & bias)
                                  {
                                    auto find_it = next_vcf.variants[0].infos.find(id);

                                    if (find_it != next_vcf.variants[0].infos.end())
                                    {
                                      std::vector<uint32_t> split_nums =
                                        split_bias_to_numbers(find_it->second);

                                      for (std::size_t i = 0; i < split_nums.size(); ++i)
                                      {
                                        if (i < strand_forward.size())
                                          bias[i] += split_nums[i];
                                        else
                                          bias.push_back(split_nums[i]);
                                      }
                                    }
                                  };

        add_to_bias_lambda("SBF", strand_forward);
        add_to_bias_lambda("SBR", strand_reverse);
        add_to_bias_lambda("SBF1", r1_strand_forward);
        add_to_bias_lambda("SBF2", r2_strand_forward);
        add_to_bias_lambda("SBR1", r1_strand_reverse);
        add_to_bias_lambda("SBR2", r2_strand_reverse);
//        add_to_bias_lambda("RACount", realignment_count);
//        add_to_bias_lambda("RADist", realignment_distance);
      }

      std::move(next_vcf.variants[0].calls.begin(),
                next_vcf.variants[0].calls.end(),
                std::back_inserter(var.calls));

      next_vcf.variants.clear();
    }

    // Add strand bias, this must happend before the INFO is generated
    var.infos["SBF"] = join_strand_bias(strand_forward);
    var.infos["SBR"] = join_strand_bias(strand_reverse);

    var.infos["SBF1"] = join_strand_bias(r1_strand_forward);
    var.infos["SBF2"] = join_strand_bias(r2_strand_forward);
    var.infos["SBR1"] = join_strand_bias(r1_strand_reverse);
    var.infos["SBR2"] = join_strand_bias(r2_strand_reverse);

//    var.infos["RACount"] = join_strand_bias(realignment_count);
//    var.infos["RADist"] = join_strand_bias(realignment_distance);

    // Add MQsquared
    var.infos["MQsquared"] = std::to_string(mapq_squared);

    // Add CR
    var.infos["CR"] = std::to_string(number_of_clipped_reads);

    // generate the rest of the INFOs
    var.generate_infos();

    if (var.calls.size() != vcf.sample_names.size())
    {
      BOOST_LOG_TRIVIAL(error) << __HERE__ << " Number of calls "
                               << "a variant had did not matches the number of samples "
                               << var.calls.size() << " vs. " << vcf.sample_names.size();
      std::exit(1);
    }

    // Only write variants if they have any alternative alleles
    assert(var.seqs.size() >= 2);
    std::string const & new_variant_type = var.determine_variant_type();

    // Check if this is duplicated ID, and if so make them unique by adding a suffix to the ID
    if (old_abs_pos != var.abs_pos || old_variant_type != new_variant_type)
    {
      vcf.write_record(var, "" /*suffix*/);
      dup = -1;
    }
    else
    {
      ++dup;
      assert(dup >= 0l);
      vcf.write_record(var, std::string(".") + std::to_string(dup) /*suffix*/);
    }

    old_abs_pos = var.abs_pos;
    old_variant_type = new_variant_type;

    var = Variant(); // Free memory
  }

  // Close all the files
  for (auto & next_vcf : next_vcfs)
    next_vcf.close_vcf_file();

  vcf.close_vcf_file();
  //vcf.write_tbi_index();
}


void
vcf_merge_and_break(std::vector<std::string> const & vcfs,
                    std::string const & output,
                    std::string const & region,
                    bool const FILTER_ZERO_QUAL,
                    bool const force_no_variant_overlapping)
{
  if (vcfs.size() == 0)
    return;

  long const ploidy = Options::const_instance()->ploidy;
  GenomicRegion genomic_region(region);
  uint32_t const region_begin = 1 + absolute_pos.get_absolute_position(genomic_region.chr,
                                                                       genomic_region.begin
                                                                       );

  uint32_t const region_end = absolute_pos.get_absolute_position(genomic_region.chr,
                                                                 genomic_region.end);

  auto const & vcf_fn = vcfs[0];
  Vcf vcf;
  load_vcf(vcf, vcf_fn, 0);
  long n_batch{1};

  // Read the entire file
  while (append_vcf(vcf, vcf_fn, n_batch))
    ++n_batch;

  vcf.open(WRITE_MODE, output); // Change to write mode
  vcf.open_for_writing();
  BOOST_LOG_TRIVIAL(debug) << __HERE__ << " Read " << vcf.variants.size() << " variants.";

  std::vector<gyper::Vcf> next_vcfs(vcfs.size() - 1);
  n_batch = 1;

  // Open all VCFs and add sample names
  for (long i = 1; i < static_cast<long>(vcfs.size()); ++i)
  {
    auto const & vcf_fn = vcfs[i];
    gyper::Vcf & next_vcf = next_vcfs[i - 1];
    load_vcf(next_vcf, vcf_fn, 0);

    assert(next_vcf.variants.size() == next_vcfs[0].variants.size()); // all bins have the same number of variants

    // Add the samples on the first batch
    vcf.sample_names.insert(vcf.sample_names.end(),
                            next_vcf.sample_names.begin(),
                            next_vcf.sample_names.end());


  }

  // We have all the samples, write the header now
  vcf.write_header();

  long reach{-1};
  std::vector<Variant> broken_vars; // broken down variants
  long v_next{0}; // number of variants analyzed in previous batches

  // lambda function for updating the reach
  auto update_reach =
    [&reach](std::vector<Variant> const & new_variants) -> void
    {
      for (auto const & v : new_variants)
      {
        long curr_reach = v.seqs[0].size() + v.abs_pos;

        if (curr_reach > reach)
          reach = curr_reach;
      }
    };

  for (long v = 0; v < static_cast<long>(vcf.variants.size()); ++v)
  {
    auto & var = vcf.variants[v];

    // CR
    uint32_t number_of_clipped_reads = 0;

    {
      auto find_it = var.infos.find("CR");

      if (find_it != var.infos.end())
        number_of_clipped_reads = std::strtoull(find_it->second.c_str(), NULL, 10);
    }

    // MQ. Keep track of the total mapping quality rooted/squared
    uint64_t mapq_squared = 0;

    {
      auto find_it = var.infos.find("MQsquared");

      if (find_it != var.infos.end())
        mapq_squared = std::strtoull(find_it->second.c_str(), NULL, 10);
    }

    // SBF and SBR
    std::vector<uint32_t> strand_forward;
    std::vector<uint32_t> strand_reverse;

    std::vector<uint32_t> r1_strand_forward;
    std::vector<uint32_t> r1_strand_reverse;
    std::vector<uint32_t> r2_strand_forward;
    std::vector<uint32_t> r2_strand_reverse;

    {
      auto find_it = var.infos.find("SBF");

      if (find_it != var.infos.end())
        strand_forward = split_bias_to_numbers(find_it->second);

      find_it = var.infos.find("SBR");

      if (find_it != var.infos.end())
        strand_reverse = split_bias_to_numbers(find_it->second);

      // Read specific strand bias
      find_it = var.infos.find("SBF1");

      if (find_it != var.infos.end())
        r1_strand_forward = split_bias_to_numbers(find_it->second);

      find_it = var.infos.find("SBF2");

      if (find_it != var.infos.end())
        r2_strand_forward = split_bias_to_numbers(find_it->second);

      find_it = var.infos.find("SBR1");

      if (find_it != var.infos.end())
        r1_strand_reverse = split_bias_to_numbers(find_it->second);

      find_it = var.infos.find("SBR2");

      if (find_it != var.infos.end())
        r2_strand_reverse = split_bias_to_numbers(find_it->second);
    }

    // Trigger read next batch
    if ((v - v_next) == static_cast<long>(next_vcfs[0].variants.size()))
    {
      BOOST_LOG_TRIVIAL(debug) << __HERE__ << " Updating next_vcfs";

      v_next += static_cast<long>(next_vcfs[0].variants.size());

      for (long i = 1; i < static_cast<long>(vcfs.size()); ++i)
      {
        auto const & vcf_fn = vcfs[i];
        Vcf & next_vcf = next_vcfs[i - 1];
        next_vcf.variants.clear(); // clear previous batch
        bool ret = append_vcf(next_vcf, vcf_fn, n_batch);

        if (!ret)
        {
          BOOST_LOG_TRIVIAL(error) << __HERE__ << " Could not find file " << vcf_fn;
          std::exit(1);
        }
      }

      ++n_batch;
    }

    for (auto & next_vcf : next_vcfs)
    {
      assert((v - v_next) >= 0);
      assert((v - v_next) < static_cast<long>(next_vcf.variants.size()));
      auto & next_vcf_var = next_vcf.variants[v - v_next];

      // Get CR
      {
        auto find_it = next_vcf_var.infos.find("CR");

        if (find_it != next_vcf_var.infos.end())
          number_of_clipped_reads += std::strtoull(find_it->second.c_str(), NULL, 10);
      }

      // Get MQ
      {
        auto find_it = next_vcf_var.infos.find("MQsquared");

        if (find_it != next_vcf_var.infos.end())
          mapq_squared += std::strtoull(find_it->second.c_str(), NULL, 10);
      }


      // Get SBF and SBR
      {
        auto add_to_bias_lambda =
          [&](std::string id, std::vector<uint32_t> & bias)
          {
            auto find_it = next_vcf_var.infos.find(id);

            if (find_it != next_vcf_var.infos.end())
            {
              std::vector<uint32_t> split_nums = split_bias_to_numbers(find_it->second);

              for (long i = 0; i < static_cast<long>(split_nums.size()); ++i)
              {
                if (i < static_cast<long>(strand_forward.size()))
                  bias[i] += split_nums[i];
                else
                  bias.push_back(split_nums[i]);
              }
            }
          };

        add_to_bias_lambda("SBF", strand_forward);
        add_to_bias_lambda("SBR", strand_reverse);
        add_to_bias_lambda("SBF1", r1_strand_forward);
        add_to_bias_lambda("SBF2", r2_strand_forward);
        add_to_bias_lambda("SBR1", r1_strand_reverse);
        add_to_bias_lambda("SBR2", r2_strand_reverse);
      }

      std::move(next_vcf_var.calls.begin(),
                next_vcf_var.calls.end(),
                std::back_inserter(var.calls));
    }

    // Add strand bias, this must happend before the INFO is generated
    var.infos["SBF"] = join_strand_bias(strand_forward);
    var.infos["SBR"] = join_strand_bias(strand_reverse);

    var.infos["SBF1"] = join_strand_bias(r1_strand_forward);
    var.infos["SBF2"] = join_strand_bias(r2_strand_forward);
    var.infos["SBR1"] = join_strand_bias(r1_strand_reverse);
    var.infos["SBR2"] = join_strand_bias(r2_strand_reverse);

    // Add MQsquared
    var.infos["MQsquared"] = std::to_string(mapq_squared);

    // Add CR
    var.infos["CR"] = std::to_string(number_of_clipped_reads);

    if (var.calls.size() != vcf.sample_names.size())
    {
      BOOST_LOG_TRIVIAL(error) << "Number of calls a variant had did not matches the number of samples "
                               << var.calls.size() << " vs. " << vcf.sample_names.size();
      std::exit(1);
    }

    // break down the merged variants
    bool const is_no_variant_overlapping{Options::const_instance()->no_variant_overlapping ||
                                         force_no_variant_overlapping};

    bool const is_all_biallelic{Options::const_instance()->is_all_biallelic};
    std::vector<Variant> new_variants = break_down_variant(std::move(var),
                                                           reach,
                                                           is_no_variant_overlapping,
                                                           is_all_biallelic);
    assert(new_variants.size() > 0);

    for (auto & new_var : new_variants)
    {
      // If we have not processed this variant before, do so now
      if (!new_var.is_info_generated)
      {
        new_var.normalize();

        if (ploidy > 2)
          new_var.update_camou_phred(ploidy);

        new_var.generate_infos();

        // Remove MQsquared
        new_var.infos.erase("MQsquared");
        new_var.infos.erase("PS");
      }
    }

    update_reach(new_variants);
    std::move(new_variants.begin(), new_variants.end(), std::back_inserter(broken_vars));

    long constexpr W = 500; // Print variants that are more than W bp before the newest one

    if ((broken_vars[0].abs_pos + 2 * W) < broken_vars[broken_vars.size() - 1].abs_pos)
    {
      // Make sure we do no print outside of the region
      uint32_t const reg_end =
        std::min(region_end, static_cast<uint32_t>(broken_vars[broken_vars.size() - 1].abs_pos - W));

      vcf.write_records(region_begin,
                        reg_end,
                        FILTER_ZERO_QUAL,
                        broken_vars);

      // Remove variants that were written (or before the region, since they will never be printed)
      broken_vars.erase(
        std::remove_if(broken_vars.begin(),
                       broken_vars.end(),
                       [&](Variant const & v)
        {
          return v.abs_pos <= reg_end || v.abs_pos < region_begin;
        }),
        broken_vars.end()
        );
    }

    var = Variant(); // Free memory
  }

  // Write the remaining broken variants
  vcf.write_records(region_begin,
                    region_end,
                    FILTER_ZERO_QUAL,
                    broken_vars);

  vcf.close_vcf_file();
  vcf.write_tbi_index();
}


void
vcf_concatenate(std::vector<std::string> const & vcfs,
                std::string const & output,
                bool const SKIP_SORT,
                bool const SITES_ONLY,
                std::string const & region)
{
  gyper::Vcf vcf;

  if (vcfs.size() == 0)
  {
    vcf.open(WRITE_MODE, output);
    vcf.open_for_writing();
    vcf.write_header();
    vcf.close_vcf_file();
    return;
  }

  if (SKIP_SORT)
  {
    vcf.open(WRITE_MODE, output);
    vcf.open_for_writing();

    if (SITES_ONLY)
    {
      vcf.sample_names.clear();

      for (auto & var : vcf.variants)
        var.calls.clear();
    }

    for (std::size_t i = 0; i < vcfs.size(); ++i)
    {
      // Skip if the filename contains '*'
      if (std::count(vcfs[i].begin(), vcfs[i].end(), '*') > 0)
        continue;

      gyper::Vcf next_vcf;
      next_vcf.open(READ_MODE, vcfs[i]);
      next_vcf.open_vcf_file_for_reading();
      next_vcf.read_samples();

      if (SITES_ONLY)
      {
        next_vcf.sample_names.clear();

        for (auto & var : next_vcf.variants)
          var.calls.clear();
      }

      // Copy sample names if this is the first VCF
      if (i == 0)
      {
        std::copy(next_vcf.sample_names.begin(), next_vcf.sample_names.end(),
                  std::back_inserter(vcf.sample_names));
        vcf.write_header();
      }
      else if (next_vcf.sample_names.size() != vcf.sample_names.size())
      {
        BOOST_LOG_TRIVIAL(error) << "[graphtyper::vcf_operations] The VCF file "
                                 << vcfs[i]
                                 << " has different amount of samples! ("
                                 << next_vcf.sample_names.size()
                                 << " but not " << vcf.sample_names.size() << ")\n";
        std::exit(1);
      }

      while (next_vcf.read_record(SITES_ONLY))
      {
        // Add variants
        assert(vcf.variants.size() == 0);
        std::move(next_vcf.variants.begin(), next_vcf.variants.end(),
                  std::back_inserter(vcf.variants));
        assert(vcf.variants.size() == 1);

        if (vcf.sample_names.size() > 0)
          vcf.variants[0].generate_infos();

        vcf.write_records();

        // Sure all is cleared now
        vcf.variants.clear();
        next_vcf.variants.clear();
      }

      next_vcf.close_vcf_file();
    }

    vcf.close_vcf_file();
  }
  else // Also sort
  {
    vcf.open(WRITE_MODE, output);
    vcf.open_for_writing();

    if (SITES_ONLY)
      vcf.sample_names.clear();

    BOOST_LOG_TRIVIAL(debug) << "[graphtyper::vcf_operations] Total number of samples read is "
                             << vcf.sample_names.size();

    for (std::size_t i = 0; i < vcfs.size(); ++i)
    {
      // Skip if the filename contains '*'
      if (std::count(vcfs[i].begin(), vcfs[i].end(), '*') > 0)
        continue;

      gyper::Vcf next_vcf;
      next_vcf.open(READ_MODE, vcfs[i]);
      next_vcf.read(SITES_ONLY);

      if (SITES_ONLY)
      {
        next_vcf.sample_names.clear();

        for (auto & var : next_vcf.variants)
          var.calls.clear();
      }
      else if (i == 0)
      {
        // First VCF gets
        vcf.sample_names = next_vcf.sample_names;
      }

      // Merge with the other VCF
      if (next_vcf.sample_names.size() != vcf.sample_names.size())
      {
        BOOST_LOG_TRIVIAL(error) << "[graphtyper::vcf_operations] The VCF file "
                                 << vcfs[i]
                                 << " has different amount of samples! ("
                                 << next_vcf.sample_names.size()
                                 << " but not " << vcf.sample_names.size() << ")";
        std::exit(1);
      }

      // Add variants
      std::move(next_vcf.variants.begin(),
                next_vcf.variants.end(),
                std::back_inserter(vcf.variants));
    }

    // Regenerate the INFO scores
    if (vcf.sample_names.size() > 0)
    {
      for (auto & var : vcf.variants)
        var.generate_infos();
    }

    vcf.write(region);
  }
}


void
vcf_break_down(std::string const & vcf, std::string const & output, std::string const & region)
{
  gyper::Vcf vcf_in;
  vcf_in.open(READ_MODE, vcf);

  gyper::Vcf vcf_out;
  vcf_out.open(WRITE_MODE, output);

  // Open the VCF files
  vcf_in.open_vcf_file_for_reading();
  vcf_out.open_for_writing();

  // Read the sample names and add them
  vcf_in.read_samples();
  BOOST_LOG_TRIVIAL(debug) << "[graphtyper::vcf_break_down] "
                           << "Total number of samples read is "
                           << vcf_in.sample_names.size();

  // Copy sample names
  std::copy(vcf_in.sample_names.begin(),
            vcf_in.sample_names.end(),
            std::back_inserter(vcf_out.sample_names)
            );

  GenomicRegion genomic_region(region);
  uint32_t const region_begin = 1 + absolute_pos.get_absolute_position(genomic_region.chr,
                                                                       genomic_region.begin
                                                                       );

  uint32_t const region_end = absolute_pos.get_absolute_position(genomic_region.chr,
                                                                 genomic_region.end
                                                                 );

  // Read first record
  bool not_at_end = vcf_in.read_record();

  // Write header to output
  vcf_out.write_header();
  BOOST_LOG_TRIVIAL(debug) << "[graphtyper::vcf_break_down] VCF header written";

  long reach = -1; // Indicates how long the previous variants reached

  auto update_reach =
    [&reach](std::vector<Variant> const & new_variants)
    {
      for (auto const & var : new_variants)
      {
        long curr_reach = var.seqs[0].size() + var.abs_pos;

        if (curr_reach > reach)
          reach = curr_reach;
      }
    };

  // Loop over the entire VCF in file, line by line
  for (; not_at_end; not_at_end = vcf_in.read_record())
  {
    vcf_in.variants[0].add_base_in_front(); // First add a single base in front
    assert(vcf_in.variants.size() == 1);

    // Make sure the number of calls matches the number of samples
    if (vcf_in.variants[0].calls.size() != vcf_in.sample_names.size())
    {
      BOOST_LOG_TRIVIAL(error) << "[graphtyper::vcf_operations::vcf_break_down]"
                               << "The number of calls, "
                               << vcf_in.variants[0].calls.size()
                               << " did not match the number of samples, "
                               << vcf_in.sample_names.size();

      std::exit(1);
    }

    bool const is_no_variant_overlapping{Options::const_instance()->no_variant_overlapping};
    bool const is_all_biallelic{Options::const_instance()->is_all_biallelic};
    std::vector<Variant> new_variants = break_down_variant(std::move(vcf_in.variants[0]),
                                                           reach,
                                                           is_no_variant_overlapping,
                                                           is_all_biallelic);

    update_reach(new_variants);
    std::move(new_variants.begin(), new_variants.end(), std::back_inserter(vcf_out.variants));
    assert(vcf_out.variants.size() > 0);

    long constexpr W = 250; // Print variants that are more than W bp before the newest one

    if ((vcf_out.variants.front().abs_pos + 3 * W) < vcf_in.variants.back().abs_pos)
    {
      // Generate infos
      for (auto & var : vcf_out.variants)
      {
        var.normalize();
        var.generate_infos();
      }

      // Make sure we do no print outside of the region
      uint32_t const reg_end = std::min(region_end,
                                        static_cast<uint32_t>(vcf_in.variants.back().abs_pos - W)
                                        );

      vcf_out.write_records(region_begin,
                            reg_end,
                            true /*FILTER_ZERO_QUAL*/,
                            vcf_out.variants
                            );

      // Remove variants that were written (or before the region, since they will never be printed)
      vcf_out.variants.erase(
        std::remove_if(vcf_out.variants.begin(),
                       vcf_out.variants.end(),
                       [&](Variant const & v)
        {
          return v.abs_pos <= reg_end;
        }),

        vcf_out.variants.end()
        );
    }

    // Write the records
    assert(vcf_out.variants[0].calls.size() == vcf_out.sample_names.size());

    // Clear all the things we don't need anymore
    vcf_in.variants.clear();
  }

  // Generate infos
  for (auto & var : vcf_out.variants)
  {
    var.normalize();
    var.generate_infos();
  }

  vcf_out.write_records(region_begin, region_end, true /*FILTER_ZERO_QUAL*/, vcf_out.variants);
  vcf_in.close_vcf_file();
  vcf_out.close_vcf_file();
  vcf_out.write_tbi_index();
}


void
vcf_update_info(std::string const & vcf, std::string const & output)
{
  gyper::Vcf vcf_in;
  vcf_in.open(READ_MODE, vcf);

  gyper::Vcf vcf_out;
  vcf_out.open(WRITE_MODE, output);

  // Open the VCF files
  vcf_in.open_vcf_file_for_reading();
  vcf_out.open_for_writing();

  // Read the sample names and add them
  vcf_in.read_samples();

  // Copy sample names
  std::copy(vcf_in.sample_names.begin(), vcf_in.sample_names.end(),
            std::back_inserter(vcf_out.sample_names));

  // Write header to output
  vcf_out.write_header();

  // Loop over the entire VCF in file, line by line
  while (true)
  {
    assert(vcf_in.variants.size() == 0);

    // Stop when we can't read any more records
    if (!vcf_in.read_record())
      break;

    assert(vcf_in.variants.size() == 1);
    assert(vcf_out.variants.size() == 0);

    std::move(vcf_in.variants.begin(),
              vcf_in.variants.end(),
              std::back_inserter(vcf_out.variants)
              );

    assert(vcf_out.variants.size() == 1);

    // Generate INFOs
    if (vcf_out.sample_names.size() > 0)
    {
      for (auto & var : vcf_out.variants)
        var.generate_infos();
    }

    // Write the record
    vcf_out.write_records();

    // Clear all the things we don't need anymore
    vcf_in.variants.clear();
    vcf_out.variants.clear();
  }

  vcf_in.close_vcf_file();
  vcf_out.close_vcf_file();
}


} // namespace gyper
