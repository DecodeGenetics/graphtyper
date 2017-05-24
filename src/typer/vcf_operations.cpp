#include <cmath> // sqrt
#include <iostream> // std::cout, std::endl
#include <string> // std::string
#include <vector> // std::vector

#include <boost/log/trivial.hpp> // BOOST_LOG_TRIVIAL

#include <graphtyper/typer/vcf_operations.hpp>
#include <graphtyper/typer/vcf.hpp> // gyper::Vcf
#include <graphtyper/typer/var_stats.hpp> // gyper::join_strand_bias, gyper::split_bias_to_strings


namespace gyper
{

void
vcf_merge(std::vector<std::string> vcfs, std::string const & output)
{
  // Skip if the filename contains '*', cause that means there was no file which match that expansion
  vcfs.erase(std::remove_if(vcfs.begin(), vcfs.end(), [](std::string const & vcf){
      return std::count(vcf.begin(), vcf.end(), '*') > 0;
    }), vcfs.end());

  if (vcfs.size() == 0)
    return;

  gyper::Vcf vcf;
  vcf.open(READ_MODE, vcfs[0]);
  vcf.read(); // Read the entire file
  vcf.open(WRITE_MODE, output); // Change to write mode
  vcf.open_for_writing();

  std::vector<gyper::Vcf> next_vcfs(vcfs.size() - 1);

  // Open all VCFs and add sample names
  for (std::size_t i = 1; i < vcfs.size(); ++i)
  {
    if (std::count(vcfs[i].begin(), vcfs[i].end(), '*') > 0)
      continue;

    gyper::Vcf & next_vcf = next_vcfs[i - 1];
    next_vcf.open(READ_MODE, vcfs[i]);

    // Open the VCF file
    next_vcf.open_vcf_file_for_reading();

    // Read the sample names and add them
    next_vcf.read_samples();

    // Add samples names (we chose to copy them for assertions when reading the records)
    vcf.sample_names.insert(vcf.sample_names.end(), next_vcf.sample_names.begin(), next_vcf.sample_names.end());
  }

  BOOST_LOG_TRIVIAL(info) << "[graphtyper::vcf_operations::vcf_merge] Total number of samples read is " << vcf.sample_names.size();
  vcf.write_header(); // Now that we know all the sample names we can write the header

  // For each variant in the first VCF, add the calls from the other VCFs
  for (auto & var : vcf.variants)
  {
    // CR
    uint32_t number_of_clipped_reads = 0;

    {
      auto find_it = var.infos.find("CR");
      if (find_it != var.infos.end())
        number_of_clipped_reads = std::strtoull(find_it->second.c_str(), NULL, 10);
    }

    // Unaligned
    uint32_t number_of_unaligned_reads = 0;

    {
      auto find_it = var.infos.find("Unaligned");
      if (find_it != var.infos.end())
        number_of_unaligned_reads = std::strtoull(find_it->second.c_str(), NULL, 10);
    }

    // MQ
    uint64_t total_mapq_root = var.get_rooted_mapq(); // Keep track of the total mapping quality rooted

    // MQ0
    uint32_t mapq_zero_count = 0;

    {
      auto find_it = var.infos.find("MQ0");
      if (find_it != var.infos.end())
        mapq_zero_count = std::strtoull(find_it->second.c_str(), NULL, 10);
    }

    // SBF and SBR
    std::vector<uint32_t> strand_forward;
    std::vector<uint32_t> strand_reverse;

    std::vector<uint32_t> r1_strand_forward;
    std::vector<uint32_t> r1_strand_reverse;
    std::vector<uint32_t> r2_strand_forward;
    std::vector<uint32_t> r2_strand_reverse;

    std::vector<uint32_t> realignment_distance;
    std::vector<uint32_t> realignment_count;

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

      // Realignment count and distance
      find_it = var.infos.find("RACount");

      if (find_it != var.infos.end())
        realignment_count = split_bias_to_numbers(find_it->second);

      find_it = var.infos.find("RADist");

      if (find_it != var.infos.end())
        realignment_distance = split_bias_to_numbers(find_it->second);
    }

    for (auto & next_vcf : next_vcfs)
    {
      assert(next_vcf.variants.size() == 0);
      bool const SUCCESS = next_vcf.read_record();

      if (!SUCCESS)
      {
        BOOST_LOG_TRIVIAL(error) << "[graphtyper::vcf_operations::vcf_merge] There was a problem reading " << next_vcf.filename << ".";
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
      total_mapq_root += next_vcf.variants[0].get_rooted_mapq();

      // Get MQ0
      {
        auto find_it = next_vcf.variants[0].infos.find("MQ0");
        if (find_it != next_vcf.variants[0].infos.end())
          mapq_zero_count += std::strtoull(find_it->second.c_str(), NULL, 10);
      }

      // Get SBF and SBR
      {
        auto add_to_bias_lambda = [&](std::string id, std::vector<uint32_t> & bias)
        {
          auto find_it = next_vcf.variants[0].infos.find(id);

          if (find_it != next_vcf.variants[0].infos.end())
          {
            std::vector<uint32_t> split_nums = split_bias_to_numbers(find_it->second);

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
        add_to_bias_lambda("RACount", realignment_count);
        add_to_bias_lambda("RADist", realignment_distance);
      }

      std::move(next_vcf.variants[0].calls.begin(), next_vcf.variants[0].calls.end(), std::back_inserter(var.calls));
      std::move(next_vcf.variants[0].phase.begin(), next_vcf.variants[0].phase.end(), std::back_inserter(var.phase));
      next_vcf.variants.clear();
    }

    // Add strand bias, this must happend before the INFO is generated
    var.infos["SBF"] = join_strand_bias(strand_forward);
    var.infos["SBR"] = join_strand_bias(strand_reverse);

    var.infos["SBF1"] = join_strand_bias(r1_strand_forward);
    var.infos["SBF2"] = join_strand_bias(r2_strand_forward);
    var.infos["SBR1"] = join_strand_bias(r1_strand_reverse);
    var.infos["SBR2"] = join_strand_bias(r2_strand_reverse);

    var.infos["RACount"] = join_strand_bias(realignment_count);
    var.infos["RADist"] = join_strand_bias(realignment_distance);

    // MQ requirsed the generate INFOs
    var.generate_infos();

    // Add MQ
    {
      assert (var.infos.count("SeqDepth") == 1);
      std::string seq_depth_str = var.infos["SeqDepth"];

      if (seq_depth_str.size() > 0)
      {
        uint64_t const seq_depth = std::strtoull(seq_depth_str.c_str(), NULL, 10);

        if (seq_depth > 0)
          var.infos["MQ"] = std::to_string(static_cast<uint16_t>(sqrt(static_cast<double>(total_mapq_root) / static_cast<double>(seq_depth))));
      }
    }

    // Add MQ0
    var.infos["MQ0"] = std::to_string(mapq_zero_count);

    // Add CR
    var.infos["CR"] = std::to_string(number_of_clipped_reads);

    // Add number of unaligned reads
    var.infos["Unaligned"] = std::to_string(number_of_unaligned_reads);

    var.remove_uncalled_alleles();

    if (var.calls.size() != vcf.sample_names.size())
    {
      BOOST_LOG_TRIVIAL(error) << "[graphtyper::vcf_operations::vcf_merge] Number of calls a variant had did not matches the number of samples "
                               << var.calls.size() << " vs. " << vcf.sample_names.size();
      std::exit(1);
    }

    // Only write variants if they have any alternative alleles, after removing the uncalled alleles
    if (var.seqs.size() >= 2)
      vcf.write_record(var);

    var = Variant(); // Free memory
  }

  // Close all the files
  for (auto & next_vcf : next_vcfs)
    next_vcf.close_vcf_file();

  vcf.close_vcf_file();
}


void
vcf_concatenate(std::vector<std::string> const & vcfs, std::string const & output, bool const SKIP_SORT, std::string const & region)
{
  if (vcfs.size() == 0)
    return;

  gyper::Vcf vcf;

  if (SKIP_SORT)
  {
    vcf.open(WRITE_MODE, output);
    vcf.open_for_writing();

    for (std::size_t i = 0; i < vcfs.size(); ++i)
    {
      // Skip if the filename contains '*', cause that means there was no file which match that expansion
      if (std::count(vcfs[i].begin(), vcfs[i].end(), '*') > 0)
        continue;

      gyper::Vcf next_vcf;
      next_vcf.open(READ_MODE, vcfs[i]);
      next_vcf.open_vcf_file_for_reading();
      next_vcf.read_samples();

      // Copy sample names if this is the first VCF
      if (i == 0)
      {
        std::copy(next_vcf.sample_names.begin(), next_vcf.sample_names.end(), std::back_inserter(vcf.sample_names));
        vcf.write_header();
      }
      else if (next_vcf.sample_names.size() != vcf.sample_names.size())
      {
        BOOST_LOG_TRIVIAL(error) << "[graphtyper::vcf_operations::vcf_concatenate] The VCF file " << vcfs[i] << " has different amount of samples!";
        std::exit(1);
      }

      while (next_vcf.read_record())
      {
        // Add variants
        assert (vcf.variants.size() == 0);
        std::move(next_vcf.variants.begin(), next_vcf.variants.end(), std::back_inserter(vcf.variants));
        assert (vcf.variants.size() == 1);
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
    vcf.open(READ_MODE, vcfs[0]);
    vcf.read();
    BOOST_LOG_TRIVIAL(info) << "[graphtyper::vcf_operations::vcf_concatenate] Total number of samples read is " << vcf.sample_names.size();

    for (std::size_t i = 1; i < vcfs.size(); ++i)
    {
      // Skip if the filename contains '*', cause that means there was no file which match that expansion
      if (std::count(vcfs[i].begin(), vcfs[i].end(), '*') > 0)
        continue;

      gyper::Vcf next_vcf;
      next_vcf.open(READ_MODE, vcfs[i]);
      next_vcf.read();

      // Merge with the other VCF
      if (next_vcf.sample_names.size() != vcf.sample_names.size())
      {
        BOOST_LOG_TRIVIAL(error) << "[graphtyper::vcf_operations::vcf_concatenate] The VCF file " << vcfs[i] << " has different amount of samples!";
        std::exit(1);
      }

      // Add variants
      std::move(next_vcf.variants.begin(), next_vcf.variants.end(), std::back_inserter(vcf.variants));
    }

    // Regenerate the INFO scores
    for (auto & var : vcf.variants)
      var.generate_infos();

    // Change the filename to whatever the user wishes
    vcf.filename = output;

    // Change the filemode to writing
    vcf.set_filemode(WRITE_MODE);
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
  BOOST_LOG_TRIVIAL(info) << "[graphtyper::vcf_operations::vcf_break_down] Total number of samples read is " << vcf_in.sample_names.size();

  // Copy sample names
  std::copy(vcf_in.sample_names.begin(), vcf_in.sample_names.end(), std::back_inserter(vcf_out.sample_names));

  // Write header to output
  vcf_out.write_header();

  // Loop over the entire VCF in file, line by line
  while (true)
  {
    assert(vcf_in.variants.size() == 0);

    // Stop when we can't read any more records
    if (!vcf_in.read_record())
      break;

    vcf_in.variants[0].add_base_in_front(); // First add a single base in front
    assert(vcf_in.variants.size() == 1);
    assert(vcf_out.variants.size() == 0);

    // Make sure the number of calls matches the number of samples
    if (vcf_in.variants[0].calls.size() != vcf_in.sample_names.size())
    {
      BOOST_LOG_TRIVIAL(error) << "[graphtyper::vcf_operations::vcf_break_down] The number of calls, " << vcf_in.variants[0].calls.size()
                               << " did not match the number of samples, " << vcf_in.sample_names.size();
      std::exit(1);
    }

    // vcf_in.variants[0].remove_uncalled_alleles();
    vcf_out.variants = break_down_variant(std::move(vcf_in.variants[0]), 1 /*THRESHOLD*/);
    assert(vcf_out.variants.size() > 0);

    // Generate INFOs for the broken down variants
    vcf_out.post_process_variants(false /*normalize*/); // Do not normalize the variants, because then we could get duplicates

    // Write the records
    vcf_out.write_records(region, true /*FILTER_ZERO_QUAL*/);

    // Clear all the things we don't need anymore
    vcf_in.variants.clear();
    vcf_out.variants.clear();
  }

  vcf_in.close_vcf_file();
  vcf_out.close_vcf_file();
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
  std::copy(vcf_in.sample_names.begin(), vcf_in.sample_names.end(), std::back_inserter(vcf_out.sample_names));

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
    std::move(vcf_in.variants.begin(), vcf_in.variants.end(), std::back_inserter(vcf_out.variants));
    assert(vcf_out.variants.size() == 1);

    // Generate INFOs for the broken down variant
    for (auto & var : vcf_out.variants)
      var.generate_infos();

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
