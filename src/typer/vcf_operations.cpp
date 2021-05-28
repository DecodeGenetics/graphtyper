#include <cmath> // sqrt
#include <string> // std::string
#include <sstream> // std::ostringstream
#include <vector> // std::vector

#include <graphtyper/utilities/logging.hpp> // BOOST_LOG_TRIVIAL

#include <graphtyper/graph/absolute_position.hpp>
#include <graphtyper/graph/genomic_region.hpp>
#include <graphtyper/typer/variant.hpp> // gyper::break_down_variant
#include <graphtyper/typer/vcf_operations.hpp>
#include <graphtyper/typer/vcf.hpp> // gyper::Vcf
#include <graphtyper/typer/var_stats.hpp> // gyper::join_strand_bias, gyper::split_bias_to_strings
#include <graphtyper/utilities/options.hpp> // gyper::Options


namespace gyper
{

// vcf should be constructed with gyper::Vcf and no functions called
void
vcf_merge_and_return(gyper::Vcf & vcf, std::vector<std::string> & vcfs, std::string const & output)
{
  if (vcfs.size() == 0)
    return;

  auto const & copts = *(Options::const_instance());
  //long const ploidy = copts.ploidy;
  auto const & vcf_fn = vcfs[0];
  load_vcf(vcf, vcf_fn, 0);
  long n_batch{1};

  // Read the entire file
  while (append_vcf(vcf, vcf_fn, n_batch))
    ++n_batch;

  vcf.open(WRITE_MODE, output); // Change to write mode
  vcf.open_for_writing(copts.threads);
  print_log(log_severity::debug, __HERE__, " Read ", vcf.variants.size(), " variants.");

  std::vector<gyper::Vcf> next_vcfs(vcfs.size() - 1);
  n_batch = 1;

  // Open all VCFs and add sample names
  for (long i = 1; i < static_cast<long>(vcfs.size()); ++i)
  {
    auto const & vcf_fn = vcfs[i];
    assert((i - 1) < static_cast<long>(next_vcfs.size()));
    gyper::Vcf & next_vcf = next_vcfs[i - 1];
    load_vcf(next_vcf, vcf_fn, 0);

    assert(next_vcf.variants.size() == next_vcfs[0].variants.size()); // all bins have the same number of variants

    // Add the samples on the first batch
    vcf.sample_names.insert(vcf.sample_names.end(),
                            next_vcf.sample_names.begin(),
                            next_vcf.sample_names.end());
  }

  // Check whether the sample names are unique
  {
    std::vector<std::string> sample_names_cp(vcf.sample_names);
    std::sort(sample_names_cp.begin(), sample_names_cp.end());
    auto uniq_it = std::unique(sample_names_cp.begin(), sample_names_cp.end());

    if (uniq_it != sample_names_cp.end())
    {
      print_log(log_severity::warning,
                __HERE__,
                " Sample names are not unique. The output VCF file will contain duplicated sample names (which is "
                "against the VCF specs).");
    }
  }

  // We have all the samples, write the header now
  vcf.write_header(false); // do not drop genotypes
  long v_next{0}; // number of variants analyzed in previous batches

  // lambda function for updating the reach
  //auto update_reach =
  //  [&reach](std::vector<Variant> const & new_variants) -> void
  //  {
  //    for (auto const & v : new_variants)
  //    {
  //      long curr_reach = v.seqs[0].size() + v.abs_pos;
  //
  //      if (curr_reach > reach)
  //        reach = curr_reach;
  //    }
  //  };

  for (long v{0}; v < static_cast<long>(vcf.variants.size()); ++v)
  {
    auto & var = vcf.variants[v];

    // Trigger read next batch
    if (next_vcfs.size() > 0 && (v - v_next) == static_cast<long>(next_vcfs[0].variants.size()))
    {
      print_log(log_severity::debug, __HERE__, " Updating next_vcfs");

      v_next += static_cast<long>(next_vcfs[0].variants.size());

      for (long i = 1; i < static_cast<long>(vcfs.size()); ++i)
      {
        auto const & vcf_fn = vcfs[i];
        Vcf & next_vcf = next_vcfs[i - 1];
        next_vcf.variants.clear(); // clear previous batch
        bool ret = append_vcf(next_vcf, vcf_fn, n_batch);

        if (!ret)
        {
          print_log(log_severity::error, __HERE__, " Could not find file ", vcf_fn);
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
      var.stats.add_stats(next_vcf_var.stats);

      std::move(next_vcf_var.calls.begin(),
                next_vcf_var.calls.end(),
                std::back_inserter(var.calls));
    }

    // Add strand bias, this must happend before the INFO is generated
    //stats.add_stats(var.infos);

    if (var.calls.size() != vcf.sample_names.size())
    {
      print_log(log_severity::error, "Number of calls a variant had did not matches the number of samples "
                              , var.calls.size(), " vs. ", vcf.sample_names.size());
      std::exit(1);
    }
  }

  //vcf.close_vcf_file();
  //vcf.write_tbi_index();
}


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

  print_log(log_severity::debug, __HERE__, " Total number of samples read is "
                          , vcf.sample_names.size());

  vcf.write_header(); // Now that we know all the sample names we can write the header

  // For each variant in the first VCF, add the calls from the other VCFs
  for (auto & var : vcf.variants)
  {
    //VarStats stats(var.seqs.size());
    //stats.read_stats(var.infos);

    for (auto & next_vcf : next_vcfs)
    {
      if (!next_vcf.bgzf_in)
        continue;

      assert(next_vcf.variants.size() == 0);
      bool const SUCCESS = next_vcf.read_record();

      if (!SUCCESS)
      {
        print_log(log_severity::error, __HERE__, " There was a problem reading "
                                , next_vcf.filename, ".");
        std::exit(1);
      }

      assert(next_vcf.variants.size() == 1);
      var.stats.add_stats(next_vcf.variants[0].stats);

      std::move(next_vcf.variants[0].calls.begin(),
                next_vcf.variants[0].calls.end(),
                std::back_inserter(var.calls));

      next_vcf.variants.clear();
    }

    //stats.add_stats(var.infos);

    // generate the rest of the INFOs
    var.generate_infos();

    if (var.calls.size() != vcf.sample_names.size())
    {
      print_log(log_severity::error,
                __HERE__,
                " Number of calls a variant had did not matches the number of samples ",
                var.calls.size(),
                " vs. ",
                vcf.sample_names.size());
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
  vcf.write_tbi_index();
}


void
vcf_merge_and_filter(std::vector<std::string> const & vcfs,
                     std::string const & output,
                     std::map<std::pair<uint16_t, uint16_t>,
                              std::map<std::pair<uint16_t, uint16_t>, int8_t> > const & ph)
{
  if (vcfs.size() == 0)
  {
    print_log(log_severity::error, __HERE__, " No input VCFs");
    std::exit(1);
  }

  auto const & copts = *(Options::const_instance());
  auto const & vcf_fn = vcfs[0];
  long var_id{0};

  Vcf vcf;
  load_vcf(vcf, vcf_fn, 0);
  long n_batch{1};
  long v_next{0}; // number of variants analyzed in previous batches

  // Read the entire file
  while (append_vcf(vcf, vcf_fn, n_batch))
    ++n_batch;

  vcf.open(WRITE_MODE, output); // Change to write mode
  vcf.open_for_writing(copts.threads);
  print_log(log_severity::debug, __HERE__, " Read ", vcf.variants.size(), " variants.");

  // Add map from hap_id to variant index
  std::unordered_map<int32_t, long> hap_id2var_id;

  for (long v{0}; v < static_cast<long>(vcf.variants.size()); ++v)
  {
    Variant const & var = vcf.variants[v];
    assert(var.hap_id >= 0);
    hap_id2var_id.emplace(var.hap_id, var_id);
    var_id += (static_cast<long>(var.seqs.size()) - 1l);
  }

  var_id = 0;   // reset var_id
  assert(hap_id2var_id.size() == vcf.variants.size()); // no duplicate hap_ids

  std::vector<gyper::Vcf> next_vcfs(vcfs.size() - 1);
  n_batch = 1;

  // Open all VCFs and add sample names
  for (long i{1}; i < static_cast<long>(vcfs.size()); ++i)
  {
    auto const & vcf_fn = vcfs[i];
    assert((i - 1) < static_cast<long>(next_vcfs.size()));
    gyper::Vcf & next_vcf = next_vcfs[i - 1];
    load_vcf(next_vcf, vcf_fn, 0);

    assert(next_vcf.variants.size() == next_vcfs[0].variants.size()); // all bins have the same number of variants
  }

  vcf.sample_names.clear(); // this should only contain variant sites
  vcf.write_header(true); // is dropping genotypes

  for (long v{0}; v < static_cast<long>(vcf.variants.size()); ++v)
  {
    auto & var = vcf.variants[v];

    // Trigger read next batch
    if (next_vcfs.size() > 0 && (v - v_next) == static_cast<long>(next_vcfs[0].variants.size()))
    {
      print_log(log_severity::debug, __HERE__, " Updating next_vcfs");

      v_next += static_cast<long>(next_vcfs[0].variants.size());

      for (long i = 1; i < static_cast<long>(vcfs.size()); ++i)
      {
        auto const & vcf_fn = vcfs[i];
        Vcf & next_vcf = next_vcfs[i - 1];
        next_vcf.variants.clear(); // clear previous batch
        bool ret = append_vcf(next_vcf, vcf_fn, n_batch);

        if (!ret)
        {
          print_log(log_severity::error, __HERE__, " Could not find file ", vcf_fn);
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
      var.stats.add_stats(next_vcf_var.stats);

      std::move(next_vcf_var.calls.begin(),
                next_vcf_var.calls.end(),
                std::back_inserter(var.calls));
    }

    std::vector<int8_t> is_good_alt = var.generate_infos();

    // write out good alleles
    for (long a{0}; a < static_cast<long>(var.seqs.size()) - 1l; ++a)
    {
      ++var_id;

      if (is_good_alt[a] == 0)
      {
        //BOOST_LOG_TRIVIAL(info) << __HERE__ << " Skipping alt " << a << " of " << var.to_string();
        continue; // skip
      }

      Variant new_var;
      new_var.abs_pos = var.abs_pos;
      new_var.seqs.reserve(2);
      new_var.seqs.push_back(var.seqs[0]);
      new_var.seqs.push_back(var.seqs[a + 1]);
      new_var.infos["GT_ID"] = std::to_string(var_id);

      // Add INFO of anti events
      {
        std::ostringstream ss_anti;
        std::ostringstream ss_hap;

        bool is_anti_empty{true};
        bool is_hap_empty{true};

        // Other alleles in this variant are anti alleles
        for (long a2{a + 1}; a2 < static_cast<long>(var.seqs.size()) - 1l; ++a2)
        {
          if (is_good_alt[a2] == 0)
            continue;

          if (!is_anti_empty)
            ss_anti << ",";

          ss_anti << (var_id + a2 - a);
          is_anti_empty = false;
        }

        // Check ph
        // ph is std::map<std::pair<uint16_t, uint16_t>, std::map<std::pair<uint16_t, uint16_t>, int8_t> >
        auto find_it = ph.find(std::make_pair<uint16_t, uint16_t>(var.hap_id, a + 1));

        if (find_it != ph.end())
        {
          std::map<std::pair<uint16_t, uint16_t>, int8_t> const & other_map = find_it->second;

          for (auto const & other : other_map)
          {
            int32_t const other_hap_id = other.first.first;
            int32_t const other_allele = other.first.second;

            if (other_allele == 0)
              continue;

            int8_t const flags = other.second;

            if (flags != IS_ANY_HAP_SUPPORT && flags != IS_ANY_ANTI_HAP_SUPPORT)
              continue;

            auto other_find_it = hap_id2var_id.find(other_hap_id);
            assert(other_find_it != hap_id2var_id.end());
            long const var_id_other = other_find_it->second + other_allele;

            if (flags == IS_ANY_HAP_SUPPORT)
            {
              if (!is_hap_empty)
                ss_hap << ",";

              ss_hap << var_id_other;
              is_hap_empty = false;
            }
            else
            {
              assert(flags == IS_ANY_ANTI_HAP_SUPPORT);

              if (!is_anti_empty)
                ss_anti << ",";

              ss_anti << var_id_other;
              is_anti_empty = false;
            }
          }
        }

        if (!is_anti_empty)
          new_var.infos["GT_ANTI_HAPLOTYPE"] = ss_anti.str();

        if (!is_hap_empty)
          new_var.infos["GT_HAPLOTYPE"] = ss_hap.str();
      }

      vcf.write_record(new_var, "", false /*FILTER_ZERO_QUAL*/, true /*is dropping genotypes*/);
    }

    var = Variant(); // Free memory
  }

  // Write the remaining broken variants
  vcf.close_vcf_file();
  vcf.write_tbi_index();
}


void
vcf_merge_and_break(std::vector<std::string> const & vcfs,
                    std::string const & output,
                    std::string const & region,
                    bool const FILTER_ZERO_QUAL,
                    bool const force_no_variant_overlapping,
                    bool const force_no_break_down)
{
  if (vcfs.size() == 0)
    return;

  auto const & copts = *(Options::const_instance());
  long const ploidy = copts.ploidy;
  GenomicRegion genomic_region(region);
  uint32_t const region_begin = 1 + absolute_pos.get_absolute_position(genomic_region.chr,
                                                                       genomic_region.begin);

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
  vcf.open_for_writing(copts.threads);
  print_log(log_severity::debug, __HERE__, " Read ", vcf.variants.size(), " variants.");

  std::vector<gyper::Vcf> next_vcfs(vcfs.size() - 1);
  n_batch = 1;

  // Open all VCFs and add sample names
  for (long i = 1; i < static_cast<long>(vcfs.size()); ++i)
  {
    auto const & vcf_fn = vcfs[i];
    assert((i - 1) < static_cast<long>(next_vcfs.size()));
    gyper::Vcf & next_vcf = next_vcfs[i - 1];
    load_vcf(next_vcf, vcf_fn, 0);

    assert(next_vcf.variants.size() == next_vcfs[0].variants.size()); // all bins have the same number of variants

    // Add the samples on the first batch
    vcf.sample_names.insert(vcf.sample_names.end(),
                            next_vcf.sample_names.begin(),
                            next_vcf.sample_names.end());
  }

  // Check whether the sample names are unique
  {
    std::vector<std::string> sample_names_cp(vcf.sample_names);
    std::sort(sample_names_cp.begin(), sample_names_cp.end());
    auto uniq_it = std::unique(sample_names_cp.begin(), sample_names_cp.end());

    if (uniq_it != sample_names_cp.end())
    {
      print_log(log_severity::warning,
                __HERE__,
                " Sample names are not unique. The output VCF file will contain duplicated sample names (which is "
                "against the VCF specs).");
    }
  }

  // We have all the samples, write the header now
  vcf.write_header(copts.is_dropping_genotypes);

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

  for (long v{0}; v < static_cast<long>(vcf.variants.size()); ++v)
  {
    auto & var = vcf.variants[v];

    // Trigger read next batch
    if (next_vcfs.size() > 0 && (v - v_next) == static_cast<long>(next_vcfs[0].variants.size()))
    {
      print_log(log_severity::debug, __HERE__, " Updating next_vcfs");

      v_next += static_cast<long>(next_vcfs[0].variants.size());

      for (long i = 1; i < static_cast<long>(vcfs.size()); ++i)
      {
        auto const & vcf_fn = vcfs[i];
        Vcf & next_vcf = next_vcfs[i - 1];
        next_vcf.variants.clear(); // clear previous batch
        bool ret = append_vcf(next_vcf, vcf_fn, n_batch);

        if (!ret)
        {
          print_log(log_severity::error, __HERE__, " Could not find file ", vcf_fn);
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
      //next_vcf_var.scan_calls();
      var.stats.add_stats(next_vcf_var.stats);

      std::move(next_vcf_var.calls.begin(),
                next_vcf_var.calls.end(),
                std::back_inserter(var.calls));
    }

    // Add strand bias, this must happend before the INFO is generated
    //stats.add_stats(var.infos);

    if (var.calls.size() != vcf.sample_names.size())
    {
      print_log(log_severity::error, "Number of calls a variant had did not matches the number of samples "
                              , var.calls.size(), " vs. ", vcf.sample_names.size());
      std::exit(1);
    }

    //assert(var.infos.count("SBF1") == 1);
    //assert(var.stats.per_allele.size() == var.seqs.size());

    // break down the merged variants
    bool const is_no_variant_overlapping{copts.no_variant_overlapping ||
                                         force_no_variant_overlapping};

    bool const is_all_biallelic{copts.is_all_biallelic};
    std::vector<Variant> new_variants;

    if (force_no_break_down)
    {
      new_variants.push_back(std::move(var));
    }
    else
    {
      new_variants = break_down_variant(std::move(var),
                                        reach,
                                        is_no_variant_overlapping,
                                        is_all_biallelic);
    }

    assert(new_variants.size() > 0);

    //BOOST_LOG_TRIVIAL(info) << __HERE__ << " generating INFOs";

    for (auto it = new_variants.begin(); it != new_variants.end();)  // no increment & new_var : new_variants)
    {
      auto & new_var = *it;
      long const normalize_distance = new_var.normalize();

      if (normalize_distance <= 200)
      {
        if (ploidy > 2)
          new_var.update_camou_phred(ploidy);

        auto is_good_alt = new_var.generate_infos();

        if (!Options::const_instance()->force_no_filter_bad_alts &&
            std::all_of(is_good_alt.begin(), is_good_alt.end(), [](int8_t is_good){
            return is_good == 0;
          }))
        {
          print_log(log_severity::debug, __HERE__, " Removed variant=", new_var.to_string(true)
                                  , " because every alt was bad");
          it = new_variants.erase(it);
        }
        else
        {
          ++it;
        }

        // Remove uneeded infos
        //new_var.infos.erase("CRal");
        //new_var.infos.erase("MMal");
        //new_var.infos.erase("MQal");
        //new_var.infos.erase("MQsquared");
        //new_var.infos.erase("MQSal");
        //new_var.infos.erase("SDal");
      }
      else
      {
        print_log(log_severity::warning, __HERE__, " Removed a variant which moved ", normalize_distance
                                  , " bp during normalization.");
        it = new_variants.erase(it);
      }
    }

    //BOOST_LOG_TRIVIAL(info) << __HERE__ << " done generating INFOs";

    if (new_variants.size() > 0)
    {
      update_reach(new_variants);
      std::move(new_variants.begin(), new_variants.end(), std::back_inserter(broken_vars));

      long constexpr W{700}; // Print variants that are more than W bp before the newest one

      auto min_max_it_pair =
        std::minmax_element(broken_vars.begin(),
                            broken_vars.end(),
                            [](Variant const & a, Variant const & b)
        {
          return a.abs_pos < b.abs_pos;
        });

      assert(min_max_it_pair.first != broken_vars.end());
      assert(min_max_it_pair.second != broken_vars.end());
      long const min_abs_pos = min_max_it_pair.first->abs_pos;
      long const max_abs_pos = min_max_it_pair.second->abs_pos;
      assert(min_abs_pos <= max_abs_pos);

      if ((min_abs_pos + 2l * W) < max_abs_pos)
      {
        // Make sure we do no print outside of the region
        long const reg_end =
          std::min(static_cast<long>(region_end), max_abs_pos - static_cast<long>(W));

        if (reg_end >= static_cast<long>(region_begin))
        {
          vcf.write_records(region_begin,
                            static_cast<uint32_t>(reg_end),
                            FILTER_ZERO_QUAL,
                            copts.is_dropping_genotypes,
                            broken_vars);

          // Remove variants that were written (or before the region, since they will never be printed)
          broken_vars.erase(
            std::remove_if(broken_vars.begin(),
                           broken_vars.end(),
                           [&](Variant const & v)
            {
              return static_cast<long>(v.abs_pos) <= reg_end;
            }),
            broken_vars.end()
            );
        }
      }
    }

    var = Variant(); // Free memory
  }

  // Write the remaining broken variants
  vcf.write_records(region_begin,
                    region_end,
                    FILTER_ZERO_QUAL,
                    copts.is_dropping_genotypes,
                    broken_vars);

  vcf.close_vcf_file();
  vcf.write_tbi_index();
}


void
vcf_concatenate(std::vector<std::string> const & vcfs,
                std::string const & output,
                bool const SKIP_SORT,
                bool const SITES_ONLY,
                bool const WRITE_TBI,
                std::string const & region)
{
  print_log(log_severity::info, __HERE__, " running vcf_concatenate.");

  gyper::Vcf vcf;

  if (vcfs.size() == 0)
  {
    vcf.open(WRITE_MODE, output);
    vcf.open_for_writing();
    vcf.write_header();
    vcf.close_vcf_file();
    print_log(log_severity::warning, __HERE__, " nothing to do.");
    return;
  }

  if (SKIP_SORT)
  {
    print_log(log_severity::info, __HERE__, " running vcf_concatenate without sort.");
    vcf.open(WRITE_MODE, output);
    vcf.open_for_writing();

    if (SITES_ONLY)
    {
      vcf.sample_names.clear();

      for (auto & var : vcf.variants)
        var.calls.clear();
    }

    for (long i{0}; i < static_cast<long>(vcfs.size()); ++i)
    {
      // Skip if the filename contains '*'
      if (std::count(vcfs[i].begin(), vcfs[i].end(), '*') > 0)
      {
        print_log(log_severity::warning, __HERE__, " skipped VCF: ", vcfs[i]);
        continue;
      }

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
        print_log(log_severity::error,
                  "[graphtyper::vcf_operations] The VCF file ",
                  vcfs[i],
                  " has different amount of samples! (",
                  next_vcf.sample_names.size(),
                  " but not ",
                  vcf.sample_names.size(), ")");
        std::exit(1);
      }

      long n_read{0};

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
        ++n_read;
      }

      print_log(log_severity::info, __HERE__, " read ", n_read, " records.");
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

    print_log(log_severity::debug, __HERE__, " Total number of samples read is ", vcf.sample_names.size());

    for (long i{0}; i < static_cast<long>(vcfs.size()); ++i)
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
        print_log(log_severity::error,
                  __HERE__,
                  " The VCF file ",
                  vcfs[i],
                  " has unexpected number of sample names! (",
                  next_vcf.sample_names.size(),
                  " but not ",
                  vcf.sample_names.size(),
                  ")");
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

  if (WRITE_TBI)
    vcf.write_tbi_index();
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
  print_log(log_severity::debug, __HERE__, " Total number of samples read is "
                          , vcf_in.sample_names.size());

  // Copy sample names
  std::copy(vcf_in.sample_names.begin(),
            vcf_in.sample_names.end(),
            std::back_inserter(vcf_out.sample_names));

  GenomicRegion genomic_region(region);
  uint32_t const region_begin = 1 + absolute_pos.get_absolute_position(genomic_region.chr,
                                                                       genomic_region.begin);

  uint32_t const region_end = absolute_pos.get_absolute_position(genomic_region.chr,
                                                                 genomic_region.end);

  // Read first record
  bool not_at_end = vcf_in.read_record();

  // Write header to output
  vcf_out.write_header();
  long reach{-1}; // Indicates how long the previous variants reached

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
      print_log(log_severity::error,
                __HERE__,
                " The number of calls, ",
                vcf_in.variants[0].calls.size(),
                " did not match the number of samples, ",
                vcf_in.sample_names.size());

      std::exit(1);
    }

    bool const is_no_variant_overlapping{Options::const_instance()->no_variant_overlapping};
    bool const is_all_biallelic{Options::const_instance()->is_all_biallelic};
    assert(vcf_in.variants[0].infos.count("SBF1") == 1);
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
                            false, // is_dropping_genotypes
                            vcf_out.variants);

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

  vcf_out.write_records(region_begin,
                        region_end,
                        true /*FILTER_ZERO_QUAL*/,
                        false /*is_dropping_genotypes*/,
                        vcf_out.variants);

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
              std::back_inserter(vcf_out.variants));

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
