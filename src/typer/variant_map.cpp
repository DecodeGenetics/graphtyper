#include <cassert>
#include <cstdint> // uint32_t
#include <unordered_map>

#include <graphtyper/graph/graph.hpp>
#include <graphtyper/graph/reference_depth.hpp>
#include <graphtyper/graph/var_record.hpp>
#include <graphtyper/typer/variant.hpp>
#include <graphtyper/typer/variant_map.hpp>
#include <graphtyper/typer/variant_support.hpp>
#include <graphtyper/typer/vcf.hpp>
#include <graphtyper/utilities/options.hpp>
#include <graphtyper/utilities/type_conversions.hpp>


namespace gyper
{


void
VariantMap::set_pn_count(std::size_t const pn_count)
{
  map_mutexes = std::vector<std::mutex>(pn_count);
  varmaps.resize(pn_count);
}


void
VariantMap::add_variants(std::vector<VariantCandidate> && vars, std::size_t const pn_index)
{
  assert(pn_index < varmaps.size());
  assert(pn_index < map_mutexes.size());

  {
    std::lock_guard<std::mutex> lock(map_mutexes[pn_index]);
    auto & varmap = varmaps[pn_index];

    for (auto && var : vars)
    {
      assert(var.seqs.size() >= 2);
      assert(var.is_normalized());
      assert(var.seqs[0].size() > 0);
      assert(var.seqs[1].size() > 0);
      auto find_it = varmap.find(var);

      if (var.is_low_qual)
      {
        if (find_it == varmap.end())
          varmap[std::move(var)] = {0u, 1u}; // 0 hq, 1 lq support
        else if (find_it->second.second < 0xFFFFul)
          ++find_it->second.second; // Increment lq support
      }
      else
      {
        if (find_it == varmap.end())
          varmap[std::move(var)] = {1u, 0u}; // 1 hq, 0 lq
        else if (find_it->second.first < 0xFFFFul)
          ++find_it->second.first; // Increment HQ support
      }
    }
  }
}


void
VariantMap::create_varmap_for_all()
{
  assert(varmaps.size() == global_reference_depth.depths.size());

  for (std::size_t i = 0; i < varmaps.size(); ++i)
  {
    for (auto map_it = varmaps[i].begin(); map_it != varmaps[i].end(); ++map_it)
    {
      uint16_t const TOTAL_SUPPORT = std::min(static_cast<uint32_t>(0xFFFFul),
                                              static_cast<uint32_t>(map_it->second.first) + static_cast<uint32_t>(map_it->second.second / 2)
                                              );

      // Check if support is above cutoff
      if (TOTAL_SUPPORT  >= Options::instance()->minimum_variant_support)
      {
        VariantSupport new_var_support(TOTAL_SUPPORT, global_reference_depth.get_read_depth(map_it->first, i));

        // Check for underflow
        if (new_var_support.depth >= TOTAL_SUPPORT)
        {
          assert (static_cast<int>(new_var_support.depth) - static_cast<int>(map_it->second.second / 2) >= 0);
          new_var_support.depth -= map_it->second.second / 2;
        }

        // Check if ratio is above cutoff
        if (new_var_support.is_ratio_above_cutoff())
        {
          // In this case, we can add the variant
          assert(new_var_support.is_above_cutoff());
          pool_varmap[map_it->first].push_back(std::move(new_var_support));
        }
      }
    }
  }
}


void
VariantMap::filter_varmap_for_all()
{
  BOOST_LOG_TRIVIAL(info) << "[graphtyper::variant_map] Number of variants above minimum cutoff is " << pool_varmap.size();

  // No need to filter no variants, in fact it will segfault
  if (pool_varmap.size() == 0)
    return;

  /** Filter variants with high error rate */
  if (varmaps.size() >= 10) // Only consider error rate if the sample size is small
  {
    for (auto map_it = pool_varmap.begin(); map_it != pool_varmap.end(); /*no increment*/)
    {
      // If we are certain that someone has this variant, don't bother calculating the error rate.
      if (std::any_of(map_it->second.begin(), map_it->second.end(), [](VariantSupport const & var_support){return var_support.is_highly_certainly_real();}))
      {
        ++map_it;
        continue;
      }

      uint64_t support_of_samples_below_ratio_cutoff = 0u;
      uint64_t depth_of_samples_below_ratio_cutoff = 0u;

      // Get support and depth of samples below cutoff with at least one supported read
      {
        for (auto const & var_support : map_it->second)
        {
          if (!var_support.is_ratio_above_cutoff())
          {
            support_of_samples_below_ratio_cutoff += var_support.support;
            depth_of_samples_below_ratio_cutoff += var_support.depth;
          }
        }
      }

      // Get coverage of all samples which do not have the variant
      {
        std::vector<uint32_t> pn_indexes_below;

        for (uint32_t i = 0; i < varmaps.size(); ++i)
        {
          if (varmaps[i].count(map_it->first) == 0)
            pn_indexes_below.push_back(i);
        }

        depth_of_samples_below_ratio_cutoff += global_reference_depth.get_total_read_depth_of_samples(map_it->first, pn_indexes_below);
      }

      // If there is either no support (error rate = 0) or low depth, we keep this variant
      if (support_of_samples_below_ratio_cutoff == 0 || depth_of_samples_below_ratio_cutoff < 500)
      {
        ++map_it;
        continue;
      }

      double const ABHOM_ERROR = static_cast<double>(support_of_samples_below_ratio_cutoff) / static_cast<double>(depth_of_samples_below_ratio_cutoff);

      // Remove the variant if the allalic balance homozygous error rate is too high
      if (ABHOM_ERROR > Options::instance()->maximum_homozygous_allele_balance)
      {
        BOOST_LOG_TRIVIAL(info) << "[graphtyper::variant_map] The variant " << map_it->first.print() << " has a too high error rate to be added ("
                                << ABHOM_ERROR << ").";
        map_it = pool_varmap.erase(map_it);
      }
      else
      {
        ++map_it;
      }
    }
  }
  else
  {
    BOOST_LOG_TRIVIAL(info) << "[graphtyper::variant_map] Too few samples to filter based on ABHom.";
  }

  uint32_t const WINDOW_NOT_STARTED = 0xFFFFFFFFULL;

  /** Limit to the number of variants in a 100 bp window */
  if (pool_varmap.size() > 0)
  {
    BOOST_LOG_TRIVIAL(info) << "[graphtyper::variant_map] Soft cap of variants in 100 bp window is " << Options::instance()->soft_cap_of_variants_in_100_bp_window;
    BOOST_LOG_TRIVIAL(info) << "[graphtyper::variant_map] Hard cap of variants in 100 bp window is " << Options::instance()->hard_cap_of_variants_in_100_bp_window;
    std::size_t const original_variant_support_cutoff = Options::instance()->minimum_variant_support;
    std::size_t const original_variant_support_ratio_cutoff = Options::instance()->minimum_variant_support_ratio;
    std::size_t const original_max_variants_in_100_bp_window = Options::instance()->soft_cap_of_variants_in_100_bp_window;
    std::vector<uint32_t> positions;

    for (auto map_it = pool_varmap.begin(); map_it != pool_varmap.end(); ++map_it)
      positions.push_back(map_it->first.abs_pos);

    assert(positions.size() == pool_varmap.size());
    assert(positions.size() > 0);
    uint32_t a = 0;
    uint32_t b = 0;
    uint32_t window_start = WINDOW_NOT_STARTED;
    std::deque<uint32_t> dq;
    std::vector<std::pair<uint32_t, uint32_t> > indexes_to_filter;

    // Search for a window with too many variants
    while (b < positions.size())
    {
      dq.push_back(positions[b]);

      // Remove from the back
      while(dq.front() + 100 < dq.back())
      {
        dq.pop_front();
        ++a;
      }

      // Check if the deque is too big
      if (window_start == WINDOW_NOT_STARTED && dq.size() > original_max_variants_in_100_bp_window)
      {
        window_start = a;
      }
      else if (window_start != WINDOW_NOT_STARTED && dq.size() < original_max_variants_in_100_bp_window)
      {
        indexes_to_filter.push_back({window_start, b - window_start});
        window_start = WINDOW_NOT_STARTED;
      }

      ++b;
    }

    if (window_start != WINDOW_NOT_STARTED)
      indexes_to_filter.push_back({window_start, positions.size() - window_start});

    BOOST_LOG_TRIVIAL(info) << "[graphtyper::variant_map] Found " << indexes_to_filter.size() << " windows with too many variants.";
    BOOST_LOG_TRIVIAL(info) << "[graphtyper::variant_map] Total variant range is " << (positions.back() - positions.front());

    // Start on the last window so the indexes are correct
    for (auto it = indexes_to_filter.rbegin(); it != indexes_to_filter.rend(); ++it)
    {
      assert(it->first < positions.size());
      assert(it->first + it->second - 1 < positions.size());
      std::size_t max_variants_in_100_bp_window = (original_max_variants_in_100_bp_window *
                                                  std::max(static_cast<uint32_t>(125), positions[it->first + it->second - 1] - positions[it->first])
                                                  ) / 125;
      auto map_it = std::next(pool_varmap.begin(), it->first);
      uint32_t window_size = 0;
      uint32_t erased = 0;

      while(map_it != pool_varmap.end() && Options::instance()->minimum_variant_support < 15)
      {
        if (window_size == it->second)
        {
          if (it->second - erased <= max_variants_in_100_bp_window)
          {
            break; // We are done, we have erased enough
          }
          else
          {
            // We need to keep increasing the cutoffs
            Options::instance()->minimum_variant_support += 1;

            if (Options::instance()->minimum_variant_support_ratio < 0.95)
              Options::instance()->minimum_variant_support_ratio += 0.033;

            if (max_variants_in_100_bp_window < Options::instance()->hard_cap_of_variants_in_100_bp_window &&
                Options::instance()->minimum_variant_support % 2 == 0
                )
            {
              max_variants_in_100_bp_window += 1; // Also increase the number of variants allowed every now and then
            }

            window_size = erased;
            assert(it->first < pool_varmap.size());
            map_it = std::next(pool_varmap.begin(), it->first);

            if (map_it == pool_varmap.end())
              break;
          }
        }

        assert(map_it != pool_varmap.end());

        if (std::any_of(map_it->second.begin(), map_it->second.end(), [](VariantSupport const & var_support){return var_support.is_above_cutoff();}))
        {
          ++map_it;
        }
        else
        {
          ++erased;
          BOOST_LOG_TRIVIAL(info) << "[graphtyper::variant_map] Due to too high graph complexity I erased the variant: " << map_it->first.print()
                                  << ". It had support in " << map_it->second.size() << " samples.";
          map_it = pool_varmap.erase(map_it);
        }

        ++window_size;
      }

      Options::instance()->minimum_variant_support = original_variant_support_cutoff;
      Options::instance()->minimum_variant_support_ratio = original_variant_support_ratio_cutoff;
    }
  }


  /** Limit to the number of non-SNPs in a 100 bp window */
  if (pool_varmap.size() > 0)
  {
    BOOST_LOG_TRIVIAL(info) << "[graphtyper::variant_map] Soft cap of non-SNPs in 100 bp window is " << Options::instance()->soft_cap_of_non_snps_in_100_bp_window;
    BOOST_LOG_TRIVIAL(info) << "[graphtyper::variant_map] Hard cap of non-SNPs in 100 bp window is " << Options::instance()->hard_cap_of_non_snps_in_100_bp_window;
    std::size_t const original_variant_support_cutoff = Options::instance()->minimum_variant_support;
    std::size_t const original_variant_support_ratio_cutoff = Options::instance()->minimum_variant_support_ratio;
    std::size_t const original_max_nonsnps_in_100_bp_window = Options::instance()->soft_cap_of_non_snps_in_100_bp_window;
    std::vector<uint32_t> positions;

    for (auto map_it = pool_varmap.begin(); map_it != pool_varmap.end(); ++map_it)
    {
      if (!map_it->first.is_snp_or_snps())
        positions.push_back(map_it->first.abs_pos);
      else
        positions.push_back(0);
    }

    if (positions.size() > original_max_nonsnps_in_100_bp_window)
    {
      uint32_t b = 0;
      uint32_t window_start = WINDOW_NOT_STARTED;
      std::deque<uint32_t> aq;
      std::deque<uint32_t> dq;
      std::vector<std::pair<uint32_t, uint32_t> > indexes_to_filter;

      // Search for a window with too many variants
      while (b < positions.size())
      {
        if (positions[b] == 0)
        {
          ++b;
          continue;
        }

        aq.push_back(b);
        dq.push_back(positions[b]);

        // Remove from the back
        while(dq.front() + 100 < dq.back())
        {
          dq.pop_front();
          aq.pop_front();
        }

        assert (dq.size() == aq.size());

        // Check if the deque is too big
        if (window_start == WINDOW_NOT_STARTED && dq.size() > original_max_nonsnps_in_100_bp_window)
        {
          window_start = aq.front();
        }
        else if (window_start != WINDOW_NOT_STARTED && dq.size() < original_max_nonsnps_in_100_bp_window)
        {
          indexes_to_filter.push_back({window_start, aq.back() + 1 - window_start});
          window_start = WINDOW_NOT_STARTED;
        }

        ++b;
      }

      if (window_start != WINDOW_NOT_STARTED)
      {
        indexes_to_filter.push_back({window_start, aq.back() + 1 - window_start});
      }

      BOOST_LOG_TRIVIAL(info) << "[graphtyper::variant_map] Found " << indexes_to_filter.size() << " windows with too many non SNPs.";
      assert (positions.size() > 0);

      // Start on the last window so the indexes are correct
      for (auto it = indexes_to_filter.rbegin(); it != indexes_to_filter.rend(); ++it)
      {
        assert (positions[it->first] != 0);
        assert (positions[it->first + it->second - 1] != 0);
        std::size_t max_nonsnps_in_100_bp_window = (original_max_nonsnps_in_100_bp_window *
                                                    std::max(static_cast<uint32_t>(125), positions[it->first + it->second - 1] - positions[it->first])
                                                    ) / 125;
        assert(it->first < positions.size());
        assert(it->first + it->second - 1 < positions.size());
        auto map_it = std::next(pool_varmap.begin(), it->first);
        uint32_t window_size = 0;
        uint32_t erased = 0;

        while(map_it != pool_varmap.end() && Options::instance()->minimum_variant_support < 15)
        {
          // Loop until we have looped the entire window
          if (window_size == it->second)
          {
            auto count_map_it = std::next(pool_varmap.begin(), it->first);

            // Count non SNPs
            if (std::count_if(count_map_it,
                              std::next(count_map_it, it->second),
                              [](std::pair<Variant, std::vector<VariantSupport> > const & var)
                              {
                                return var.first.is_snp_or_snps() == false;
                              }
                              ) - erased <= static_cast<long>(max_nonsnps_in_100_bp_window)
               )
            {
              break; // We are done, we have erased enough
            }
            else
            {
              // We need to keep increasing the cutoffs
              Options::instance()->minimum_variant_support += 1;

              if (Options::instance()->minimum_variant_support_ratio < 0.95)
                Options::instance()->minimum_variant_support_ratio += 0.033;

              if (max_nonsnps_in_100_bp_window < Options::instance()->hard_cap_of_non_snps_in_100_bp_window &&
                  Options::instance()->minimum_variant_support % 2 == 0
                  )
              {
                max_nonsnps_in_100_bp_window += 1; // Also increase the number of variants allowed every now and then (but to some cap)
              }

              window_size = erased;
              map_it = std::next(pool_varmap.begin(), it->first);

              if (map_it == pool_varmap.end())
                break;
            }
          }

          assert(map_it != pool_varmap.end());

          if (map_it->first.is_snp_or_snps())
          {
            // std::cout << "kept snp " << map_it->first.abs_pos << std::endl;
            ++map_it;
            ++window_size;
            continue;
          }

          if (std::any_of(map_it->second.begin(), map_it->second.end(), [](VariantSupport const & var_support){return var_support.is_above_cutoff();}))
          {
            // std::cout << "kept non-snp " << map_it->first.abs_pos << std::endl;
            ++map_it;
          }
          else
          {
            ++erased;
            BOOST_LOG_TRIVIAL(info) << "[graphtyper::variant_map] Due to too high non-SNP graph complexity I erased the variant: " << map_it->first.print()
                                    << ". It had support in " << map_it->second.size() << " samples with cutoffs " << Options::instance()->minimum_variant_support << " " << Options::instance()->minimum_variant_support_ratio;
            map_it = pool_varmap.erase(map_it);
          }

          ++window_size;
        }

        Options::instance()->minimum_variant_support = original_variant_support_cutoff;
        Options::instance()->minimum_variant_support_ratio = original_variant_support_ratio_cutoff;
      }
    } // if (positions.size() > 0)
  }

  /** Break down variants to check if they are duplicates */
  std::vector<VariantCandidate> broken_vars_to_add;
  std::vector<std::vector<VariantSupport> > supports_to_add;

  for (auto map_it = pool_varmap.begin(); map_it != pool_varmap.end(); /*no increment*/)
  {
    VariantCandidate var(map_it->first);
    std::size_t const EXTRA_BASES_TO_ADD = 5;

    for (std::size_t i = 0; i < EXTRA_BASES_TO_ADD; ++i)
      if (!var.add_base_in_front(false /*add_N*/))
        break;

    for (std::size_t i = 0; i < EXTRA_BASES_TO_ADD; ++i)
      if (!var.add_base_in_back(false /*add_N*/))
        break;

    std::size_t const THRESHOLD = 1;

    //
    Variant var_cp;
    var_cp.abs_pos = var.abs_pos;
    var_cp.seqs = var.seqs;
    std::vector<Variant> new_broken_down_vars = break_down_variant(Variant(var_cp), THRESHOLD);
    assert(new_broken_down_vars.size() != 0);

    if (new_broken_down_vars.size() == 1)
    {
      ++map_it; // The variant could not be broken down, move on
      continue;
    }

    for (auto & broken_var : new_broken_down_vars)
      broken_var.normalize();

    // Change Variant -> VariantCandidate
    std::vector<VariantCandidate> new_broken_down_var_candidates(new_broken_down_vars.size());

    for (unsigned i = 0; i < new_broken_down_vars.size(); ++i)
    {
      new_broken_down_var_candidates[i].abs_pos = new_broken_down_vars[i].abs_pos;
      new_broken_down_var_candidates[i].seqs = std::move(new_broken_down_vars[i].seqs);
      new_broken_down_var_candidates[i].is_low_qual = var.is_low_qual;
    }

    // Remove broken variants which are already in the varmap
    new_broken_down_var_candidates.erase(
      std::remove_if(new_broken_down_var_candidates.begin(),
                     new_broken_down_var_candidates.end(),
                     [&](VariantCandidate const & broken_var){
        return pool_varmap.find(broken_var) != pool_varmap.end();
      }), new_broken_down_var_candidates.end());

    // Remove broken variants that are already in the graph
    new_broken_down_var_candidates.erase(
      std::remove_if(new_broken_down_var_candidates.begin(),
                     new_broken_down_var_candidates.end(),
                     [&](VariantCandidate const & broken_var){
        return graph.is_variant_in_graph(Variant(broken_var));
      }), new_broken_down_var_candidates.end());

    if (new_broken_down_var_candidates.size() <= 1)
    {
      // Add this variant (if there is one) and remove the old one
      if (new_broken_down_var_candidates.size() == 1)
      {
        broken_vars_to_add.push_back(std::move(new_broken_down_var_candidates[0]));
        supports_to_add.push_back(map_it->second);
      }

      map_it = pool_varmap.erase(map_it); // Erase the old variant
    }
    else
    {
      ++map_it; // Don't erase the old variant, just move on to the next one
    }
  }

  assert(broken_vars_to_add.size() == supports_to_add.size());

  for (std::size_t i = 0; i < broken_vars_to_add.size(); ++i)
    pool_varmap[broken_vars_to_add[i]] = supports_to_add[i];
}


void
VariantMap::write_vcf(std::string const & output_name)
{
  Vcf new_variant_vcf(WRITE_BGZF_MODE, output_name);

  for (auto map_it = pool_varmap.begin(); map_it != pool_varmap.end(); ++map_it)
    new_variant_vcf.variants.push_back(map_it->first);

  new_variant_vcf.write();
}


VariantMap global_varmap;

} // namespace gyper
