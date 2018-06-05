#include <cassert>
#include <cstdint> // uint32_t
#include <ios>
#include <iomanip>
#include <unordered_map>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/unordered_map.hpp>

#include <graphtyper/graph/absolute_position.hpp>
#include <graphtyper/graph/graph.hpp>
#include <graphtyper/graph/reference_depth.hpp>
#include <graphtyper/graph/var_record.hpp>
#include <graphtyper/typer/variant.hpp>
#include <graphtyper/typer/variant_map.hpp>
#include <graphtyper/typer/variant_support.hpp>
#include <graphtyper/typer/vcf.hpp>
#include <graphtyper/utilities/io.hpp>
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

  std::lock_guard<std::mutex> lock(map_mutexes[pn_index]);
  auto & varmap = varmaps[pn_index];

  for (auto && var : vars)
  {
    assert(var.seqs.size() >= 2);
    assert(var.is_normalized());
    assert(var.seqs[0].size() > 0);
    assert(var.seqs[1].size() > 0);

    // Extract all required information
    bool const IS_LOW_QUAL = var.is_low_qual;
    bool const IS_IN_PROPER_PAIR = var.is_in_proper_pair;
    bool const IS_MAPQ0 = var.is_mapq0;
    bool const IS_UNALIGNED = var.is_unaligned;
    bool const IS_CLIPPED = var.is_clipped;
    bool const IS_FIRST_IN_PAIR = var.is_first_in_pair;
    bool const IS_SEQ_REVERSED = var.is_seq_reversed;
    uint32_t const ORIGINAL_POS = var.original_pos;

    auto && new_item = std::make_pair<VariantCandidate, VariantSupport>(std::move(var), VariantSupport());
    auto it = varmap.insert(std::forward<std::pair<VariantCandidate, VariantSupport> >(new_item)).first;

    // Add the discovered variant candidate
    if (IS_LOW_QUAL)
      ++it->second.lq_support; // Increment LQ support
    else
      ++it->second.hq_support; // Increment HQ support

    ++it->second.depth; // Increment depth

    if (IS_IN_PROPER_PAIR)
      ++it->second.proper_pairs; // Increment number of proper pairs

    if (IS_MAPQ0)
      ++it->second.mapq0; // Increment number of reads with MAPQ==0

    if (IS_UNALIGNED)
      ++it->second.unaligned; // Increment number of unaligned reads

    if (IS_CLIPPED)
      ++it->second.clipped;

    if (IS_FIRST_IN_PAIR)
      ++it->second.first_in_pairs;

    if (IS_SEQ_REVERSED)
      ++it->second.sequence_reversed;

    it->second.unique_positions.insert(ORIGINAL_POS);
  }
}


void
VariantMap::create_varmap_for_all()
{
  assert(varmaps.size() == global_reference_depth.depths.size());
  unsigned const NUM_SAMPLES = varmaps.size();

  for (unsigned i = 0; i < NUM_SAMPLES; ++i)
  {
    auto & varmap = varmaps[i];

    for (auto map_it = varmap.begin(); map_it != varmap.end(); ++map_it)
    {
      VariantSupport & new_var_support = map_it->second;

      // Check if support is above cutoff and some reasonable hard filters
      if (new_var_support.is_support_above_cutoff())
      {
        new_var_support.set_depth(global_reference_depth.get_read_depth(map_it->first, i));

        // Check if ratio is above cutoff
        if (new_var_support.is_ratio_above_cutoff())
        {
          // In this case, we can add the variant
          new_var_support.pn_index = i;
          pool_varmap[map_it->first].push_back(new_var_support);
        }
      }
    }
  }
}


void
VariantMap::filter_varmap_for_all()
{
  BOOST_LOG_TRIVIAL(info) << "[graphtyper::variant_map] Number of variants above minimum cutoff is "
                          << pool_varmap.size();

  // No need to filter no variants, in fact it will segfault
  if (pool_varmap.size() == 0)
    return;

  uint32_t const WINDOW_NOT_STARTED = 0xFFFFFFFFULL;

  /** Limit to the number of variants in a 100 bp window */
  if (pool_varmap.size() > 0)
  {
    BOOST_LOG_TRIVIAL(info) << "[graphtyper::variant_map] Soft cap of variants in 100 bp window is "
                            << Options::instance()->soft_cap_of_variants_in_100_bp_window;
    BOOST_LOG_TRIVIAL(info) << "[graphtyper::variant_map] Hard cap of variants in 100 bp window is "
                            << Options::instance()->hard_cap_of_variants_in_100_bp_window;

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

  if (samples.size() == 0 || Options::instance()->stats.size() == 0)
  {
    for (auto map_it = pool_varmap.begin(); map_it != pool_varmap.end(); ++map_it)
      new_variant_vcf.variants.push_back(Variant(map_it->first));
  }
  else
  {
    auto const & pn = samples[0];

    // Determine filename
    std::string const discovery_fn =
      Options::instance()->stats + "/" + pn + "_discovery.tsv.gz";

    // Open the stream
    std::stringstream discovery_ss;

    // Write header to the stream
    discovery_ss << "CHROM\tPOS\tREF\tALT\tFORMAT";

    for (auto const & sample_name : samples)
      discovery_ss << "\t" << sample_name;

    discovery_ss << '\n';

    for (auto map_it = pool_varmap.begin(); map_it != pool_varmap.end(); ++map_it)
    {
      Variant var(map_it->first);
      assert (var.seqs.size() == 2);
      auto contig_pos = absolute_pos.get_contig_position(var.abs_pos);
      discovery_ss << contig_pos.first << "\t" << contig_pos.second << "\t";

      // REF
      std::copy(var.seqs[0].begin(), var.seqs[0].end(), std::ostream_iterator<char>(discovery_ss));
      discovery_ss << '\t';

      // ALT
      std::copy(var.seqs[1].begin(), var.seqs[1].end(), std::ostream_iterator<char>(discovery_ss));

      // FORMAT
      discovery_ss << "\tnaiveGT:AD:SeqDepth:HQsupport:LQsupport:ProperPairs:"
                      "MAPQ0:Unaligned:Clipped:"
                      "readUniquePositions:ratio:naivePL\t";

      // Add naive genotype calls
      std::sort(map_it->second.begin(),
                map_it->second.end(),
                [](VariantSupport const & a, VariantSupport const & b) -> bool
                {
                  return a.pn_index < b.pn_index;
                }
        );

      assert (Options::instance()->stats.size() > 0);
      std::size_t m = 0; // Index in map->second

      // Loop over all samples
      for (std::size_t p = 0; p < varmaps.size(); ++p)
      {
        std::vector<uint16_t> coverage(2);
        std::vector<uint16_t> phred(3);

        if (m >= map_it->second.size() || p < map_it->second[m].pn_index)
        {
          // Add empty call
          discovery_ss << ".\t";
        }
        else
        {
          auto const & var_support = map_it->second[m];
          assert (p == var_support.pn_index);

          coverage[1] = var_support.lq_support + var_support.hq_support;
          coverage[0] = var_support.depth - coverage[1];
          uint32_t phred0 = coverage[1] * 25u;
          uint32_t phred1 = var_support.depth * 3u;
          uint32_t phred2 = coverage[0] * 25u;
          uint32_t const normalizer = std::min(phred0, std::min(phred1, phred2));
          phred0 -= normalizer;
          phred1 -= normalizer;
          phred2 -= normalizer;
          assert (phred0 == 0 || phred1 == 0 || phred2 == 0);
          phred[0] = phred0 < 255 ? phred0 : 255;
          phred[1] = phred1 < 255 ? phred1 : 255;
          phred[2] = phred2 < 255 ? phred2 : 255;

          long const gt = std::distance(phred.begin(),
                                               std::find(phred.begin(),phred.end(), 0u)
            );

          if (gt == 0)
            discovery_ss << "0/0:";
          else if (gt == 1)
            discovery_ss << "0/1:";
          else
            discovery_ss << "1/1:";

          discovery_ss.setf(std::ios::fixed);
          discovery_ss.precision(2);
          discovery_ss << coverage[0] << ',' << coverage[1] << ':'
                       << var_support.depth << ':'
                       << var_support.hq_support << ':'
                       << var_support.lq_support << ':'
                       << var_support.proper_pairs << ':'
                       << var_support.mapq0 << ':'
                       << var_support.unaligned << ':'
                       << var_support.clipped << ':'
                       << var_support.unique_positions.size() << ':'
                       << var_support.get_ratio() << ':'
                       << phred[0] << ',' << phred[1] << ',' << phred[2]
                       << '\t';

          ++m;
        }
      }

      discovery_ss << '\n';
      new_variant_vcf.variants.push_back(std::move(var));
    }

    write_gzipped_to_file(discovery_ss, discovery_fn, false /*append*/);
  }

  new_variant_vcf.write();
}


void
save_variant_map(std::string const & path, VariantMap const & map)
{
  std::ofstream ofs(path.c_str(), std::ios::binary);

  if (!ofs.is_open())
  {
    BOOST_LOG_TRIVIAL(fatal) << "[graphtyper::variant_map] Could not save graph at '"
                             << path
                             << "'";

    std::exit(1);
  }

  boost::archive::binary_oarchive oa(ofs);
  oa << map;
}


VariantMap
load_variant_map(std::string const & path)
{
  VariantMap map;
  std::ifstream ifs(path.c_str(), std::ios::binary);

  if (!ifs.is_open())
  {
    BOOST_LOG_TRIVIAL(fatal) << "[graphtyper::constructor] Could not load graph at '"
                             << path
                             << "'";

    std::exit(1);
  }

  boost::archive::binary_iarchive ia(ifs);
  ia >> map;
  return map;
}


void
VariantMap::set_samples(std::vector<std::string> const & new_samples)
{
  samples = new_samples;
  this->set_pn_count(samples.size());
}


void
VariantMap::load_many_variant_maps(std::string const & path)
{
  this->samples.clear();
  this->pool_varmap.clear();  // Should be empty initially

  std::ifstream ifs(path.c_str());

  for (std::string line; std::getline(ifs, line);)
  {
    VariantMap new_map = load_variant_map(line);

    // Sample index offset
    long sample_index_offset = samples.size();

    // Copy samples
    std::copy(new_map.samples.begin(), new_map.samples.end(), std::back_inserter(samples));

    for (auto it = new_map.pool_varmap.begin(); it != new_map.pool_varmap.end(); ++it)
    {
      // Change all pn_indexes with offset
      for (auto variant_support : it->second)
        variant_support.pn_index += sample_index_offset;

      auto find_it = pool_varmap.find(it->first);

      if (find_it == pool_varmap.end())
      {
        // Not found
        pool_varmap[it->first] = it->second;
      }
      else
      {
        // Found
        std::copy(it->second.begin(), it->second.end(), std::back_inserter(find_it->second));
      }
    }
  }
}


VariantMap global_varmap;

} // namespace gyper
