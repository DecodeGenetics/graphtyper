#include <cassert>
#include <cstdint> // uint32_t
#include <cstdlib> // abs
#include <iomanip>
#include <ios>
#include <unordered_map>
#include <utility>

#include <cereal/archives/binary.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>

#include <graphtyper/graph/absolute_position.hpp>
#include <graphtyper/graph/graph.hpp>
#include <graphtyper/graph/reference_depth.hpp>
#include <graphtyper/graph/var_record.hpp>
#include <graphtyper/typer/variant.hpp>
#include <graphtyper/typer/variant_map.hpp>
#include <graphtyper/typer/variant_support.hpp>
#include <graphtyper/typer/vcf.hpp>
#include <graphtyper/utilities/io.hpp>
#include <graphtyper/utilities/logging.hpp>
#include <graphtyper/utilities/options.hpp>
#include <graphtyper/utilities/type_conversions.hpp>

namespace
{
std::pair<double, long> get_abs_dev_strand_bias(std::vector<gyper::VariantSupport> const & supports)
{
  long total_depth = 0;
  long total_reversed = 0;

  for (auto const & sup : supports)
  {
    total_depth += sup.lq_support + sup.hq_support;
    total_reversed += sup.sequence_reversed;
  }

  assert(total_depth > 0);
  double abs_dev_sb = std::abs(static_cast<double>(total_reversed) / static_cast<double>(total_depth) - 0.5);
  return std::pair<double, long>(abs_dev_sb, total_depth);
}

std::pair<double, long> get_abs_dev_read_bias(std::vector<gyper::VariantSupport> const & supports)
{
  long total_depth = 0;
  long total_first_read = 0;

  for (auto const & sup : supports)
  {
    total_depth += sup.lq_support + sup.hq_support;
    total_first_read += sup.first_in_pairs;
  }

  double abs_dev_rb = std::abs(static_cast<double>(total_first_read) / static_cast<double>(total_depth) - 0.5);
  return std::pair<double, long>(abs_dev_rb, total_depth);
}

} // namespace

namespace gyper
{
void VariantMap::add_variants(std::vector<VariantCandidate> && vars, long const sample_index)
{
  assert(varmaps.size() > 0);
  assert(sample_index < static_cast<long>(varmaps.size()));
  auto & varmap = varmaps[sample_index];

  for (auto & var : vars)
  {
    assert(var.seqs.size() >= 2);
    assert(var.is_normalized());
    assert(var.seqs[0].size() > 0);
    assert(var.seqs[1].size() > 0);

    // Extract all required information
    uint32_t const ORIGINAL_POS = var.original_pos;

    auto it = varmap.find(var);

    if (it == varmap.end())
    {
      it =
        varmap.insert(std::make_pair<VariantCandidate, VariantSupport>(VariantCandidate(var), VariantSupport())).first;

      // Expand to learn the true size
      it->second.is_indel = var.seqs[0].size() != var.seqs[1].size();
      long const old_size = std::max(var.seqs[0].size(), var.seqs[1].size()) - 1;
      var.expanded_normalized();
      assert(var.seqs[0].size() > 0);
      assert(var.seqs[1].size() > 0);
      it->second.var_size = std::max(var.seqs[0].size(), var.seqs[1].size()) - 1;
      it->second.growth = std::max(0l, static_cast<long>(it->second.var_size) - old_size);
    }

    ++it->second.depth; // Increment depth
    it->second.lq_support += ((var.flags & IS_LOW_BASE_QUAL) != 0);
    it->second.hq_support += ((var.flags & IS_LOW_BASE_QUAL) == 0);
    it->second.proper_pairs += ((var.flags & IS_PROPER_PAIR) != 0);

    if ((var.flags & IS_MAPQ_BAD) == 0)
      it->second.is_any_mapq_good = true;

    it->second.first_in_pairs += ((var.flags & IS_FIRST_IN_PAIR) != 0);
    it->second.sequence_reversed += ((var.flags & IS_SEQ_REVERSED) != 0);
    it->second.clipped += ((var.flags & IS_CLIPPED) != 0);
    it->second.unique_positions.insert(ORIGINAL_POS);
  }
}

void VariantMap::create_varmap_for_all(ReferenceDepth const & reference_depth)
{
  assert(varmaps.size() == reference_depth.depths.size());
  long const NUM_SAMPLES = static_cast<long>(varmaps.size());

  for (long i = 0; i < NUM_SAMPLES; ++i)
  {
    auto & varmap = varmaps[i];
    long new_min_support = minimum_variant_support;

    if (varmap.size() > 50)
    {
      while (new_min_support < 15)
      {
        long transitions = 0;
        long transversions = 0;
        long num_above_threshold = 0;

        for (auto it = varmap.cbegin(); it != varmap.cend(); ++it)
        {
          if (it->second.is_support_above_cutoff(new_min_support))
          {
            ++num_above_threshold;
            int ret = it->first.is_transition_or_transversion();

            if (ret == 1)
              ++transitions;
            else if (ret == 2)
              ++transversions;
          }
        }

        if ((num_above_threshold > 50 && transversions > transitions) ||
            (num_above_threshold > 75 &&
             (static_cast<double>(transitions) / static_cast<double>(transversions)) < 1.2) ||
            (num_above_threshold > 150 &&
             (static_cast<double>(transitions) / static_cast<double>(transversions)) < 1.3))
        {
          ++new_min_support;
        }
        else
        {
          break;
        }
      } // while(true)

#ifndef NDEBUG
      if (new_min_support > minimum_variant_support)
        print_log(log_severity::debug, "[graphtyper::variant_map] Min support increased to ", new_min_support);
#endif // NDEBUG
    }  // if (varmap.size() > 50)

    double const min_ratio = minimum_variant_support_ratio;

    for (auto map_it = varmap.begin(); map_it != varmap.end(); ++map_it)
    {
      VariantCandidate const & var = map_it->first;
      VariantSupport & new_var_support = map_it->second;

      // Check if support is above cutoff and some reasonable hard filters
      if (new_var_support.is_support_above_cutoff(new_min_support))
      {
        new_var_support.set_depth(reference_depth.get_read_depth(map_it->first, i));

        // Check if ratio is above cutoff
        if (new_var_support.is_ratio_above_cutoff(min_ratio))
        {
// In this case, we can add the variant
#ifndef NDEBUG
          new_var_support.pn_index = static_cast<uint32_t>(i);
#endif // NDEBUG

          pool_varmap[var].push_back(new_var_support);
        }
      }
    }
  } // for (long i = 0; i < NUM_SAMPLES; ++i)

#ifndef NDEBUG
  print_log(log_severity::debug, "Pool varmap size ", pool_varmap.size());

  for (auto it = pool_varmap.begin(); it != pool_varmap.end(); ++it)
  {
    print_log(log_severity::debug, it->first.print());
  }
#endif // NDEBUG
}

void VariantMap::filter_varmap_for_all()
{
  print_log(log_severity::debug,
            "[graphtyper::variant_map] Number of variants above minimum cutoff is ",
            pool_varmap.size());

  // No need to filter no variants, in fact it will segfault
  if (pool_varmap.size() == 0)
    return;

  // Filter on strand bias
  //*
  for (auto it = pool_varmap.begin(); it != pool_varmap.end();) // no increment
  {
    // Skip read bias on indels
    bool is_any_hq = false;

    for (auto const & sup : it->second)
    {
      if (sup.hq_support >= 5 && sup.proper_pairs >= 5)
      {
        is_any_hq = true;
        break;
      }
    }

    assert(it->second.size() > 0);
    bool const is_indel = it->second[0].is_indel;

    // Filter on strand bias
    if (!is_any_hq)
    {
      std::pair<double, long> abs_dev_sb_depth = get_abs_dev_strand_bias(it->second);
      double abs_dev_sb = abs_dev_sb_depth.first;
      long depth = abs_dev_sb_depth.second;

      if (is_indel && abs_dev_sb > 0.07)
        abs_dev_sb -= 0.07;

      if ((abs_dev_sb > 0.49999) || (abs_dev_sb > 0.45 && depth > 30) || (abs_dev_sb > 0.40 && depth > 80) ||
          (abs_dev_sb > 0.37 && depth > 200) || (abs_dev_sb > 0.34 && depth > 500))
      {
#ifndef NDEBUG
        print_log(log_severity::debug, "Strand bias removed ", it->first.print(), " ", abs_dev_sb, " ", depth);
#endif // NDEBUG
        it = pool_varmap.erase(it);
        continue;
      }

      // Filter on read bias
      if (!is_indel)
      {
        std::pair<double, long> abs_dev_rb_depth = get_abs_dev_read_bias(it->second);
        double abs_dev_rb = abs_dev_rb_depth.first;
        long depth = abs_dev_rb_depth.second;

        if ((abs_dev_rb > 0.49999 && depth > 10) || (abs_dev_rb > 0.45 && depth > 40) ||
            (abs_dev_rb > 0.40 && depth > 100) || (abs_dev_rb > 0.35 && depth > 500))
        {
#ifndef NDEBUG
          print_log(log_severity::debug, "Read bias removed ", it->first.print(), " ", abs_dev_rb, " ", depth);
#endif // NDEBUG
          it = pool_varmap.erase(it);
          continue;
        }
      }
    }

    ++it;
  }
  //*/

  /// Limit to the number of variants in a 100 bp window
  if (static_cast<long>(pool_varmap.size()) > Options::const_instance()->soft_cap_of_variants_in_100_bp_window)
  {
    long const max_variants_in_100bp_window = Options::const_instance()->soft_cap_of_variants_in_100_bp_window;
    assert(max_variants_in_100bp_window > 0);

    print_log(log_severity::debug,
              "[graphtyper::variant_map] Soft cap of variants in 100 bp window is ",
              Options::const_instance()->soft_cap_of_variants_in_100_bp_window);

    // Gather how many variants fall into each bucket (which covers 100 bp)
    std::vector<long> max_scores_in_bucket;
    auto window_it = pool_varmap.begin();
    long current_bucket = window_it->first.abs_pos / 100l;

    auto filter_window = [&max_scores_in_bucket, max_variants_in_100bp_window, this](PoolVarMap::iterator window_it)
    {
      print_log(log_severity::debug, "[graphtyper::variant_map] Too many variants! ", max_scores_in_bucket.size());

      // Find the minimum max score for passing a variant in this region
      std::vector<long> sorted_max_scores_in_bucket(max_scores_in_bucket);
      std::sort(sorted_max_scores_in_bucket.begin(), sorted_max_scores_in_bucket.end());
      long const index_with_min_score = sorted_max_scores_in_bucket.size() - max_variants_in_100bp_window;
      assert(index_with_min_score <= static_cast<long>(sorted_max_scores_in_bucket.size()));
      long const min_score_pass = std::min(50l, sorted_max_scores_in_bucket[index_with_min_score]);

      for (long s{0}; s < static_cast<long>(max_scores_in_bucket.size()); ++s)
      {
        if (max_scores_in_bucket[s] >= min_score_pass)
        {
          ++window_it;
        }
        else
        {
#ifndef NDEBUG
          print_log(log_severity::debug,
                    "Due to too high graph complexity I erased the variant: ",
                    window_it->first.print(),
                    ". It had support in ",
                    window_it->second.size(),
                    " samples and score of ",
                    max_scores_in_bucket[s]);
#endif // NDEBUG
          window_it = pool_varmap.erase(window_it);
        }
      }

      return window_it; // return the new iterator
    };

    for (auto map_it = pool_varmap.begin(); map_it != pool_varmap.end(); ++map_it)
    {
      long const bucket = map_it->first.abs_pos / 100l;
      assert(bucket >= current_bucket); // We should never go backwards in order

      if (bucket != current_bucket) // Let's continue looping through the variants until we reach the next bucket
      {
        assert(bucket > current_bucket);

        // We have reached a new bucket, check if there are too many variants in the bucket
        if (static_cast<long>(max_scores_in_bucket.size()) > max_variants_in_100bp_window)
        {
          assert(std::distance(window_it, map_it) == static_cast<long>(max_scores_in_bucket.size()));
          window_it = filter_window(window_it);
          map_it = window_it;
        }
        else
        {
          window_it = map_it;
        }

        // Update current bucket
        current_bucket = bucket;
        max_scores_in_bucket.clear();
      }

      // Add the max variant score
      long max_score{0};

      // Get the max score for this variant
      for (VariantSupport const & supp : map_it->second)
      {
        long const score = supp.get_score();

        if (score > max_score)
          max_score = score;
      }

      max_scores_in_bucket.push_back(max_score);
    }

    // check final bucket
    if (window_it != pool_varmap.end() && static_cast<long>(max_scores_in_bucket.size()) > max_variants_in_100bp_window)
    {
      filter_window(window_it);
    }
  }

  /** Break down variants to check if they are duplicates */
  // std::vector<VariantCandidate> broken_vars_to_add;
  // std::vector<std::vector<VariantSupport> > supports_to_add;

  for (auto map_it = pool_varmap.begin(); map_it != pool_varmap.end();) // no increment
  {
    VariantCandidate var(map_it->first);

    {
      long constexpr EXTRA_BASES_TO_ADD{5};

      for (long i{0}; i < EXTRA_BASES_TO_ADD; ++i)
        if (!var.add_base_in_front(false)) // false is add_N
          break;

      for (long i{0}; i < EXTRA_BASES_TO_ADD; ++i)
        if (!var.add_base_in_back(false)) // false is add_N
          break;
    }

    Variant var_cp;
    var_cp.abs_pos = var.abs_pos;
    var_cp.seqs = var.seqs;

    // Force that no variant is overlapping
    long const reach{-1};
    bool const is_no_variant_overlapping{true}; // Force it to be true so we wont use skyr
    bool const is_all_biallelic{false};
    Variant new_var(var_cp);
    new_var.infos["SBF1"] = ".";

    std::vector<Variant> new_broken_down_vars =
      break_down_variant(std::move(new_var), reach, is_no_variant_overlapping, is_all_biallelic);

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

    for (long i = 0; i < static_cast<long>(new_broken_down_vars.size()); ++i)
    {
      new_broken_down_var_candidates[i].abs_pos = new_broken_down_vars[i].abs_pos;
      new_broken_down_var_candidates[i].seqs = std::move(new_broken_down_vars[i].seqs);
      new_broken_down_var_candidates[i].flags = var.flags;
    }

    // Remove broken variants which are already in the varmap
    new_broken_down_var_candidates.erase(std::remove_if(new_broken_down_var_candidates.begin(),
                                                        new_broken_down_var_candidates.end(),
                                                        [&](VariantCandidate const & broken_var)
                                                        { return pool_varmap.find(broken_var) != pool_varmap.end(); }),
                                         new_broken_down_var_candidates.end());

    // Remove broken variants that are already in the graph
    new_broken_down_var_candidates.erase(std::remove_if(new_broken_down_var_candidates.begin(),
                                                        new_broken_down_var_candidates.end(),
                                                        [&](VariantCandidate const & broken_var)
                                                        { return graph.is_variant_in_graph(Variant(broken_var)); }),
                                         new_broken_down_var_candidates.end());

    if (new_broken_down_var_candidates.size() == 0)
    {
#ifndef NDEBUG
      print_log(log_severity::debug, __HERE__, " Erased ", map_it->first.print());
#endif                                    // NDEBUG
      map_it = pool_varmap.erase(map_it); // Erase the old variant
    }
    else
    {
      ++map_it; // Don't erase the old variant, just move on to the next one
    }
  }
  //*/
}

void VariantMap::clear()
{
  samples.clear();
  pool_varmap.clear();
  varmaps.clear();
}

#ifndef NDEBUG
void VariantMap::write_stats(std::string const & prefix)
{
  auto const & pn = samples[0];

  // Determine filename
  std::string const discovery_fn = Options::const_instance()->stats + "/" + prefix + "_" + pn + "_discovery.tsv.gz";

  // Open the stream
  std::stringstream discovery_ss;

  // Write header to the stream
  discovery_ss << "CHROM\tPOS\tREF\tALT\tSAMPLE\tGT\tAD\tDP\tHQ\tLQ\tPP\tMQ0\tCR\tVS\tG\tUP\tREV\t1\tSUP\tSUPR\tPL\n";

  // for (auto const & sample_name : samples)
  //  discovery_ss << "\t" << sample_name;

  // discovery_ss << '\n';

  for (auto map_it = pool_varmap.begin(); map_it != pool_varmap.end(); ++map_it)
  {
    Variant var(map_it->first);
    assert(var.seqs.size() == 2);
    auto contig_pos = absolute_pos.get_contig_position(var.abs_pos, gyper::graph.contigs);
    discovery_ss << contig_pos.first << "\t" << contig_pos.second << "\t";

    // REF
    std::copy(var.seqs[0].begin(), var.seqs[0].end(), std::ostream_iterator<char>(discovery_ss));
    discovery_ss << '\t';

    // ALT
    std::copy(var.seqs[1].begin(), var.seqs[1].end(), std::ostream_iterator<char>(discovery_ss));
    discovery_ss << '\t';

    // Add naive genotype calls
    std::sort(map_it->second.begin(),
              map_it->second.end(),
              [](VariantSupport const & a, VariantSupport const & b) -> bool { return a.pn_index < b.pn_index; });

    assert(Options::const_instance()->stats.size() > 0);
    std::size_t m = 0; // Index in map->second

    // Loop over all samples
    for (std::size_t p = 0; p < varmaps.size(); ++p)
    {
      // SAMPLE
      discovery_ss << samples[p] << '\t';

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
        assert(p == var_support.pn_index);

        coverage[1] = var_support.lq_support + var_support.hq_support;
        coverage[0] = var_support.depth - coverage[1];
        uint32_t phred0 = coverage[1] * 25u;
        uint32_t phred1 = var_support.depth * 3u;
        uint32_t phred2 = coverage[0] * 25u;
        uint32_t const normalizer = std::min(phred0, std::min(phred1, phred2));
        phred0 -= normalizer;
        phred1 -= normalizer;
        phred2 -= normalizer;
        assert(phred0 == 0 || phred1 == 0 || phred2 == 0);
        phred[0] = phred0 < 255 ? phred0 : 255;
        phred[1] = phred1 < 255 ? phred1 : 255;
        phred[2] = phred2 < 255 ? phred2 : 255;

        long const gt = std::distance(phred.begin(), std::find(phred.begin(), phred.end(), 0u));

        if (gt == 0)
          discovery_ss << "0/0\t";
        else if (gt == 1)
          discovery_ss << "0/1\t";
        else
          discovery_ss << "1/1\t";

        discovery_ss.setf(std::ios::fixed);
        discovery_ss.precision(2);
        discovery_ss << coverage[0] << ',' << coverage[1] << '\t' << var_support.depth << '\t' << var_support.hq_support
                     << '\t' << var_support.lq_support << '\t' << var_support.proper_pairs
                     << '\t'
                     //                     << var_support.mapq0 << '\t'
                     //                     << var_support.clipped << '\t'
                     << var_support.var_size << '\t' << var_support.growth << '\t'
                     << var_support.unique_positions.size() << '\t' << var_support.sequence_reversed << '\t'
                     << var_support.first_in_pairs << '\t' << var_support.get_corrected_support() << '\t'
                     << var_support.get_ratio() << '\t' << phred[0] << ',' << phred[1] << ',' << phred[2] << '\t';

        ++m;
      }
    }

    discovery_ss << '\n';
    // new_variant_vcf.variants.push_back(std::move(var));
  }

  write_gzipped_to_file(discovery_ss, discovery_fn, false /*append*/);
}

#endif // NDEBUG

void VariantMap::get_vcf(Vcf & new_variant_vcf, std::string const & output_name)
{
  new_variant_vcf.open(WRITE_BGZF_MODE, output_name);

#ifndef NDEBUG
  print_log(log_severity::debug, "[graphtyper::variant_map] Writing ", pool_varmap.size(), " variants.");
#endif // NDEBUG

  for (auto map_it = pool_varmap.begin(); map_it != pool_varmap.end(); ++map_it)
  {
    new_variant_vcf.variants.push_back(Variant(map_it->first));
  }
}

void save_variant_map(std::string const & path, VariantMap const & map)
{
  std::ofstream ofs(path.c_str(), std::ios::binary);

  if (!ofs.is_open())
  {
    print_log(log_severity::error, "[graphtyper::variant_map] Could not save variant_map at '", path, "'");

    std::exit(1);
  }

  cereal::BinaryOutputArchive oa(ofs);
  oa << map;
}

VariantMap load_variant_map(std::string const & path)
{
  VariantMap map;
  std::ifstream ifs(path.c_str(), std::ios::binary);

  if (!ifs.is_open())
  {
    print_log(log_severity::error, "[graphtyper::constructor] Could not load variant_map at '", path, "'");

    std::exit(1);
  }

  cereal::BinaryInputArchive ia(ifs);
  ia >> map;
  return map;
}

void VariantMap::set_samples(std::vector<std::string> const & new_samples)
{
  samples = new_samples;
  varmaps.resize(samples.size());
}

void VariantMap::load_many_variant_maps(std::string const & path)
{
  std::vector<std::string> paths;
  std::ifstream ifs(path.c_str());

  for (std::string line; std::getline(ifs, line);)
    paths.push_back(line);

  this->load_many_variant_maps(paths);
}

void VariantMap::load_many_variant_maps(std::vector<std::string> const & paths)
{
  this->pool_varmap.clear(); // Should be empty initially

  for (long p = 0; p < static_cast<long>(paths.size()); ++p)
  {
    auto const & file_path = paths[p];
    VariantMap new_map = load_variant_map(file_path);

    if (p == 0)
    {
      minimum_variant_support = new_map.minimum_variant_support;
      minimum_variant_support_ratio = new_map.minimum_variant_support_ratio;
    }

#ifndef NDEBUG
    // Sample index offset
    long sample_index_offset = samples.size();
#endif // NDEBUG

    // Copy samples
    std::copy(new_map.samples.begin(), new_map.samples.end(), std::back_inserter(samples));

    for (auto it = new_map.pool_varmap.begin(); it != new_map.pool_varmap.end(); ++it)
    {
#ifndef NDEBUG
      // Change all pn_indexes with offset
      for (auto variant_support : it->second)
        variant_support.pn_index += sample_index_offset;
#endif // NDEBUG

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

} // namespace gyper
