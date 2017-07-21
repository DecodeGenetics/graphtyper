#include <algorithm>
#include <assert.h>
#include <cstdlib>
#include <iterator>
#include <iostream>
#include <set>
#include <time.h>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/log/trivial.hpp>

#include <stations/join.hpp>
#include <stations/split.hpp>
#include <stations/station.hpp>

#include <graphtyper/graph/absolute_position.hpp>
#include <graphtyper/graph/graph_serialization.hpp>
#include <graphtyper/graph/haplotype_calls.hpp>
#include <graphtyper/graph/haplotype_extractor.hpp>
#include <graphtyper/typer/read_stats.hpp>
#include <graphtyper/typer/vcf_writer.hpp>
#include <graphtyper/utilities/graph_help_functions.hpp>
#include <graphtyper/utilities/options.hpp>
#include <graphtyper/utilities/vcf_help_functions.hpp>


namespace
{


bool
compare_gts(gyper::Genotype const & a, gyper::Genotype const & b)
{
  return a.id < b.id;
}


bool
are_genotype_paths_good(gyper::GenotypePaths const & geno)
{
  bool const fully_aligned = geno.all_paths_fully_aligned();

  if (geno.paths.size() == 0 || (!fully_aligned && (!geno.all_paths_unique() || geno.paths[0].size() < 95)))
    return false;

  if (gyper::Options::instance()->hq_reads)
  {
    if (!fully_aligned) // Require reads to be fully aligned
      return false;

    // Any path overlapping a variant must also not have too many mismatches
    for (auto const & path : geno.paths)
    {
      if (path.var_order.size() > 0 && path.mismatches > 2)
        return false;
    }
  }

  return true;
}


} // anon namespace


namespace gyper
{

VcfWriter::VcfWriter(std::vector<std::string> const & samples, uint32_t variant_distance)
  : pn(samples[0])
{
  BOOST_LOG_TRIVIAL(info) << "[graphtyper::vcf_writer] Getting all haplotypes.";
  BOOST_LOG_TRIVIAL(info) << "[graphtyper::vcf_writer] Number of variant nodes in graph "
                          << graph.var_nodes.size();

  haplotypes = gyper::graph.get_all_haplotypes(variant_distance);
  BOOST_LOG_TRIVIAL(info) << "[graphtyper::vcf_writer] Got "
                          << haplotypes.size()
                          << " haplotypes.";

  // Insert items in the id to haplotype map
  for (unsigned i = 0; i < haplotypes.size(); ++i)
  {
    haplotypes[i].clear_and_resize_samples(samples.size());
    std::vector<uint32_t> gt_ids = haplotypes[i].get_genotype_ids();

    for (unsigned j = 0; j < gt_ids.size(); ++j)
      id2hap[gt_ids[j]] = {i, j};
  }
}


void
VcfWriter::update_haplotype_scores_from_paths(std::vector<GenotypePaths> & genos,
                                              std::size_t const pn_index
  )
{
  std::lock_guard<std::mutex> lock(haplotype_mutex);
  assert(!Options::instance()->is_segment_calling);

  for (auto & geno : genos)
  {
    if (are_genotype_paths_good(geno))
      update_haplotype_scores_from_path(geno, pn_index, 0 /*unpaired*/);
  }
}


void
VcfWriter::update_haplotype_scores_from_paths(
  std::vector<std::pair<GenotypePaths, GenotypePaths> > & genos,
  std::size_t const pn_index
  )
{
  std::lock_guard<std::mutex> lock(haplotype_mutex);

  for (auto & geno : genos)
  {
    // if (Options::instance()->is_segment_calling)
    // {
    //   if (geno.first.paths.size() > 0 && geno.first.all_paths_fully_aligned() &&
    //       geno.second.paths.size() > 0 && geno.second.all_paths_fully_aligned()
    //       )
    //   {
    //     update_haplotype_scores_from_path(geno.first, pn_index);
    //     update_haplotype_scores_from_path(geno.second, pn_index);
    //   }
    // }
    // else
    //{
    bool const READ1_IS_GOOD = are_genotype_paths_good(geno.first);
    bool const READ2_IS_GOOD = are_genotype_paths_good(geno.second);

    if (Options::instance()->hq_reads)
    {
      // Make sure both read pairs support the same haplotypes
      if (READ1_IS_GOOD && READ2_IS_GOOD)
      {
        update_haplotype_scores_from_path(geno.first, pn_index, 1 /*first in pair*/);
        update_haplotype_scores_from_path(geno.second, pn_index, 2 /*second in pair*/);
      }
    }
    else
    {
      if (READ1_IS_GOOD && READ2_IS_GOOD)
      {
        update_haplotype_scores_from_path(geno.first, pn_index, 1 /*first in pair*/);
        update_haplotype_scores_from_path(geno.second, pn_index, 2 /*second in pair*/);

        if (Options::instance()->stats.size() > 0)
        {
          for (auto & hap : haplotypes)
            hap.hap_samples[pn_index].stats->pair_stats.clear();
        }
      }
      else if (READ1_IS_GOOD)
      {
        update_haplotype_scores_from_path(geno.first, pn_index, 0 /*unpaired*/);
      }
      else if (READ2_IS_GOOD)
      {
        update_haplotype_scores_from_path(geno.second, pn_index, 0 /*unpaired*/);
      }
    }
  }
}


void
VcfWriter::update_statistics(GenotypePaths & geno, std::size_t const pn_index, unsigned const read_pair)
{
  // If the alignment is unique, log the sequence successor and connected haplotypes
  if (geno.paths.size() == 1)
  {
    auto const & path = geno.paths[0];

    for (unsigned i = 0; i < path.var_order.size(); ++i)
    {
      if (path.nums[i].count() != 1)
        continue;

      uint32_t b = 0;

      while (!path.nums[i].test(b))
      {
        ++b;
        assert (b < path.nums[i].size());
      }

      uint32_t v = 0;

      while (graph.var_nodes[v].get_label().order != path.var_order[i] || graph.var_nodes[v].get_label().variant_num != b)
      {
        ++v;
        assert (v < graph.var_nodes.size());
      }

      uint32_t const reach = graph.var_nodes[v].get_label().reach();
      assert(id2hap.count(path.var_order[i]) == 1);
      std::pair<uint32_t, uint32_t> const type_ids = id2hap.at(path.var_order[i]); //hap_id = first, gen_id = second
      assert(type_ids.first < haplotypes.size());
      assert(pn_index < haplotypes[type_ids.first].hap_samples.size());
      assert(haplotypes[type_ids.first].hap_samples[pn_index].stats);
      auto const & stats = haplotypes[type_ids.first].hap_samples[pn_index].stats;

      // Read pair information
      if (read_pair == 1)
      {
        stats->pair_stats.push_back({type_ids.first, b}); // Add this pair of haplotype ID and genotype ID
      }
      else if (read_pair == 2)
      {
        // Check all other haplotypes
        for (auto & hap : haplotypes)
        {
          auto const & stats_first = hap.hap_samples[pn_index].stats;

          // Add read pair information to stats
          for (auto const & first_read : stats_first->pair_stats)
          {
            assert (first_read.gt_id < haplotypes[first_read.hap_id].hap_samples[pn_index].stats->pair_info.size());
            assert (b < haplotypes[type_ids.first].hap_samples[pn_index].stats->pair_info.size());
            haplotypes[first_read.hap_id].hap_samples[pn_index].stats->pair_info[first_read.gt_id].push_back({type_ids.first, b});
            haplotypes[type_ids.first].hap_samples[pn_index].stats->pair_info[b].push_back({first_read.hap_id, first_read.gt_id});
          }
        }
      }

      uint32_t constexpr minW = 1;
      uint32_t constexpr W = 1;

      if (path.end_correct_pos() >= reach + minW && path.end_correct_pos() - reach < geno.read.size())
      {
        assert (b < stats->successor.size());
        auto start_it = geno.read.end() - path.end_correct_pos() + reach;
        assert (std::distance(start_it, geno.read.end()) >= minW);

        if (std::distance(start_it, geno.read.end()) > W)
          stats->successor[b].push_back(std::vector<char>(start_it, start_it + W));
        else
          stats->successor[b].push_back(std::vector<char>(start_it, geno.read.end()));

        auto qual_start_it = geno.qual.end() - path.end_correct_pos() + reach;

        // Add it again if all bases are high quality
        if (std::count(qual_start_it, qual_start_it + minW, 33) == 0)
        {
          if (std::distance(start_it, geno.read.end()) > W)
            stats->successor[b].push_back(std::vector<char>(start_it, start_it + W));
          else
            stats->successor[b].push_back(std::vector<char>(start_it, geno.read.end()));
        }
      }

      // Code for predecessor
      // if (graph.var_nodes[v].get_label().order >= path.start_correct_pos() + minW &&
      //     graph.var_nodes[v].get_label().order - path.start_correct_pos() < geno.read.size()
      //    )
      // {
      //   assert (b < stats->predecessor.size());
      //   auto end_it = geno.read.begin() + graph.var_nodes[v].get_label().order - path.start_correct_pos();
      //   auto qual_end_it = geno.qual.begin() + graph.var_nodes[v].get_label().order - path.start_correct_pos();
      //
      //   if (std::count(qual_end_it - minW, qual_end_it, 33) == 0)
      //   {
      //     if (std::distance(geno.read.begin(), end_it) > W)
      //       stats->predecessor[b].push_back(std::vector<char>(end_it - W, end_it));
      //     else
      //       stats->predecessor[b].push_back(std::vector<char>(geno.read.begin(), end_it));
      //   }
      // }
    }
  }
}


void
VcfWriter::update_haplotype_scores_from_path(GenotypePaths & geno, std::size_t const pn_index, unsigned const read_pair)
{
  assert (are_genotype_paths_good(geno));

  // Quality metrics
  bool const fully_aligned = geno.all_paths_fully_aligned();
  bool const non_unique_paths = !geno.all_paths_unique();
  std::size_t const mismatches = geno.paths[0].mismatches;
  assert (geno.read.size() >= 63);

  if (Options::instance()->stats.size() > 0) // Update statistics if '--stats' was passed
    this->update_statistics(geno, pn_index, read_pair);

  std::vector<uint32_t> recent_ids;
  bool has_low_quality_snp = false;

  for (auto p_it = geno.paths.begin(); p_it != geno.paths.end(); ++p_it)
  {
    assert(p_it->var_order.size() == p_it->nums.size());

    for (std::size_t i = 0; i < p_it->var_order.size(); ++i)
    {
      assert(id2hap.count(p_it->var_order[i]) == 1);
      std::pair<uint32_t, uint32_t> const type_ids = id2hap.at(p_it->var_order[i]); // hap_id = first, gen_id = second
      assert(type_ids.first < haplotypes.size());

      // Check if the genotype is a SNP
      if (!has_low_quality_snp && graph.is_snp(type_ids.second))
      {
        uint32_t const offset = p_it->var_order[i] - p_it->start_correct_pos();

        if (offset < geno.qual.size())
        {
          assert(geno.qual[offset] >= 33);
          has_low_quality_snp = geno.qual[offset] - 33 <= 25;
        }
      }

      haplotypes[type_ids.first].add_explanation(type_ids.second, p_it->nums[i]);
      recent_ids.push_back(type_ids.first);

      // Add coverage if the explanation is unique
      if (p_it->nums[i].count() == 1)
      {
        // Check which bit is set
        uint16_t b = 0;

        while (not p_it->nums[i].test(b))
          ++b;

        haplotypes[type_ids.first].add_coverage(type_ids.second, b);
      }
      else /* Otherwise the coverage is ambigous => 0xFFFEu */
      {
        haplotypes[type_ids.first].add_coverage(type_ids.second, 0xFFFEu);
      }
    }
  }

  std::sort(recent_ids.begin(), recent_ids.end());
  auto last = std::unique(recent_ids.begin(), recent_ids.end());

  // After each read, move the "explain" to the score vector.
  for (auto it = recent_ids.begin(); it != last; ++it)
  {
    assert(*it < haplotypes.size());
    assert(pn_index < haplotypes[*it].hap_samples.size());
    assert (geno.paths.size() > 0);

    // Move INFO to variant statistics. This needs to called before 'coverage_to_gts', because it uses the coverage of each variant.
    haplotypes[*it].clipped_reads_to_stats(geno.all_paths_fully_aligned());
    haplotypes[*it].mapq_to_stats(geno.mapq);
    haplotypes[*it].strand_to_stats(geno.forward_strand, geno.is_first_in_pair);
    haplotypes[*it].realignment_to_stats(geno.is_originally_unaligned,
                                         geno.is_originally_clipped,
                                         geno.original_pos /*original_pos*/,
                                         absolute_pos.get_contig_position(geno.paths[0].start_correct_pos()).second /*new_pos*/
                                         );

    haplotypes[*it].graph_complexity_to_stats();

    // Update the likelihood scores
    haplotypes[*it].explain_to_score(pn_index, has_low_quality_snp, non_unique_paths, geno.mapq, fully_aligned, mismatches);
    haplotypes[*it].coverage_to_gts(pn_index); // Update the coverage
  }
}


THapCalls
VcfWriter::get_haplotype_calls() const
{
  THapCalls hap_calls;
  hap_calls.reserve(haplotypes.size());

  for (auto const & hap : haplotypes)
    hap_calls.push_back(std::make_pair(hap.get_haplotype_calls(), hap.gts));

  return hap_calls;
}


std::vector<Genotype>
VcfWriter::get_gts() const
{
  // Get all genotypes
  std::vector<Genotype> gts;

  for (auto hap_it = haplotypes.cbegin(); hap_it != haplotypes.cend(); ++hap_it)
  {
    std::vector<Genotype> new_gts(hap_it->gts);

    if (new_gts.size() == 0)
      BOOST_LOG_TRIVIAL(error) << "A haplotype contained no genotypes.";

    std::move(new_gts.begin(), new_gts.end(), std::back_inserter(gts));
  }

  BOOST_LOG_TRIVIAL(info) << "[graphtyper::vcf_writer] Number of markers are " << gts.size();
  std::sort(gts.begin(), gts.end(), compare_gts);
  return gts;
}


std::vector<uint32_t>
VcfWriter::explain_map_to_haplotype_scores(std::size_t const pn_index, ExplainMap const & explain_map)
{
  std::shared_ptr<ExplainMap> shared_explain_map = std::make_shared<ExplainMap>(explain_map);
  std::shared_ptr<std::size_t> shared_pn_index = std::make_shared<std::size_t>(pn_index);
  std::size_t const num = explain_map.cbegin()->second.size();
  std::vector<uint32_t> hap_scores((num + 1) * num / 2, 0u);
  auto hap_scores_split = stations::split(hap_scores, Options::instance()->threads);

  auto update_scores = [this](std::shared_ptr<ExplainMap> explain_map,
                              std::shared_ptr<std::vector<uint32_t> > hap_scores,
                              std::shared_ptr<uint32_t> i_start,
                              std::shared_ptr<std::size_t> pn_index
                              )
  {
    for (auto it = explain_map->cbegin(); it != explain_map->cend(); ++it)
    {
      std::pair<uint16_t, uint16_t> xy = to_pair(*i_start);
      uint16_t x = xy.first;
      uint16_t y = xy.second;

      for (uint32_t i = 0; i < hap_scores->size(); ++i)
      {
        // std::cout << "x, y = " << x << " " << y << std::endl;
        assert (x < it->second.size());
        assert (y < it->second.size());
        // assert (i < hap_scores.size());

        (*hap_scores)[i] += this->haplotypes[it->first].best_score_of_a_path(*pn_index, it->second[x], it->second[y]);

        if (x == y)
        {
          x = 0;
          ++y;
        }
        else
        {
          ++x;
        }
      }
    }
  };

  {
    stations::Station hap_scores_station(Options::instance()->threads);
    uint32_t i = 0;

    for (unsigned j = 0; j < hap_scores_split.size(); ++j)
    {
      hap_scores_station.add_to_thread(j /*thread id*/,
                                       update_scores /*work*/,
                                       shared_explain_map,
                                       hap_scores_split[j],
                                       std::make_shared<uint32_t>(i),
                                       shared_pn_index
                                     );
      i += hap_scores_split[j]->size();
    }
  }

  stations::join(hap_scores, hap_scores_split);
  return hap_scores;
}


uint32_t
VcfWriter::explain_map_specific_indexes_to_haplotype_scores(
  std::size_t const pn_index,
  std::pair<uint32_t, uint32_t> const index,
  ExplainMap const & explain_map
) const
{
  uint32_t hap_score = 0;

  for (auto it = explain_map.cbegin(); it != explain_map.cend(); ++it)
  {
    // uint32_t const cnum = haplotypes[it->first].get_genotype_num();
    // std::vector<uint16_t> max_scores_per_y(cnum, 0u);
    //
    // for (uint32_t y = 0; y < cnum; ++y)
    // {
    //   auto score_begin = haplotypes[it->first].hap_samples[pn_index].log_score.begin() + to_index(0, y);
    //   auto score_end = haplotypes[it->first].hap_samples[pn_index].log_score.begin() + to_index(0, y + 1);
    //   auto max_it = std::max_element(score_begin, score_end);
    //   assert(max_it != score_end);
    //   // std::cout << "Max iterator is " << *max_it << std::endl;
    //   max_scores_per_y[y] = *max_it;
    // }

    hap_score += haplotypes[it->first].best_score_of_a_path(pn_index, it->second[index.first], it->second[index.second]);
  }


  return hap_score;
}


void
VcfWriter::find_path_explanation(GenotypePaths const & gt_path,
                                 std::vector<std::pair<uint32_t, std::bitset<MAX_NUMBER_OF_HAPLOTYPES> > > & ids_and_path_explain
                                 )
{
  // ids_and_path_explain is a vector of haplotype_id and haplotype_explain pairs
  assert (ids_and_path_explain.size() == 0);
  std::vector<uint32_t> recent_ids;

  for (auto p_it = gt_path.paths.cbegin(); p_it != gt_path.paths.cend(); ++p_it)
  {
    assert (p_it->var_order.size() == p_it->nums.size());

    for (unsigned i = 0; i < p_it->var_order.size(); ++i)
    {
      std::pair<uint32_t, uint32_t> type_ids = id2hap.at(p_it->var_order[i]); //hap_id = first, gen_id = second
      assert (type_ids.first < haplotypes.size());
      haplotypes[type_ids.first].add_explanation(type_ids.second, p_it->nums[i]);
      recent_ids.push_back(type_ids.first);
    }
  }

  std::sort(recent_ids.begin(), recent_ids.end());
  auto last_it = std::unique(recent_ids.begin(), recent_ids.end());

  for (auto it = recent_ids.begin(); it != last_it; ++it)
  {
    assert (*it < haplotypes.size());
    std::bitset<MAX_NUMBER_OF_HAPLOTYPES> path_explain = haplotypes[*it].explain_to_path_explain();
    ids_and_path_explain.push_back(std::make_pair(*it, path_explain));
    assert (ids_and_path_explain.back().second.any());
  }
}


void
VcfWriter::generate_statistics(std::vector<std::string> const & pns)
{
  if (Options::instance()->stats.size() == 0 || pns.size() == 0)
    return;

  std::string const hap_stats = Options::instance()->stats + "/" + pns[0] + "_haplotypes.tsv";
  BOOST_LOG_TRIVIAL(info) << "[graphtyper::vcf_writer] Generating haplotype statistics to " << hap_stats;
  std::ofstream hap_file(hap_stats.c_str());

  // Loop samples
  for (unsigned i = 0; i < pns.size(); ++i)
  {
    auto const & pn = pns[i];

    // Haplotype statistics
    for (unsigned ps = 0; ps < haplotypes.size(); ++ps)
    {
      auto const & hap = haplotypes[ps];
      auto const & gts = hap.gts;
      assert (hap.hap_samples.size() == pns.size());
      assert (i < hap.hap_samples.size());
      auto const & hap_sample = hap.hap_samples[i];
      assert (hap_sample.stats);
      uint32_t const cnum = hap.get_genotype_num();

      for (uint32_t c = 0; c < cnum; ++c)
      {
        std::vector<uint32_t> gt_calls;
        uint32_t call = c;
        uint32_t q = cnum;

        for (uint32_t i = 0; i < gts.size(); ++i)
        {
          q /= gts[i].num;
          assert (q != 0);
          gt_calls.push_back(call / q);
          call %= q;
        }

        hap_file << pn << '\t' << ps << '\t' << c; // << "\t" << gt_calls[0];

        //for (uint32_t l = 1; l < gt_calls.size(); ++l)
        //  hap_file << ">" << gt_calls[l];

        hap_file << '\t' << static_cast<uint16_t>(hap_sample.stats->hap_coverage[c])
                 << '\t' << static_cast<uint16_t>(hap_sample.stats->hap_unique_coverage[c])
                 << '\t';

        if (hap_sample.stats->successor[c].size() > 0)
        {
          hap_file << std::string(hap_sample.stats->successor[c][0].begin(), hap_sample.stats->successor[c][0].end());

          for (auto it = hap_sample.stats->successor[c].begin() + 1; it != hap_sample.stats->successor[c].end(); ++it)
            hap_file << "," << std::string(it->begin(), it->end());
        }
        else
        {
          hap_file << "N/A";
        }

        hap_file << '\t';

        if (hap_sample.stats->pair_info[c].size() > 0)
        {
          hap_file << "(" << hap_sample.stats->pair_info[c][0].first
                   << "," << hap_sample.stats->pair_info[c][0].second
                   << ")";

          for (auto it = hap_sample.stats->pair_info[c].begin() + 1; it != hap_sample.stats->pair_info[c].end(); ++it)
            hap_file << ",(" << it->first << "," << it->second << ")";
        }
        else
        {
          hap_file << "N/A";
        }

        // if (hap_sample.stats->predecessor[c].size() > 0)
        // {
        //   hap_file << std::string(hap_sample.stats->predecessor[c][0].begin(), hap_sample.stats->predecessor[c][0].end());
        //
        //   for (auto it = hap_sample.stats->predecessor[c].begin() + 1; it != hap_sample.stats->predecessor[c].end(); ++it)
        //     hap_file << "," << std::string(it->begin(), it->end());
        // }
        // else
        // {
        //   hap_file << "N/A";
        // }

        // hap_file << "\tPredecessors";
        //
        // for (auto const & p : hap_sample.stats->predecessor[c])
        //   hap_file << ":" << std::string(p.begin(), p.end());

        hap_file << "\n";
      }
    }
  }
}


} // namespace gyper
