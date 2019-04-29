#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iterator>
#include <iostream>
#include <set>
#include <time.h>
#include <utility>

#include <boost/log/trivial.hpp>

#include <paw/station.hpp>

#include <graphtyper/graph/absolute_position.hpp>
#include <graphtyper/graph/graph_serialization.hpp>
#include <graphtyper/graph/haplotype_calls.hpp>
#include <graphtyper/graph/haplotype_extractor.hpp>
#include <graphtyper/typer/read_stats.hpp>
#include <graphtyper/typer/vcf_writer.hpp>
#include <graphtyper/utilities/graph_help_functions.hpp>
#include <graphtyper/utilities/io.hpp>
#include <graphtyper/utilities/options.hpp>
#include <graphtyper/utilities/vcf_help_functions.hpp>


namespace
{


bool
are_genotype_paths_perfect(gyper::GenotypePaths const & geno)
{
  if (geno.paths.size() == 0 || !geno.all_paths_fully_aligned() || geno.paths[0].mismatches > 0)
    return false;

  return true;
}

bool
are_genotype_paths_good(gyper::GenotypePaths const & geno)
{
  if (geno.paths.size() == 0)
    return false;

  bool const fully_aligned = geno.all_paths_fully_aligned();

  if (!fully_aligned && (!geno.all_paths_unique() || geno.paths[0].size() < 63))
    return false;

  double const mismatch_ratio =
    static_cast<double>(geno.paths[0].mismatches) / static_cast<double>(geno.paths[0].size());

  if (mismatch_ratio > 0.05)
    return false;

  if (!fully_aligned && mismatch_ratio > 0.025)
    return false;

  if (gyper::graph.is_sv_graph)
  {
    if (!fully_aligned || geno.paths[0].size() < 90 || mismatch_ratio > 0.03)
      return false;
  }

  if (gyper::Options::instance()->hq_reads)
  {
    if (!fully_aligned || geno.paths[0].size() < 90 || mismatch_ratio > 0.025)
      return false;
  }

  return true;
}


} // anon namespace


namespace gyper
{

VcfWriter::VcfWriter(std::vector<std::string> const & samples, uint32_t variant_distance)
  : pns(samples)
  , pn(samples[0])
  , haplotypes(gyper::graph.get_all_haplotypes(variant_distance))
{
  BOOST_LOG_TRIVIAL(info) << "[graphtyper::vcf_writer] Number of variant nodes in graph "
                          << graph.var_nodes.size();
  BOOST_LOG_TRIVIAL(info) << "[graphtyper::vcf_writer] Got "
                          << haplotypes.size()
                          << " haplotypes.";

  // Insert items in the id to haplotype map
  for (unsigned i = 0; i < haplotypes.size(); ++i)
  {
    haplotypes[i].clear_and_resize_samples(samples.size());
    std::vector<uint32_t> gt_ids = haplotypes[i].get_genotype_ids();

    for (unsigned j = 0; j < gt_ids.size(); ++j)
    {
      id2hap[gt_ids[j]] =
        std::make_pair<uint32_t, uint32_t>(static_cast<uint32_t>(i), static_cast<uint32_t>(j));
    }
  }
}


void
VcfWriter::update_haplotype_scores_from_paths(std::vector<GenotypePaths> & genos,
                                              std::size_t const pn_index
  )
{
  if (Options::instance()->is_perfect_alignments_only)
  {
    // Perfect alignments
    std::lock_guard<std::mutex> lock(haplotype_mutex);

    for (auto & geno : genos)
    {
      if (are_genotype_paths_perfect(geno))
        update_haplotype_scores_from_path(geno, pn_index);
    }
  }
  else
  {
    // Good alignments (default)
    std::lock_guard<std::mutex> lock(haplotype_mutex);

    for (auto & geno : genos)
    {
      if (are_genotype_paths_good(geno))
        update_haplotype_scores_from_path(geno, pn_index);
    }
  }

  if (Options::instance()->stats.size() > 0)
    update_statistics(genos, pn_index);
}


void
VcfWriter::update_haplotype_scores_from_paths(
  std::vector<std::pair<GenotypePaths, GenotypePaths> > & genos,
  std::size_t const pn_index
  )
{
  if (Options::instance()->is_perfect_alignments_only)
  {
    std::lock_guard<std::mutex> lock(haplotype_mutex);

    for (auto & geno : genos)
    {
      bool const READ1_IS_GOOD = are_genotype_paths_good(geno.first);
      bool const READ2_IS_GOOD = are_genotype_paths_good(geno.second);

      // Require both reads to be at least good
      if (READ1_IS_GOOD && READ2_IS_GOOD)
      {
        if (are_genotype_paths_perfect(geno.first))
          update_haplotype_scores_from_path(geno.first, pn_index);

        if (are_genotype_paths_perfect(geno.second))
          update_haplotype_scores_from_path(geno.second, pn_index);
      }
    }
  }
  else
  {
    std::lock_guard<std::mutex> lock(haplotype_mutex);

    for (auto & geno : genos)
    {
      bool const READ1_IS_GOOD = are_genotype_paths_good(geno.first);
      bool const READ2_IS_GOOD = are_genotype_paths_good(geno.second);

      if (Options::instance()->hq_reads)
      {
        // Make sure both read pairs support the same haplotypes
        if (READ1_IS_GOOD && READ2_IS_GOOD)
        {
          update_haplotype_scores_from_path(geno.first, pn_index);
          update_haplotype_scores_from_path(geno.second, pn_index);
        }
      }
      else
      {
        if (READ1_IS_GOOD && READ2_IS_GOOD)
        {
          update_haplotype_scores_from_path(geno.first, pn_index);
          update_haplotype_scores_from_path(geno.second, pn_index);
        }
        else if (READ1_IS_GOOD)
        {
          update_haplotype_scores_from_path(geno.first, pn_index);
        }
        else if (READ2_IS_GOOD)
        {
          update_haplotype_scores_from_path(geno.second, pn_index);
        }
      }
    }
  }

  if (Options::instance()->stats.size() > 0)
    update_statistics(genos, pn_index);
}


void
VcfWriter::print_variant_group_details() const
{
  std::string const hap_stats_fn = Options::instance()->stats + "/" + pns[0] + "_variant_group_details.tsv.gz";
  BOOST_LOG_TRIVIAL(info) << "[graphtyper::vcf_writer] Generating variant group info statistics to " << hap_stats_fn;
  std::stringstream hap_file;

  // Write header file
  hap_file << "groupID\talleleID\tcontig\tposition\tsequence\n";

  for (unsigned ps = 0; ps < haplotypes.size(); ++ps)
  {
    auto const & hap = haplotypes[ps];
    uint32_t const cnum = hap.get_genotype_num();

    for (uint32_t c = 0; c < cnum; ++c)
    {
      assert (hap.gts.size() > 0);
      uint32_t const abs_pos = hap.gts[0].id;
      std::vector<char> seq = graph.get_sequence_of_a_haplotype_call(hap.gts, c);
      assert (seq.size() > 1);
      auto contig_pos = absolute_pos.get_contig_position(abs_pos);

      hap_file << ps << "\t" << c << "\t"
               << contig_pos.first << "\t" << contig_pos.second << "\t"
               << std::string(seq.begin() + 1, seq.end())
               << "\n";
     }
  }

  std::lock_guard<std::mutex> lock(io_mutex);
  write_gzipped_to_file(hap_file, hap_stats_fn);
}


void
VcfWriter::print_variant_details() const
{
  std::string const variant_details_fn =
    Options::instance()->stats + "/" + pns[0] + "_variant_details.tsv.gz";

  BOOST_LOG_TRIVIAL(info) << "[graphtyper::vcf_writer] Generating variant info statistics to "
                          << variant_details_fn;

  std::stringstream variant_file;

  // Write header file
  variant_file << "variantID\tcontig\tposition\tallele_num\tsequence\n";

  for (std::size_t v = 0; v < graph.var_nodes.size(); ++v)
  {
    auto const & label = graph.var_nodes[v].get_label();
    auto contig_pos = absolute_pos.get_contig_position(label.order);

    variant_file << v << "\t"
                 << contig_pos.first << "\t"
                 << contig_pos.second << "\t"
                 << label.variant_num << "\t"
                 << std::string(label.dna.begin(), label.dna.end())
                 << "\n";
  }

  std::lock_guard<std::mutex> lock(io_mutex);
  write_gzipped_to_file(variant_file, variant_details_fn);
}


void
VcfWriter::print_statistics_headers() const
{
  assert (Options::instance()->stats.size() > 0);

  // Get filenames
  std::string const read_details_fn =
    Options::instance()->stats + "/" + pn + "_read_details.tsv.gz";

  std::string const path_details_fn =
    Options::instance()->stats + "/" + pn + "_read_path_details.tsv.gz";

  std::stringstream read_ss;
  std::stringstream path_ss;
  read_ss << "query\t" << "read_group\t" << "sample\t" << "read\t" << "read_qual\t"
          << "alignment_length\t" << "mapping_quality\t" << "original_mapped_pos\t"
          << "ml_insert_size\t" << "is_originally_unaligned\t"  << "is_originally_clipped" << "\n";

  path_ss << "query\t" << "pathID\t" << "read_start_index\t" << "read_end_index\t"
          << "num_mismatches\t" << "strand\t" << "contig\t" << "alignment_begin\t"
          << "alignment_end\t" << "overlapping_variant_nodes\n";

  std::lock_guard<std::mutex> lock(io_mutex);
  write_gzipped_to_file(read_ss, read_details_fn, false /*append*/);
  write_gzipped_to_file(path_ss, path_details_fn, false /*append*/);
}


void
VcfWriter::print_geno_statistics(std::stringstream & read_ss,
                                 std::stringstream & path_ss,
                                 GenotypePaths const & geno,
                                 std::size_t const pn_index
  )
{
  std::stringstream id;
  assert (geno.details);
  id << pns[pn_index] << "_" << geno.details->query_name << "/" << (geno.is_first_in_pair ? 1 : 2);

  read_ss << id.str() << "\t"
          << geno.details->read_group << "\t"
          << pns[pn_index] << "\t"
          << std::string(geno.read.begin(), geno.read.end()) << "\t"
          << std::string(geno.qual.begin(), geno.qual.end()) << "\t"
          << geno.longest_path_length << "\t"
          << static_cast<uint16_t>(geno.mapq) << "\t"
          << geno.original_pos << "\t";

  if (geno.ml_insert_size == 0x7FFFFFFFl)
    read_ss << "NA";
  else
    read_ss << geno.ml_insert_size;

  read_ss << "\t"
          << (geno.is_originally_unaligned ? "Y" : "N") << "\t"
          << (geno.is_originally_clipped ? "Y" : "N") << "\n";

  for (std::size_t p = 0; p < geno.paths.size(); ++p)
  {
    auto const & path = geno.paths[p];
    uint32_t const ref_reach_start = path.start_ref_reach_pos();
    uint32_t const ref_reach_end = path.end_ref_reach_pos();

    auto const contig_pos_start = absolute_pos.get_contig_position(ref_reach_start);
    auto const contig_pos_end = absolute_pos.get_contig_position(ref_reach_end);

    std::vector<std::size_t> overlapping_vars;

    // Find
    for (std::size_t i = 0; i < path.var_order.size(); ++i)
    {
      auto const & var_order = path.var_order[i];
      auto const & num = path.nums[i];

      auto find_it = id2hap.find(var_order);
      assert (find_it->second.first < haplotypes.size());
      auto const & hap = haplotypes[find_it->second.first];
      assert (find_it->second.second < hap.gts.size());
      auto const & gt = hap.gts[find_it->second.second];

      for (uint32_t g = 0; g < gt.num; ++g)
      {
        if (num.test(g))
          overlapping_vars.push_back(gt.first_variant_node + g);
      }
    }

    path_ss << id.str() << "\t"
            << p << "\t"
            << path.read_start_index << "\t"
            << path.read_end_index << "\t"
            << path.mismatches << "\t"
            << (geno.forward_strand ? "F" : "B") << "\t"
            << contig_pos_start.first << "\t"
            << contig_pos_start.second << "\t"
            << contig_pos_end.second << "\t";

    std::sort(overlapping_vars.begin(), overlapping_vars.end());

    if (overlapping_vars.size() == 0)
    {
      path_ss << "NA";
    }
    else
    {
      path_ss << overlapping_vars[0];

      for (std::size_t o = 1; o < overlapping_vars.size(); ++o)
        path_ss << "," << overlapping_vars[o];
    }

    path_ss << "\n";
  }
}


void
VcfWriter::update_statistics(std::vector<GenotypePaths> & genos, std::size_t const pn_index)
{
  // Get filenames
  std::string const read_details_fn =
    Options::instance()->stats + "/" + pn + "_read_details.tsv.gz";

  std::string const path_details_fn =
    Options::instance()->stats + "/" + pn + "_read_path_details.tsv.gz";

  // Write reads and paths to these streams and then write those in a gzipped file
  std::stringstream read_ss;
  std::stringstream path_ss;

  for (auto const & geno : genos)
    print_geno_statistics(read_ss, path_ss, geno, pn_index);

  std::lock_guard<std::mutex> lock(io_mutex);
  write_gzipped_to_file(read_ss, read_details_fn, true /*append*/);
  write_gzipped_to_file(path_ss, path_details_fn, true /*append*/);
}


void
VcfWriter::update_statistics(std::vector<std::pair<GenotypePaths, GenotypePaths> > & genos,
                             std::size_t pn_index
  )
{
  // Get filenames
  std::string const read_details_fn =
    Options::instance()->stats + "/" + pn + "_read_details.tsv.gz";

  std::string const path_details_fn =
    Options::instance()->stats + "/" + pn + "_read_path_details.tsv.gz";

  // Write reads and paths to these streams and then write those in a gzipped file
  std::stringstream read_ss;
  std::stringstream path_ss;

  for (auto const & geno : genos)
  {
    print_geno_statistics(read_ss, path_ss, geno.first, pn_index);
    print_geno_statistics(read_ss, path_ss, geno.second, pn_index);
  }

  std::lock_guard<std::mutex> lock(io_mutex);
  write_gzipped_to_file(read_ss, read_details_fn, true /*append*/);
  write_gzipped_to_file(path_ss, path_details_fn, true /*append*/);
}


void
VcfWriter::update_haplotype_scores_from_path(GenotypePaths & geno,
                                             std::size_t const pn_index
  )
{
  assert(are_genotype_paths_good(geno));

  // Quality metrics
  bool const fully_aligned = geno.all_paths_fully_aligned();
  bool const non_unique_paths = !geno.all_paths_unique();
  std::size_t const mismatches = geno.paths[0].mismatches;
  std::vector<uint32_t> recent_ids;

  for (auto p_it = geno.paths.begin(); p_it != geno.paths.end(); ++p_it)
  {
    assert(p_it->var_order.size() == p_it->nums.size());

    for (std::size_t i = 0; i < p_it->var_order.size(); ++i)
    {
      assert(id2hap.count(p_it->var_order[i]) == 1);
      std::pair<uint32_t, uint32_t> const type_ids = id2hap.at(p_it->var_order[i]); // hap_id = first, gen_id = second

      assert(type_ids.first < haplotypes.size());
      assert(p_it->nums[i].any());

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
      else /* Otherwise set the coverage is ambiguous */
      {
        haplotypes[type_ids.first].add_coverage(type_ids.second, 1);

        if (p_it->nums[i].test(0))
          haplotypes[type_ids.first].add_coverage(type_ids.second, 0);
        else
          haplotypes[type_ids.first].add_coverage(type_ids.second, 2);
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
    assert(geno.paths.size() > 0);

    auto & haplotype = haplotypes[*it];

    // Move INFO to variant statistics. This needs to called before 'coverage_to_gts', because it uses the coverage of each variant.
    haplotype.clipped_reads_to_stats(geno.all_paths_fully_aligned());
    haplotype.mapq_to_stats(geno.mapq);
    haplotype.strand_to_stats(geno.forward_strand, geno.is_first_in_pair);
    haplotype.realignment_to_stats(geno.is_originally_unaligned,
                                   geno.is_originally_clipped,
                                   geno.original_pos /*original_pos*/,
                                   absolute_pos.get_contig_position(geno.paths[0].start_correct_pos()).second /*new_pos*/
                                   );

    haplotype.graph_complexity_to_stats();

    // Update the likelihood scores
    haplotype.explain_to_score(pn_index, non_unique_paths, geno.mapq, fully_aligned, mismatches);
    haplotype.coverage_to_gts(pn_index, geno.is_proper_pair()); // Update the coverage
  }
}


std::vector<HaplotypeCall>
VcfWriter::get_haplotype_calls() const
{
  std::vector<HaplotypeCall> hap_calls;
  hap_calls.reserve(haplotypes.size());

  for (auto const & hap : haplotypes)
  {
    HaplotypeCall new_hap_call(hap.get_haplotype_calls(), hap.gts);
    hap_calls.push_back(std::move(new_hap_call));
  }

  return hap_calls;
}


std::vector<uint32_t>
VcfWriter::explain_map_to_haplotype_scores(std::size_t const pn_index, ExplainMap const & explain_map)
{
  std::shared_ptr<ExplainMap> shared_explain_map = std::make_shared<ExplainMap>(explain_map);
  std::shared_ptr<std::size_t> shared_pn_index = std::make_shared<std::size_t>(pn_index);
  std::size_t const num = explain_map.cbegin()->second.size();
  std::vector<uint32_t> hap_scores((num + 1) * num / 2, 0u);
  auto hap_scores_split = paw::split(hap_scores, Options::instance()->threads);

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
    paw::Station hap_scores_station(Options::instance()->threads);
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

  paw::join(hap_scores, hap_scores_split);
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
    hap_score += haplotypes[it->first].best_score_of_a_path(pn_index, it->second[index.first], it->second[index.second]);

  return hap_score;
}


void
VcfWriter::find_path_explanation(GenotypePaths const & gt_path,
                                 std::vector<std::pair<uint32_t, std::bitset<MAX_NUMBER_OF_HAPLOTYPES> > > & ids_and_path_explain
                                 )
{
  // ids_and_path_explain is a vector of haplotype_id and haplotype_explain pairs
  assert(ids_and_path_explain.size() == 0);
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


} // namespace gyper
