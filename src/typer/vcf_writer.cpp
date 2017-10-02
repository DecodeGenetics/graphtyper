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
#include <graphtyper/utilities/io.hpp>
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

    if (geno.paths[0].mismatches > 4)
      return false;

    //// Any path overlapping a variant must also not have too many mismatches
    //for (auto const & path : geno.paths)
    //{
    //  if (path.var_order.size() > 0 && path.mismatches > 2)
    //    return false;
    //}
  }

  return true;
}


} // anon namespace


namespace gyper
{

VcfWriter::VcfWriter(std::vector<std::string> const & samples, uint32_t variant_distance)
  : pns(samples)
  , pn(samples[0])
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

  if (Options::instance()->stats.size() > 0)
  {
    BOOST_LOG_TRIVIAL(info) << "[graphtyper::vcf_writer] Gathering read statistics";
    this->read_stats = std::unique_ptr<ReadStats>(new ReadStats(samples.size()));
  }
}


void
VcfWriter::update_haplotype_scores_from_paths(std::vector<GenotypePaths> & genos,
                                              std::size_t const pn_index
  )
{
  {
    std::lock_guard<std::mutex> lock(haplotype_mutex);

    for (auto & geno : genos)
    {
      if (are_genotype_paths_good(geno))
        update_haplotype_scores_from_path(geno, pn_index, 0 /*unpaired*/);
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

  if (Options::instance()->stats.size() > 0)
    update_statistics(genos, pn_index);
}


//void
//VcfWriter::update_statistics(GenotypePaths & geno, std::size_t const pn_index, unsigned const read_pair)
//{
//  // Update read statistics
//  if (geno.ml_insert_size != 0x7FFFFFFFl)
//  {
//    assert (read_stats);
//    assert (pn_index < read_stats->insert_sizes.size());
//    assert (0 < geno.paths.size());
//    assert (pn_index < read_stats->mismatches.size());
//    assert (geno.read_pair_mismatches != static_cast<uint8_t>(255u));
//
//    // Only log positive values, so each pair only registered once
//    if (geno.ml_insert_size > 0)
//    {
//      read_stats->insert_sizes[pn_index].push_back(static_cast<uint32_t>(geno.ml_insert_size));
//      read_stats->mismatches[pn_index].push_back(static_cast<uint8_t>(geno.read_pair_mismatches));
//    }
//  }
//
//  // Write read details to file
//  {
//    std::string const read_details_fn = Options::instance()->stats + "/reads_details.tsv.gz";
//    std::stringstream read_details;
//    read_details << geno.query_name << "\t"
//                 << std::string(geno.read.begin(), geno.read.end()) << "\t"
//                 << std::string(geno.qual.begin(), geno.qual.end()) << "\t"
//                 << geno.longest_path_length << "\t"
//                 << static_cast<uint16_t>(geno.mapq) << "\t"
//                 << geno.original_pos << "\t"
//                 << (geno.ml_insert_size == 0x7FFFFFFFl ? -1 : geno.ml_insert_size) << "\t"
//                 << static_cast<uint16_t>(geno.read_pair_mismatches) << "\t"
//                 << geno.is_first_in_pair << "\t"
//                 << geno.forward_strand << "\t"
//                 << geno.is_originally_unaligned << "\t"
//                 << geno.is_originally_clipped << "\t"
//                 << "\n";
//
//    write_gzipped_to_file(read_details, read_details_fn, true /*append*/);
//  }
//
//  // Insert all paths
//  //for (auto const & path : geno.paths)
//  //{
//  //
//  //}
//
//  // Select a path "randomly" and update coverage
//  {
//    assert (read_stats);
//    assert (geno.paths.size() > 0);
//    uint32_t const RND = 3 * geno.paths[0].start + 5 * geno.paths[0].end + 7 * geno.paths[0].mismatches;
//    auto const & path = geno.paths.at(RND % geno.paths.size());
//
//    // Check if the path is purely reference
//    if (path.nums.size () > 0)
//    {
//      bool is_pure_reference = true;
//      bool is_pure_alternative = true;
//
//      for (auto const & num : path.nums)
//      {
//        bool const REF_COV = num.test(0);
//        bool const ALT_COV = (num.count() - num.test(0)) > 0;
//
//        if (!REF_COV)
//        {
//          is_pure_reference = false;
//
//          if (!is_pure_alternative)
//            break;
//        }
//
//        if (!ALT_COV)
//        {
//          is_pure_alternative = false;
//
//          if (!is_pure_reference)
//            break;
//        }
//      }
//
//      if (is_pure_reference ^ is_pure_alternative)
//      {
//        if (is_pure_reference)
//          read_stats->ref_read_abs_pos[pn_index].push_back({path.start_correct_pos(), path.end_correct_pos()});
//        else
//          read_stats->alt_read_abs_pos[pn_index].push_back({path.start_correct_pos(), path.end_correct_pos()});
//      }
//      else
//      {
//        // Pick randomly
//        if ((RND + 11 * geno.paths[0].start) % 2 == 0)
//          read_stats->ref_read_abs_pos[pn_index].push_back({path.start_correct_pos(), path.end_correct_pos()});
//        else
//          read_stats->alt_read_abs_pos[pn_index].push_back({path.start_correct_pos(), path.end_correct_pos()});
//      }
//    }
//    else
//    {
//      // The path does not overlap any variants
//      read_stats->other_read_abs_pos[pn_index].push_back({path.start_correct_pos(), path.end_correct_pos()});
//    }
//  }
//
//  // If the alignment is unique, log the sequence successor and linked haplotypes
//  if (Options::instance()->haplotype_statistics && geno.paths.size() == 1)
//  {
//    auto const & path = geno.paths[0];
//
//    for (unsigned i = 0; i < path.var_order.size(); ++i)
//    {
//      if (path.nums[i].count() != 1)
//        continue;
//
//      uint32_t b = 0;
//
//      while (!path.nums[i].test(b))
//      {
//        ++b;
//        assert (b < path.nums[i].size());
//      }
//
//      uint32_t v = 0;
//
//      while (graph.var_nodes[v].get_label().order != path.var_order[i] || graph.var_nodes[v].get_label().variant_num != b)
//      {
//        ++v;
//        assert (v < graph.var_nodes.size());
//      }
//
//      uint32_t const reach = graph.var_nodes[v].get_label().reach();
//      assert(id2hap.count(path.var_order[i]) == 1);
//      std::pair<uint32_t, uint32_t> const type_ids = id2hap.at(path.var_order[i]); //hap_id = first, gen_id = second
//      assert(type_ids.first < haplotypes.size());
//      assert(pn_index < haplotypes[type_ids.first].hap_samples.size());
//      assert(haplotypes[type_ids.first].hap_samples[pn_index].stats);
//      auto const & stats = haplotypes[type_ids.first].hap_samples[pn_index].stats;
//
//      // Read pair information
//      if (read_pair == 1)
//      {
//        stats->pair_stats.push_back({type_ids.first, b}); // Add this pair of haplotype ID and genotype ID
//      }
//      else if (read_pair == 2)
//      {
//        // Check all other haplotypes
//        for (auto & hap : haplotypes)
//        {
//          auto const & stats_first = hap.hap_samples[pn_index].stats;
//
//          // Add read pair information to stats
//          for (auto const & first_read : stats_first->pair_stats)
//          {
//            assert (first_read.gt_id < haplotypes[first_read.hap_id].hap_samples[pn_index].stats->pair_info.size());
//            assert (b < haplotypes[type_ids.first].hap_samples[pn_index].stats->pair_info.size());
//            haplotypes[first_read.hap_id].hap_samples[pn_index].stats->pair_info[first_read.gt_id].push_back({type_ids.first, b});
//            haplotypes[type_ids.first].hap_samples[pn_index].stats->pair_info[b].push_back({first_read.hap_id, first_read.gt_id});
//          }
//        }
//      }
//
//      uint32_t constexpr minW = 1;
//      uint32_t constexpr W = 1;
//
//      if (path.end_correct_pos() >= reach + minW && path.end_correct_pos() - reach < geno.read.size())
//      {
//        assert (b < stats->successor.size());
//        auto start_it = geno.read.end() - path.end_correct_pos() + reach;
//        assert (std::distance(start_it, geno.read.end()) >= minW);
//
//        if (std::distance(start_it, geno.read.end()) > W)
//          stats->successor[b].push_back(std::vector<char>(start_it, start_it + W));
//        else
//          stats->successor[b].push_back(std::vector<char>(start_it, geno.read.end()));
//
//        auto qual_start_it = geno.qual.end() - path.end_correct_pos() + reach;
//
//        // Add it again if all bases are high quality
//        if (std::count(qual_start_it, qual_start_it + minW, 33) == 0)
//        {
//          if (std::distance(start_it, geno.read.end()) > W)
//            stats->successor[b].push_back(std::vector<char>(start_it, start_it + W));
//          else
//            stats->successor[b].push_back(std::vector<char>(start_it, geno.read.end()));
//        }
//      }
//
//      // Code for predecessor
//      // if (graph.var_nodes[v].get_label().order >= path.start_correct_pos() + minW &&
//      //     graph.var_nodes[v].get_label().order - path.start_correct_pos() < geno.read.size()
//      //    )
//      // {
//      //   assert (b < stats->predecessor.size());
//      //   auto end_it = geno.read.begin() + graph.var_nodes[v].get_label().order - path.start_correct_pos();
//      //   auto qual_end_it = geno.qual.begin() + graph.var_nodes[v].get_label().order - path.start_correct_pos();
//      //
//      //   if (std::count(qual_end_it - minW, qual_end_it, 33) == 0)
//      //   {
//      //     if (std::distance(geno.read.begin(), end_it) > W)
//      //       stats->predecessor[b].push_back(std::vector<char>(end_it - W, end_it));
//      //     else
//      //       stats->predecessor[b].push_back(std::vector<char>(geno.read.begin(), end_it));
//      //   }
//      // }
//    }
//  }
//}


void
VcfWriter::print_haplotype_details() const
{
  std::string const hap_stats_fn = Options::instance()->stats + "/" + pns[0] + "_haplotype_details.tsv.gz";
  BOOST_LOG_TRIVIAL(info) << "[graphtyper::vcf_writer] Generating haplotype info statistics to " << hap_stats_fn;
  std::stringstream hap_file;

  // Write header file
  hap_file << "#haplotypeID\talleleID\tcontig\tposition\tsequence\n";

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
               << std::string(seq.begin() + 1, seq.end());
      hap_file << "\n";
     }
  }

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
  variant_file << "#variantID\tcontig\tposition\tallele_num\tsequence\n";

  for (std::size_t v = 0; v < graph.var_nodes.size(); ++v)
  {
    auto const & label = graph.var_nodes[v].get_label();
    auto contig_pos = absolute_pos.get_contig_position(label.order);

    variant_file << v << "\t"
                 << contig_pos.first << "\t"
                 << contig_pos.second << "\t"
//                 << label.order << "\t"
                 << label.variant_num << "\t"
                 << std::string(label.dna.begin(), label.dna.end())
                 << "\n";
  }

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
  read_ss << "#query\t" << "read_group\t" << "sample\t" << "read\t" << "read_qual\t"
          << "alignment_length\t" << "mapping_quality\t" << "original_mapped_pos\t"
          << "ml_insert_size\t" << "is_originally_unaligned\t"  << "is_originally_clipped" << "\n";

  path_ss << "#query\t" << "pathID\t" << "read_start_index\t" << "read_end_index\t"
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
  id << geno.details->query_name << "/" << (geno.is_first_in_pair ? 1 : 2);

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
                             std::size_t const pn_index
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
                                             std::size_t const pn_index,
                                             unsigned const /*read_pair*/
  )
{
  assert (are_genotype_paths_good(geno));

  // Quality metrics
  bool const fully_aligned = geno.all_paths_fully_aligned();
  bool const non_unique_paths = !geno.all_paths_unique();
  std::size_t const mismatches = geno.paths[0].mismatches;
  assert (geno.read.size() >= 63);

  //if (Options::instance()->stats.size() > 0) // Update statistics if '--stats' was passed
  //  this->update_statistics(geno, pn_index, read_pair);

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
VcfWriter::generate_statistics()
{
  if (Options::instance()->stats.size() == 0 || pns.size() == 0)
    return;

  // Read statistics
  if (read_stats && read_stats->mismatches.size() > 0)
  {
    assert (read_stats->mismatches.size() == read_stats->insert_sizes.size());

    for (std::size_t p = 0; p < pns.size(); ++p)
    {
      std::string const read_stats_fn = Options::instance()->stats + "/" + pns[p] + "_reads_stats.tsv.gz";
      BOOST_LOG_TRIVIAL(info) << "[graphtyper::vcf_writer] Generating read statistics to " << read_stats_fn;
      std::stringstream read_file;
      read_file << "#InsertSize (IS)\n";

      // Options
      uint32_t constexpr MIN_REPORTED_INSERT_SIZE = 100;
      uint32_t constexpr MAX_REPORTED_INSERT_SIZE = 500;
      uint8_t constexpr MAX_REPORTED_MISMATCHES = 6;

      // Insert sizes with different number of mismatches
      {
        std::unordered_map<uint64_t, uint32_t> insert_size_mismatches_counts;
        std::unordered_map<uint32_t, uint32_t> insert_size_counts;
        std::unordered_map<uint8_t, uint32_t> mismatches_counts;

        for (std::size_t j = 0; j < read_stats->insert_sizes[p].size(); ++j)
        {
          uint32_t insert_size = read_stats->insert_sizes[p][j];
          uint8_t mismatches = read_stats->mismatches[p][j];

          if (insert_size <= 0)
            continue;

          insert_size = std::max(MIN_REPORTED_INSERT_SIZE, static_cast<uint32_t>(insert_size));
          insert_size = std::min(MAX_REPORTED_INSERT_SIZE, static_cast<uint32_t>(insert_size));
          mismatches = std::min(MAX_REPORTED_MISMATCHES, static_cast<uint8_t>(mismatches));

          ++insert_size_counts[insert_size];
          ++mismatches_counts[mismatches];
          ++insert_size_mismatches_counts[static_cast<uint64_t>(insert_size) | (static_cast<uint64_t>(mismatches) << 32)];
        }

        // Insert size 0 should be half to only count each read pair once
        read_file << "IS\t" << MIN_REPORTED_INSERT_SIZE << '\t' << (insert_size_counts[MIN_REPORTED_INSERT_SIZE]/2);

        for (uint32_t is = MIN_REPORTED_INSERT_SIZE + 1; is <= MAX_REPORTED_INSERT_SIZE; ++is)
          read_file << "\nIS\t" << is << '\t' << insert_size_counts[is];

        read_file << '\n';
        read_file << "#MisMatches (MM)\n";
        read_file << "MM\t0\t" << mismatches_counts[0];

        for (uint32_t m = 1; m <= MAX_REPORTED_MISMATCHES; ++m)
          read_file << '\n' << "MM\t" << m << '\t' << mismatches_counts[m];

        read_file << '\n';

        // Combine insert size and mismatch in a metric called read likelihood
        read_file << "#Insert size with different number of read pair Mismatches (IM)\n";
        read_file << "IM\t" << MIN_REPORTED_INSERT_SIZE;

        for (uint64_t mismatches = 0; mismatches <= MAX_REPORTED_MISMATCHES; ++mismatches)
          read_file << '\t' << insert_size_mismatches_counts[mismatches << 32];

        for (uint64_t is = MIN_REPORTED_INSERT_SIZE + 1; is <= MAX_REPORTED_INSERT_SIZE; ++is)
        {
          read_file << "\nIM\t" << is;

          for (uint64_t m = 0; m <= MAX_REPORTED_MISMATCHES; ++m)
          {
            uint64_t const key = is | (m << 32);
            read_file << '\t' << insert_size_mismatches_counts[key];
          }
        }
      }

      read_file << '\n';

      // Coverage
      {
        using Tread_cov_map = std::map<std::pair<std::string, uint32_t>, uint32_t>;
        uint32_t constexpr INTERVAL_SIZE = 5;

        auto add_rounded_pos_lambda = [INTERVAL_SIZE]
          (Tread_cov_map & read_cov_map, std::pair<uint32_t, uint32_t> const abs_pos_range)
        {
          std::pair<std::string, uint32_t> contig_pos1 = absolute_pos.get_contig_position(abs_pos_range.first);
          contig_pos1.second -= contig_pos1.second % INTERVAL_SIZE;
          std::pair<std::string, uint32_t> contig_pos2 = absolute_pos.get_contig_position(abs_pos_range.second);
          contig_pos2.second -= contig_pos2.second % INTERVAL_SIZE;

          // Make sure both are on the same chromosome
          if (contig_pos1.first == contig_pos2.first)
          {
            while (contig_pos1.second <= contig_pos2.second)
            {
              ++read_cov_map[contig_pos1];
              contig_pos1.second += INTERVAL_SIZE;
            }
          }
        };

        read_file << "#Coverage of interval (COV). Fields show reference allele coverage,"
                  << " alternative allele coverage, other coverage and total coverage.\n";

        // Map of ref coverage by a range of position
        Tread_cov_map ref_coverage;

        for (auto const abs_pos : read_stats->ref_read_abs_pos[p])
          add_rounded_pos_lambda(ref_coverage, abs_pos);

        // Map of ref coverage by a range of position
        Tread_cov_map alt_coverage;

        for (auto const abs_pos : read_stats->alt_read_abs_pos[p])
          add_rounded_pos_lambda(alt_coverage, abs_pos);

        // Map of other coverage by a range of position
        Tread_cov_map other_coverage;

        for (auto const abs_pos : read_stats->other_read_abs_pos[p])
          add_rounded_pos_lambda(other_coverage, abs_pos);

        auto a_it = alt_coverage.begin();
        auto r_it = ref_coverage.begin();

        auto print_to_read_file_lambda = [&](Tread_cov_map::const_iterator it, uint32_t const ref_count, uint32_t const alt_count)
        {
          auto find_it = other_coverage.find(it->first);
          uint32_t const other_count = find_it == other_coverage.end() ? 0 : find_it->second;

          read_file << "COV\t"
                    << it->first.first << "\t["
                    << it->first.second << "-"
                    << (it->first.second + INTERVAL_SIZE - 1) << "]\t"
                    << ref_count << "\t"
                    << alt_count << "\t"
                    << other_count << "\t"
                    << (ref_count + alt_count + other_count) << "\n";
        };

        while (a_it != alt_coverage.end() || r_it != ref_coverage.end())
        {
          // Check if either one is at the end
          if (a_it == alt_coverage.end())
          {
            assert (r_it != ref_coverage.end());
            print_to_read_file_lambda(r_it, r_it->second, 0);
            ++r_it;
            continue;
          }

          if (r_it == ref_coverage.end())
          {
            assert (a_it != alt_coverage.end());
            print_to_read_file_lambda(a_it, 0, a_it->second);
            ++a_it;
            continue;
          }

          // Both are not at the end
          if (r_it->first < a_it->first)
          {
            print_to_read_file_lambda(r_it, r_it->second, 0);
            ++r_it;
          }
          else if (a_it->first < r_it->first)
          {
            print_to_read_file_lambda(a_it, 0, a_it->second);
            ++a_it;
          }
          else
          {
            print_to_read_file_lambda(r_it, r_it->second, a_it->second);
            ++r_it;
            ++a_it;
          }
        }
      } // Coverage ends

      write_gzipped_to_file(read_file, read_stats_fn.c_str());
    }
  }
  else
  {
    BOOST_LOG_TRIVIAL(info) << "[graphtyper::vcf_writer] Skipped generating read statistics";
  }

  // Haplotypes
  {
    std::string const hap_stats_fn = Options::instance()->stats + "/" + pns[0] + "_haplotypes_info.tsv.gz";
    BOOST_LOG_TRIVIAL(info) << "[graphtyper::vcf_writer] Generating haplotype info statistics to " << hap_stats_fn;
    std::stringstream hap_file;

    // Write header file
    hap_file << "#haplotypeID\talleleID\tcontig\tpos\tsequence\n";

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
        hap_file << ps << "\t" << c << "\t" << contig_pos.first << "\t" << contig_pos.second << "\t"
                 << std::string(seq.begin() + 1, seq.end());
        hap_file << "\n";
      }
    }

    write_gzipped_to_file(hap_file, hap_stats_fn);
  }

  // Haplotype statistics
  {
    std::string const hap_stats_fn = Options::instance()->stats + "/" + pns[0] + "_haplotypes_stats.tsv.gz";
    BOOST_LOG_TRIVIAL(info) << "[graphtyper::vcf_writer] Generating haplotype statistics to " << hap_stats_fn;
    std::stringstream hap_file;

    // Write header file
    hap_file << "#sample\thaplotypeID\talleleID\tcoverage\tunique_coverage\n";

    // Loop samples
    for (unsigned i = 0; i < pns.size(); ++i)
    {
      auto const & pn = pns[i];

      // Haplotype statistics
      for (unsigned ps = 0; ps < haplotypes.size(); ++ps)
      {
        auto const & hap = haplotypes[ps];
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

          for (auto const & gt : hap.gts)
          {
            q /= gt.num;
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

          /*
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
          */

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

    write_gzipped_to_file(hap_file, hap_stats_fn);
    //std::ofstream compressed(hap_stats_fn.c_str(), std::ofstream::binary);
    //boost::iostreams::filtering_streambuf<boost::iostreams::input> out;
    //out.push(boost::iostreams::gzip_compressor());
    //out.push(hap_file);
    //boost::iostreams::copy(out, compressed);
  }
}


} // namespace gyper
