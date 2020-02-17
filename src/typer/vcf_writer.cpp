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


namespace
{

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

#ifndef NDEBUG
  if (gyper::Options::instance()->hq_reads)
  {
    if (!fully_aligned || geno.paths[0].size() < 90 || mismatch_ratio > 0.025)
      return false;
  }
#endif // NDEBUG

  return true;
}


} // anon namespace


namespace gyper
{

VcfWriter::VcfWriter(uint32_t variant_distance)
{
  haplotypes = gyper::graph.get_all_haplotypes(variant_distance);
  BOOST_LOG_TRIVIAL(debug) << "[graphtyper::vcf_writer] Number of variant nodes in graph "
                           << graph.var_nodes.size();
  BOOST_LOG_TRIVIAL(debug) << "[graphtyper::vcf_writer] Got "
                           << haplotypes.size()
                           << " haplotypes.";
}


void
VcfWriter::set_samples(std::vector<std::string> const & samples)
{
  pns = std::vector<std::string>(samples);
  long const NUM_SAMPLES = pns.size();
  assert(NUM_SAMPLES > 0);

  // Insert items in the id to haplotype map
  for (long i = 0; i < static_cast<long>(haplotypes.size()); ++i)
  {
    auto & haplotype = haplotypes[i];
    haplotype.clear_and_resize_samples(NUM_SAMPLES);
    std::vector<uint32_t> gt_ids = haplotype.get_genotype_ids();

    for (long j = 0; j < static_cast<long>(gt_ids.size()); ++j)
    {
      id2hap[gt_ids[j]] =
        std::make_pair<uint32_t, uint32_t>(static_cast<uint32_t>(i), static_cast<uint32_t>(j));
    }
  }
}


void
VcfWriter::update_haplotype_scores_geno(GenotypePaths & geno, long const pn_index)
{
  if (are_genotype_paths_good(geno))
  {
    push_to_haplotype_scores(geno, pn_index);

#ifndef NDEBUG
    if (Options::instance()->stats.size() > 0)
      update_statistics(geno, pn_index);
  }
#else
  }
#endif //NDEBUG
}


#ifndef NDEBUG
void
VcfWriter::print_variant_group_details() const
{
  assert(pns.size() > 0);
  std::string const hap_stats_fn = Options::instance()->stats + "/" + pns[0] + "_variant_group_details.tsv.gz";
  BOOST_LOG_TRIVIAL(debug) << "[graphtyper::vcf_writer] Generating variant group info statistics to " << hap_stats_fn;
  std::stringstream hap_file;

  // Write header file
  hap_file << "groupID\talleleID\tcontig\tposition\tsequence\n";

  for (unsigned ps = 0; ps < haplotypes.size(); ++ps)
  {
    auto const & hap = haplotypes[ps];
    uint32_t const cnum = hap.get_genotype_num();

    for (uint32_t c = 0; c < cnum; ++c)
    {
      assert(hap.gts.size() > 0);
      uint32_t const abs_pos = hap.gts[0].id;
      std::vector<char> seq = graph.get_sequence_of_a_haplotype_call(hap.gts, c);
      assert(seq.size() > 1);
      auto contig_pos = absolute_pos.get_contig_position(abs_pos, gyper::graph);

      hap_file << ps << "\t" << c << "\t"
               << contig_pos.first << "\t" << contig_pos.second << "\t"
               << std::string(seq.begin() + 1, seq.end())
               << "\n";
    }
  }

  //std::lock_guard<std::mutex> lock(io_mutex);
  write_gzipped_to_file(hap_file, hap_stats_fn);
}


void
VcfWriter::print_variant_details() const
{
  assert(pns.size() > 0);
  std::string const variant_details_fn =
    Options::instance()->stats + "/" + pns[0] + "_variant_details.tsv.gz";

  BOOST_LOG_TRIVIAL(debug) << "[graphtyper::vcf_writer] Generating variant info statistics to "
                           << variant_details_fn;

  std::stringstream variant_file;

  // Write header file
  variant_file << "variantID\tcontig\tposition\tallele_num\tsequence\tSV\n";

  for (long v = 0; v < static_cast<long>(graph.var_nodes.size()); ++v)
  {
    long sv_id = -1; // -1 means not an SV
    auto const & label = graph.var_nodes[v].get_label();
    auto contig_pos = absolute_pos.get_contig_position(label.order, gyper::graph);
    auto const & seq = label.dna;
    auto find_it = std::find(seq.cbegin(), seq.cend(), '<');

    if (find_it != seq.cend() && std::distance(find_it, seq.cend()) > 11)
    {
      // It is an SV
      std::istringstream ss{std::string(find_it + 4, find_it + 11)};
      ss >> sv_id;

      // If we can't parse correctly the SV ID we ignore it
      assert(ss.eof());
      assert(sv_id < static_cast<long>(graph.SVs.size()));
    }

    variant_file << v << '\t'
                 << contig_pos.first << '\t'
                 << contig_pos.second << '\t'
                 << label.variant_num << '\t'
                 << std::string(label.dna.begin(), label.dna.end()) << '\t';

    if (sv_id == -1)
      variant_file << ".";
    else
      variant_file << graph.SVs[sv_id].get_allele_with_model();

    variant_file << '\n';
  }

  //std::lock_guard<std::mutex> lock(io_mutex);
  write_gzipped_to_file(variant_file, variant_details_fn);
}


void
VcfWriter::print_statistics_headers() const
{
  assert(pns.size() > 0);
  assert(Options::instance()->stats.size() > 0);

  // Get filenames
  std::string const read_details_fn =
    Options::instance()->stats + "/" + pns[0] + "_read_details.tsv.gz";

  std::string const path_details_fn =
    Options::instance()->stats + "/" + pns[0] + "_read_path_details.tsv.gz";

  std::stringstream read_ss;
  std::stringstream path_ss;
  read_ss << "query\t" << "sample\t" << "read\t" << "read_qual\t"
          << "alignment_length\t" << "mapping_quality\t" << "original_mapped_pos\t"
          << "ml_insert_size\n";

  path_ss << "query\t" << "pathID\t" << "read_start_index\t" << "read_end_index\t"
          << "num_mismatches\t" << "strand\t" << "contig\t" << "alignment_begin\t"
          << "alignment_end\t" << "overlapping_variant_nodes\n";

  //std::lock_guard<std::mutex> lock(io_mutex);
  write_gzipped_to_file(read_ss, read_details_fn, false /*append*/);
  write_gzipped_to_file(path_ss, path_details_fn, false /*append*/);
}


void
VcfWriter::print_geno_statistics(std::stringstream & read_ss,
                                 std::stringstream & path_ss,
                                 GenotypePaths const & geno,
                                 long pn_index
                                 )
{
  std::stringstream id;
  assert(geno.details);
  id << pns[pn_index] << "_" << geno.details->query_name << "/" << ((geno.flags & IS_FIRST_IN_PAIR) != 0 ? 1 : 2);

  read_ss << id.str() << "\t"
          << pns[pn_index] << "\t"
          << std::string(geno.read2.begin(), geno.read2.end()) << "\t"
          << std::string(geno.qual2.begin(), geno.qual2.end()) << "\t"
          << geno.longest_path_length << "\t"
          << geno.original_pos << "\t";

  if (geno.ml_insert_size == 0x7FFFFFFFl)
    read_ss << ".";
  else
    read_ss << geno.ml_insert_size;

  read_ss << "\n";

  for (std::size_t p = 0; p < geno.paths.size(); ++p)
  {
    auto const & path = geno.paths[p];
    uint32_t const ref_reach_start = path.start_ref_reach_pos();
    uint32_t const ref_reach_end = path.end_ref_reach_pos();

    auto const contig_pos_start = absolute_pos.get_contig_position(ref_reach_start, gyper::graph);
    auto const contig_pos_end = absolute_pos.get_contig_position(ref_reach_end, gyper::graph);

    std::vector<std::size_t> overlapping_vars;

    for (long i = 0; i < static_cast<long>(path.var_order.size()); ++i)
    {
      auto const & var_order = path.var_order[i];
      auto const & num = path.nums[i];

      auto find_it = id2hap.find(var_order);
      assert(find_it->second.first < haplotypes.size());
      auto const & hap = haplotypes[find_it->second.first];
      assert(find_it->second.second < hap.gts.size());
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
            << (((geno.flags & IS_SEQ_REVERSED) == 0) ? "F" : "B") << "\t"
            << contig_pos_start.first << "\t"
            << contig_pos_start.second << "\t"
            << contig_pos_end.second << "\t";

    std::sort(overlapping_vars.begin(), overlapping_vars.end());

    if (overlapping_vars.size() == 0)
    {
      path_ss << ".";
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
VcfWriter::update_statistics(GenotypePaths const & geno, long const pn_index)
{
  // Get filenames
  std::string const read_details_fn =
    Options::instance()->stats + "/" + pns[0] + "_read_details.tsv.gz";

  std::string const path_details_fn =
    Options::instance()->stats + "/" + pns[0] + "_read_path_details.tsv.gz";

  // Write reads and paths to these streams and then write those in a gzipped file
  std::stringstream read_ss;
  std::stringstream path_ss;

  print_geno_statistics(read_ss, path_ss, geno, pn_index);

  //std::lock_guard<std::mutex> lock(io_mutex);
  write_gzipped_to_file(read_ss, read_details_fn, true /*append*/);
  write_gzipped_to_file(path_ss, path_details_fn, true /*append*/);
}


void
VcfWriter::update_statistics(std::vector<GenotypePaths> & genos, long const pn_index)
{
  // Get filenames
  std::string const read_details_fn =
    Options::instance()->stats + "/" + pns[0] + "_read_details.tsv.gz";

  std::string const path_details_fn =
    Options::instance()->stats + "/" + pns[0] + "_read_path_details.tsv.gz";

  // Write reads and paths to these streams and then write those in a gzipped file
  std::stringstream read_ss;
  std::stringstream path_ss;

  for (auto const & geno : genos)
    print_geno_statistics(read_ss, path_ss, geno, pn_index);

  //std::lock_guard<std::mutex> lock(io_mutex);
  write_gzipped_to_file(read_ss, read_details_fn, true /*append*/);
  write_gzipped_to_file(path_ss, path_details_fn, true /*append*/);
}


void
VcfWriter::update_statistics(std::vector<std::pair<GenotypePaths, GenotypePaths> > & genos,
                             long pn_index
                             )
{
  // Get filenames
  std::string const read_details_fn =
    Options::instance()->stats + "/" + pns[0] + "_read_details.tsv.gz";

  std::string const path_details_fn =
    Options::instance()->stats + "/" + pns[0] + "_read_path_details.tsv.gz";

  // Write reads and paths to these streams and then write those in a gzipped file
  std::stringstream read_ss;
  std::stringstream path_ss;

  for (auto const & geno : genos)
  {
    print_geno_statistics(read_ss, path_ss, geno.first, pn_index);
    print_geno_statistics(read_ss, path_ss, geno.second, pn_index);
  }

  //std::lock_guard<std::mutex> lock(io_mutex);
  write_gzipped_to_file(read_ss, read_details_fn, true /*append*/);
  write_gzipped_to_file(path_ss, path_details_fn, true /*append*/);
}


#endif // NDEBUG

void
VcfWriter::push_to_haplotype_scores(GenotypePaths & geno, long const pn_index)
{
  assert(are_genotype_paths_good(geno));

  // Quality metrics
  bool const fully_aligned = geno.all_paths_fully_aligned();
  bool const non_unique_paths = !geno.all_paths_unique();
  std::size_t const mismatches = geno.paths[0].mismatches;
  bool has_low_quality_snp = false;

  std::unordered_map<uint32_t, bool> recent_ids; // maps to is_overlapping

  for (auto p_it = geno.paths.begin(); p_it != geno.paths.end(); ++p_it)
  {
    assert(p_it->var_order.size() == p_it->nums.size());

    for (long i = 0; i < static_cast<long>(p_it->var_order.size()); ++i)
    {
      assert(id2hap.count(p_it->var_order[i]) == 1);
      std::pair<uint32_t, uint32_t> const type_ids = id2hap[p_it->var_order[i]]; // hap_id = first, gen_id = second

      assert(type_ids.first < haplotypes.size());
      assert(type_ids.second < haplotypes[type_ids.first].gts.size());
      assert(p_it->nums[i].any());

      auto & hap = haplotypes[type_ids.first];
      auto & num = p_it->nums[i];

      long constexpr MIN_OFFSET = 3;
      bool is_overlapping = (p_it->start_ref_reach_pos() + MIN_OFFSET) <= p_it->var_order[i] &&
                            (p_it->end_ref_reach_pos() - MIN_OFFSET) > p_it->var_order[i];
      recent_ids[type_ids.first] |= is_overlapping;

      if (!has_low_quality_snp && graph.is_snp(hap.gts[type_ids.second]))
      {
        long const offset = p_it->var_order[i] - p_it->start_correct_pos();

        if (offset < static_cast<long>(geno.qual2.size()))
        {
          assert(geno.qual2[offset] >= 33);
          has_low_quality_snp = (geno.qual2[offset] - 33) < 25;
        }
      }

      // Add explanation
      hap.add_explanation(type_ids.second, num);

      // Add coverage if the explanation is unique
      if (num.count() == 1)
      {
        // Check which bit is set
        uint16_t b = 0;

        while (not num.test(b))
          ++b;

        hap.add_coverage(type_ids.second, b);
      }
      else /* Otherwise set the coverage is ambiguous */
      {
        hap.add_coverage(type_ids.second, 1);

        if (num.test(0))
          hap.add_coverage(type_ids.second, 0);
        else
          hap.add_coverage(type_ids.second, 2);
      }
    }
  }

  // After each read, move the "explain" to the score vector.
  for (auto it = recent_ids.begin(); it != recent_ids.end(); ++it)
  {
    assert(it->first < haplotypes.size());
    assert(pn_index < static_cast<long>(haplotypes[it->first].hap_samples.size()));
    assert(geno.paths.size() > 0);

    auto & haplotype = haplotypes[it->first];

    // Move INFO to variant statistics. This needs to called before 'coverage_to_gts', because it uses the coverage of each variant.
    haplotype.clipped_reads_to_stats(fully_aligned);
    haplotype.mapq_to_stats(geno.mapq);
    haplotype.strand_to_stats(geno.flags);
    //haplotype.realignment_to_stats(geno.original_pos /*original_pos*/,
    //                               absolute_pos.get_contig_position(geno.paths[0].start_correct_pos()).second /*new_pos*/);

    // Update the likelihood scores
    haplotype.explain_to_score(pn_index,
                               non_unique_paths,
                               geno.flags,
                               fully_aligned,
                               it->second, // is read overlapping
                               has_low_quality_snp,
                               mismatches);

    // Update the coverage values
    haplotype.coverage_to_gts(pn_index, geno.is_proper_pair());
  }
}


std::vector<HaplotypeCall>
VcfWriter::get_haplotype_calls() const
{
  std::vector<HaplotypeCall> hap_calls;
  hap_calls.reserve(haplotypes.size());

  for (Haplotype const & hap : haplotypes)
  {
    HaplotypeCall new_hap_call(hap); //.get_haplotype_calls(), hap.gts);
    hap_calls.push_back(std::move(new_hap_call));
  }

  return hap_calls;
}


} // namespace gyper
