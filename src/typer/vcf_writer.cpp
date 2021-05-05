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
#include <graphtyper/typer/primers.hpp>
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

  if (gyper::Options::const_instance()->hq_reads)
  {
    if (!fully_aligned || geno.paths[0].size() < 90 || mismatch_ratio > 0.035)
      return false;
  }

  return true;
}


} // anon namespace


namespace gyper
{

VcfWriter::VcfWriter(uint32_t variant_distance)
{
  haplotypes = gyper::graph.get_all_haplotypes(variant_distance);
  BOOST_LOG_TRIVIAL(debug) << __HERE__ << " Number of variant nodes in graph " << graph.var_nodes.size();
  BOOST_LOG_TRIVIAL(debug) << __HERE__ << " Got " << haplotypes.size() << " haplotypes.";
}


void
VcfWriter::set_samples(std::vector<std::string> const & samples)
{
  pns = samples;
  long const NUM_SAMPLES = pns.size();
  assert(NUM_SAMPLES > 0);

  // Insert items in the id to haplotype map
  for (long i = 0; i < static_cast<long>(haplotypes.size()); ++i)
  {
    auto & haplotype = haplotypes[i];
    haplotype.clear_and_resize_samples(NUM_SAMPLES);
    id2hap[haplotype.gt.id] = static_cast<uint32_t>(i);
  }
}


void
VcfWriter::update_haplotype_scores_geno(GenotypePaths & geno, long const pn_index, Primers const * primers)
{
#ifndef NDEBUG
  if (Options::const_instance()->stats.size() > 0)
    update_statistics(geno, pn_index);
#endif

  if (Options::const_instance()->is_segment_calling)
    return;

  if (are_genotype_paths_good(geno))
  {
    if (primers)
      primers->check(geno);

#ifdef GT_DEV
    std::map<std::pair<uint16_t, uint16_t>, std::vector<std::pair<uint16_t, uint16_t> > > merged_connections;
    auto con1 = push_to_haplotype_scores(geno, pn_index);

    // add to merged connections
    for (auto it1 = con1.begin(); it1 != con1.end(); ++it1)
    {
      #ifndef NDEBUG
      for (auto const & hap_gt : it1->second)
      {
        assert(hap_gt.first > it1->first.first);
      }
      #endif

      merged_connections.insert(*it1);
    }

    // add new connections to hap_samples
    for (auto const & hap_cov_id1_2 : merged_connections)
    {
      auto const & hap_cov_id1 = hap_cov_id1_2.first;
      assert(hap_cov_id1.first < static_cast<long>(haplotypes.size()));
      Haplotype & hap1 = haplotypes[hap_cov_id1.first];
      assert(pn_index < static_cast<long>(hap1.hap_samples.size()));
      assert(hap_cov_id1.second < haplotypes[hap_cov_id1.first].hap_samples[pn_index].connections.size());
      auto & conn = hap1.hap_samples[pn_index].connections[hap_cov_id1.second];

      for (auto const & hap_cov_id2 : hap_cov_id1_2.second)
      {
        assert(hap_cov_id2.first > hap_cov_id1.first);
        Haplotype const & hap2 = haplotypes[hap_cov_id2.first];
        assert(hap_cov_id2.second < hap2.gt.num);

        auto insert_it = conn.insert({hap_cov_id2.first, std::vector<uint16_t>(hap2.gt.num)});
        assert(hap_cov_id2.second < insert_it.first->second.size());
        ++insert_it.first->second[hap_cov_id2.second];
      }
    }
#else
    push_to_haplotype_scores(geno, pn_index);
#endif // GT_DEV
  }
}


void
VcfWriter::update_haplotype_scores_geno(std::pair<GenotypePaths *, GenotypePaths *> & geno_paths,
                                        long const pn_index,
                                        Primers const * primers)
{
  assert(geno_paths.first);
  assert(geno_paths.second);

  auto & geno1 = *geno_paths.first;
  auto & geno2 = *geno_paths.second;

#ifndef NDEBUG
  if (Options::instance()->stats.size() > 0)
  {
    update_statistics(geno1, pn_index);
    update_statistics(geno2, pn_index);
  }
#endif

#ifdef GT_DEV
  std::map<std::pair<uint16_t, uint16_t>, std::vector<std::pair<uint16_t, uint16_t> > > con1;
  std::map<std::pair<uint16_t, uint16_t>, std::vector<std::pair<uint16_t, uint16_t> > > con2;
#endif

  bool const is_good1 = are_genotype_paths_good(geno1);
  bool const is_good2 = are_genotype_paths_good(geno2);

  if (Options::const_instance()->is_segment_calling && (!is_good1 || !is_good2))
    return;

  if (is_good1)
  {
    if (primers)
      primers->check(geno1);

#ifdef GT_DEV
    con1 = push_to_haplotype_scores(geno1, pn_index);
#else
    push_to_haplotype_scores(geno1, pn_index);
#endif // GT_DEV
  }

  if (is_good2)
  {
    if (primers)
      primers->check(geno2);

#ifdef GT_DEV
    con2 = push_to_haplotype_scores(geno2, pn_index);
#else
    push_to_haplotype_scores(geno2, pn_index);
#endif // GT_DEV
  }

#ifdef GT_DEV
  std::map<std::pair<uint16_t, uint16_t>, std::vector<std::pair<uint16_t, uint16_t> > > merged_connections;

  if (con1.size() > 0 || con2.size() > 0)
  {
    // merge the connections, start with connections1
    for (auto it1 = con1.begin(); it1 != con1.end(); ++it1)
    {
      #ifndef NDEBUG
      for (auto const & hap_gt : it1->second)
      {
        assert(hap_gt.first > it1->first.first);
      }
      #endif

      auto insert_it = merged_connections.insert(*it1);
      assert(insert_it.second);

      for (auto it2 = con2.begin(); it2 != con2.end(); ++it2)
      {
        if (it2->first.first > it1->first.first)
          insert_it.first->second.push_back(it2->first);
      }
    }

    // moving on to connections2
    for (auto it2 = con2.begin(); it2 != con2.end(); ++it2)
    {
      auto insert_it = merged_connections.insert(*it2);

      if (!insert_it.second)
      {
        // Nothing was inserted
        std::copy(it2->second.begin(),
                  it2->second.end(),
                  std::back_inserter(insert_it.first->second));
      }

      for (auto it1 = con1.begin(); it1 != con1.end(); ++it1)
      {
        if (it1->first.first > it2->first.first)
          insert_it.first->second.push_back(it1->first);
      }
    }
  } // else do nothing, merged_connections will be empty as well

  // add new connections to hap_samples
  for (auto const & hap_cov_id1_2 : merged_connections)
  {
    auto const & hap_cov_id1 = hap_cov_id1_2.first;
    assert(hap_cov_id1.first < static_cast<long>(haplotypes.size()));
    Haplotype & hap1 = haplotypes[hap_cov_id1.first];
    assert(pn_index < static_cast<long>(hap1.hap_samples.size()));
    assert(hap_cov_id1.second < haplotypes[hap_cov_id1.first].hap_samples[pn_index].connections.size());
    auto & conn = hap1.hap_samples[pn_index].connections[hap_cov_id1.second];

    for (auto const & hap_cov_id2 : hap_cov_id1_2.second)
    {
      assert(hap_cov_id2.first > hap_cov_id1.first);
      Haplotype const & hap2 = haplotypes[hap_cov_id2.first];
      assert(hap_cov_id2.second < hap2.gt.num);

      auto insert_it = conn.insert({hap_cov_id2.first, std::vector<uint16_t>(hap2.gt.num)});
      assert(hap_cov_id2.second < insert_it.first->second.size());
      ++insert_it.first->second[hap_cov_id2.second];
    }
  }
#endif // GT_DEV
}


#ifndef NDEBUG
/*
void
VcfWriter::print_variant_group_details() const
{
  assert(pns.size() > 0);
  std::string const hap_stats_fn = Options::instance()->stats + "/" + pns[0] + "_variant_group_details.tsv.gz";
  BOOST_LOG_TRIVIAL(debug) << __HERE__ << " Generating variant group info statistics to " << hap_stats_fn;
  std::stringstream hap_file;

  // Write header file
  hap_file << "groupID\talleleID\tcontig\tposition\tsequence\n";

  for (long ps{0}; ps < static_cast<long>(haplotypes.size()); ++ps)
  {
    auto const & hap = haplotypes[ps];
    uint32_t const cnum = hap.get_genotype_num();

    for (uint32_t c = 0; c < cnum; ++c)
    {
      assert(hap.gts.size() > 0);
      uint32_t const abs_pos = hap.gts[0].id;
      std::vector<char> seq = graph.get_sequence_of_a_haplotype_call(hap.gts, c);
      assert(seq.size() > 1);
      auto contig_pos = absolute_pos.get_contig_position(abs_pos, gyper::graph.contigs);

      hap_file << ps << "\t" << c << "\t"
               << contig_pos.first << "\t" << contig_pos.second << "\t"
               << std::string(seq.begin() + 1, seq.end())
               << "\n";
    }
  }

  write_gzipped_to_file(hap_file, hap_stats_fn);
}
*/


void
VcfWriter::print_variant_details() const
{
  assert(pns.size() > 0);
  std::string const variant_details_fn =
    Options::instance()->stats + "/" + pns[0] + "_variant_details.tsv.gz";

  BOOST_LOG_TRIVIAL(debug) << __HERE__ << " Generating variant info statistics to "
                           << variant_details_fn;

  std::stringstream variant_file;

  // Write header file
  variant_file << "variantID\tcontig\tposition\tallele_num\tsequence\tSV\n";

  for (long v = 0; v < static_cast<long>(graph.var_nodes.size()); ++v)
  {
    long sv_id = -1; // -1 means not an SV
    auto const & label = graph.var_nodes[v].get_label();
    auto contig_pos = absolute_pos.get_contig_position(label.order, gyper::graph.contigs);
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
                                 long pn_index)
{
  std::stringstream id;
  assert(geno.details);
  id << pns[pn_index] << "_" << geno.details->query_name << "/"
     << ((geno.flags & IS_FIRST_IN_PAIR) != 0 ? 1 : 2);

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

    auto const contig_pos_start = absolute_pos.get_contig_position(ref_reach_start,
                                                                   gyper::graph.contigs);

    auto const contig_pos_end = absolute_pos.get_contig_position(ref_reach_end,
                                                                 gyper::graph.contigs);

    std::vector<std::size_t> overlapping_vars;

    for (long i = 0; i < static_cast<long>(path.var_order.size()); ++i)
    {
      auto const & var_order = path.var_order[i];
      auto const & num = path.nums[i];

      auto find_it = id2hap.find(var_order);
      assert(find_it->second < haplotypes.size());
      auto const & hap = haplotypes[find_it->second];
      auto const & gt = hap.gt;

      for (uint32_t g = 0; g < gt.num; ++g)
      {
        if (num.count(g))
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
                             long pn_index)
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

#ifdef GT_DEV
std::map<std::pair<uint16_t, uint16_t>, std::vector<std::pair<uint16_t, uint16_t> > >
#else // GT_DEV
void
#endif
// GT_DEV
VcfWriter::push_to_haplotype_scores(GenotypePaths & geno, long const pn_index)
{
  assert(are_genotype_paths_good(geno));

  // Quality metrics
  assert(geno.longest_path_length <= geno.read_length);
  int const clipped_bp = geno.read_length - geno.longest_path_length;
  bool const fully_aligned = clipped_bp == 0;
  assert(fully_aligned == geno.all_paths_fully_aligned());

  bool const non_unique_paths = !geno.all_paths_unique();
  std::size_t const mismatches = geno.paths[0].mismatches;
  bool has_low_quality_snp{false};

#ifdef GT_DEV
  std::map<uint32_t, bool> recent_ids; // maps to is_overlapping
  std::map<std::pair<uint16_t, uint16_t>, std::vector<std::pair<uint16_t, uint16_t> > > new_connections;
#else
  std::unordered_map<uint32_t, bool> recent_ids; // maps to is_overlapping
#endif // GT_DEV

  for (auto p_it = geno.paths.begin(); p_it != geno.paths.end(); ++p_it)
  {
    assert(p_it->var_order.size() == p_it->nums.size());

    for (long i = 0; i < static_cast<long>(p_it->var_order.size()); ++i)
    {
#ifndef NDEBUG
      if (id2hap.count(p_it->var_order[i]) == 0)
      {
        BOOST_LOG_TRIVIAL(error) << __HERE__ << " Did not find var order=" << p_it->var_order[i];

        for (auto id_hap : id2hap)
        {
          BOOST_LOG_TRIVIAL(info) << id_hap.first << " " << id_hap.second << " " << 0;
        }
      }

      assert(id2hap.count(p_it->var_order[i]) == 1);
#endif // NDEBUG

      uint32_t const hap_id = id2hap[p_it->var_order[i]]; // hap_id = first, gen_id = second

      assert(hap_id < haplotypes.size());
      assert(p_it->nums[i].size() > 0);

      auto & hap = haplotypes[hap_id];
      auto & num = p_it->nums[i];

      long constexpr MIN_OFFSET{3};
      bool is_overlapping = (p_it->start_ref_reach_pos() + MIN_OFFSET) <= p_it->var_order[i] &&
                            (p_it->end_ref_reach_pos() - MIN_OFFSET) > p_it->var_order[i];
      recent_ids[hap_id] |= is_overlapping;

      if (!has_low_quality_snp && graph.is_snp(hap.gt))
      {
        long const offset = p_it->var_order[i] - p_it->start_correct_pos();

        if (offset < static_cast<long>(geno.qual2.size()))
        {
          assert(geno.qual2[offset] >= 33);
          has_low_quality_snp = (geno.qual2[offset] - 33) < 25;
        }
      }

      // Add explanation
      hap.explains.insert(num.begin(), num.end());

      // Add coverage if the explanation is unique
      if (num.size() == 1)
      {
        hap.add_coverage(0, *num.begin());
      }
      else /* Otherwise set the coverage is ambiguous */
      {
        hap.add_coverage(0, 1);

        if (num.count(0))
          hap.add_coverage(0, 0);
        else
          hap.add_coverage(0, 2);
      }
    }
  }

#ifdef GT_DEV
  // check connections
  for (auto it = recent_ids.begin(); it != recent_ids.end(); ++it)
  {
    auto & haplotype1 = haplotypes[it->first];
    long const hap1_explains_count = haplotype1.explains.size();

    if (hap1_explains_count == 0 || hap1_explains_count > 64)
      continue;

    long hap1_counter{0};

    for (long b1{0}; hap1_counter < hap1_explains_count; ++b1)
    {
      assert(b1 < static_cast<long>(haplotype1.gt.num));

      if (!haplotype1.explains.count(b1))
        continue;

      ++hap1_counter;
      auto insert_it = new_connections.insert({{it->first, b1}, std::vector<std::pair<uint16_t, uint16_t> >(0)});
      auto & conn = insert_it.first->second;

      for (auto it2 = std::next(it); it2 != recent_ids.end(); ++it2)
      {
        auto & haplotype2 = haplotypes[it2->first];
        long const hap2_explains_count = haplotype2.explains.size();

        if (hap2_explains_count == 0 || hap2_explains_count > 64)
          continue;

        long const weight = hap1_explains_count * hap2_explains_count;
        assert(weight > 0);

        long hap2_counter{0};

        for (long b2{0}; hap2_counter < hap2_explains_count; ++b2)
        {
          assert(b2 < static_cast<long>(haplotype2.gt.num));

          if (!haplotype2.explains.count(b2))
            continue;

          ++hap2_counter;
          long const repeat = weight >= 3 ? 6 / weight : 1; // unique read gets 6, others have a lower weight
          assert(it->first < it2->first);

          for (long r{0}; r < repeat; ++r)
            conn.push_back({it2->first, b2});
        }
      }
    }
  }
#endif // GT_DEV

  // After each read, move the "explain" to the score vector
  for (auto it = recent_ids.begin(); it != recent_ids.end(); ++it)
  {
    assert(it->first < haplotypes.size());
    assert(pn_index < static_cast<long>(haplotypes[it->first].hap_samples.size()));
    assert(geno.paths.size() > 0);

    auto & haplotype = haplotypes[it->first];

    // Move INFO to variant statistics.
    // This needs to called before 'coverage_to_gts', because it uses the coverage of each variant.
    haplotype.clipped_reads_to_stats(clipped_bp, geno.read_length);
    haplotype.mapq_to_stats(geno.mapq);
    haplotype.strand_to_stats(geno.flags);
    haplotype.mismatches_to_stats(mismatches, geno.read_length);
    haplotype.score_diff_to_stats(geno.score_diff);
    // TODO display these new stats

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

    // Reset coverage and clear explains bitset
    haplotype.coverage = Haplotype::NO_COVERAGE;
    haplotype.explains.clear();
  }

#ifdef GT_DEV
  return new_connections;

#endif // GT_DEV
}


/*
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
*/


} // namespace gyper
