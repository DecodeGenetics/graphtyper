#include <assert.h>
#include <bitset>
#include <iostream>
#include <map>

#include <boost/log/trivial.hpp>

#include <graphtyper/graph/absolute_position.hpp>
#include <graphtyper/graph/graph.hpp>
#include <graphtyper/typer/alignment.hpp>
#include <graphtyper/typer/genotype_paths.hpp>
#include <graphtyper/typer/segment_calling.hpp>
#include <graphtyper/typer/vcf_writer.hpp>
#include <graphtyper/typer/vcf.hpp>
#include <graphtyper/utilities/sam_reader.hpp>
#include <graphtyper/utilities/io.hpp>
#include <graphtyper/utilities/graph_help_functions.hpp>


namespace
{


void
print_explain_map(std::vector<std::string> const & hap_ids,
                  std::map<uint32_t, std::vector<std::bitset<gyper::MAX_NUMBER_OF_HAPLOTYPES> > > const & explain_map,
                  int32_t index = -1,
                  int32_t test_index = -1
                 )
{
  // Index == -1 means all indexes will be printed
  if (index == -1)
  {
    for (auto it = explain_map.begin(); it != explain_map.end(); ++it)
    {
      // Type of it->second is std::vector<std::bitset<gyper::MAX_NUMBER_OF_HAPLOTYPES> >
      std::cout.width(5);
      std::cout << it->first << ": ";

      for (auto const & explain_bitset : it->second)
      {
        // for (unsigned e = 0; e < explain_bitset.size(); ++e)
        // {
        //   if (explain_bitset.test(e))
        //   {
        //     std::cout << e;
        //     break;
        //   }
        //
        // }
        std::cout << explain_bitset.any();
      }
      std::cout << "\n";
    }
  }
  else if (test_index == -1)
  {
    std::cout << index << ": ";

    for (auto it = explain_map.begin(); it != explain_map.end(); ++it)
    {
      // Type of it->second is std::vector<std::bitset<gyper::MAX_NUMBER_OF_HAPLOTYPES> >
      std::cout << it->second[index].any() << " ";
    }

    std::cout << std::endl;
  }
  else
  {
    // std::cout << index << ": ";

    for (auto it = explain_map.begin(); it != explain_map.end(); ++it)
    {
      if (static_cast<int64_t>(it->first) == index)
      {
        std::cout << test_index << ": ";

        for (unsigned i = 0; i < it->second.size(); ++i)
        {
          if (it->second[i].test(test_index))
          {
            std::cout << hap_ids[i] << " ";
          }
        }

        std::cout << std::endl;
        break;
      }
      // Type of it->second is std::vector<std::bitset<gyper::MAX_NUMBER_OF_HAPLOTYPES> >
      // std::cout << it->second[index].any() << " ";
    }

    std::cout << "\n";
  }
}


void
insert_into_explain_map(std::map<uint32_t, std::vector<std::bitset<gyper::MAX_NUMBER_OF_HAPLOTYPES> > > & explain_map,
                        std::pair<uint32_t, std::bitset<gyper::MAX_NUMBER_OF_HAPLOTYPES> > const & var_explanation,
                        unsigned i,
                        std::size_t var_num
                        )
{
  auto find_it = explain_map.find(var_explanation.first);

  if (find_it == explain_map.end())
  {
    // Not found
    // assert (j = 0);
    // std::cout << "Was not found! i, j = " << i << "," << j << std::endl;
    // std::cout << "[caller] INFO: Inserting new variant " << var_explanation.first << std::endl;
    std::vector<std::bitset<gyper::MAX_NUMBER_OF_HAPLOTYPES> > new_vec(var_num);
    new_vec[i] = var_explanation.second;
    explain_map[var_explanation.first] = std::move(new_vec);
  }
  else
  {
    // Was found
    // assert (j > 0);
    // std::cout << "Was found! i, j = " << i << "," << j << std::endl;
    find_it->second[i] |= var_explanation.second;
  }
}


void
add_start_on_explain_map(std::map<uint32_t, std::vector<std::bitset<gyper::MAX_NUMBER_OF_HAPLOTYPES> > > & explain_map)
{
  std::vector<uint8_t> has_started(explain_map.begin()->second.size(), 0);

  for (auto it = explain_map.begin(); it != explain_map.end(); ++it)
  {
    // Type of it->second is std::vector<std::bitset<gyper::MAX_NUMBER_OF_HAPLOTYPES> >
    for (unsigned i = 0; i < it->second.size(); ++i)
    {
      assert(it->second.size() == explain_map.begin()->second.size());

      if (has_started[i])
      {
        continue;
      }
      else if (it->second[i].any())
      {
        has_started[i] = 1;
      }
      else
      {
        // Set all as true
        it->second[i].set();
      }
    }
  }
}


void
remove_insignificant_variants(std::map<uint32_t, std::vector<std::bitset<gyper::MAX_NUMBER_OF_HAPLOTYPES> > > & explain_map)
{
  for (auto it = explain_map.cbegin(); it != explain_map.cend(); )
  {
    // Type of it->second is std::vector<std::bitset<gyper::MAX_NUMBER_OF_HAPLOTYPES> >
    unsigned coverage = 0;

    for (auto explain_bitset : it->second)
    {
      if (explain_bitset.any())
      {
        ++coverage;
      }
    }

    double static const FILTER = 0.2;

    if (static_cast<double>(coverage) / static_cast<double>(it->second.size()) < FILTER)
    {
      // Remove if fraction of coverage is lower than FILTER
      explain_map.erase(it++);
    }
    else
    {
      ++it;
    }
  }
}


void
remove_out_of_order_variants(std::map<uint32_t, std::vector<std::bitset<gyper::MAX_NUMBER_OF_HAPLOTYPES> > > & exon_explain_map,
                             std::map<uint32_t, std::vector<std::bitset<gyper::MAX_NUMBER_OF_HAPLOTYPES> > > & intron_explain_map
                             )
{
  if (exon_explain_map.size() == 0)
    return;

  BOOST_LOG_TRIVIAL(info) << "[graphtyper::segment_calling] Removing out of order variants.";

  // Find all unique variants
  std::vector<uint32_t> uniq_variants;

  for (auto it = exon_explain_map.cbegin(); it != exon_explain_map.cend(); ++it)
  {
    auto find_it = std::find(uniq_variants.begin(), uniq_variants.end(), it->first);

    if (find_it == uniq_variants.end())
    {
      uniq_variants.push_back(it->first);
    }
  }

  for (auto it = intron_explain_map.cbegin(); it != intron_explain_map.cend(); ++it)
  {
    auto find_it = std::find(uniq_variants.begin(), uniq_variants.end(), it->first);

    if (find_it == uniq_variants.end())
    {
      uniq_variants.push_back(it->first);
    }
  }

  // Sort the unique variants
  std::sort(uniq_variants.begin(), uniq_variants.end());
  assert(uniq_variants.size() > 0);

  // Find the longest sequence of consecutive variants
  unsigned i = 0;
  uint32_t max_start_i = 0;
  uint32_t max_end_i = 0;

  while (i < uniq_variants.size() - 1)
  {
    uint32_t start_i = i;
    uint32_t end_i = i + 1;

    while (uniq_variants[i] + 1 == uniq_variants[i + 1] or uniq_variants[i] + 2 == uniq_variants[i + 1])
    {
      if (uniq_variants[i] + 2 == uniq_variants[i + 1])
      {
        ++end_i;
        ++i;
      }

      ++end_i;
      ++i;
    }

    if (end_i - start_i > max_end_i - max_start_i)
    {
      max_end_i = end_i;
      max_start_i = start_i;
    }

    ++i;
  }

  std::vector<uint32_t> longest_uniq_variants(uniq_variants.begin() + max_start_i, uniq_variants.begin() + max_end_i);

  for (auto it = exon_explain_map.cbegin(); it != exon_explain_map.cend(); )
  {
    // Type of it->second is std::vector<std::bitset<gyper::MAX_NUMBER_OF_HAPLOTYPES> >
    if (std::find(longest_uniq_variants.begin(), longest_uniq_variants.end(), it->first) == longest_uniq_variants.end())
    {
      // Remove the variant if it is not found in the longest unique variants vector
      BOOST_LOG_TRIVIAL(info) << "[graphtyper::segment_calling] Removing from exon " << it->first;
      exon_explain_map.erase(it++);
    }
    else
    {
      ++it;
    }
  }

  for (auto it = intron_explain_map.cbegin(); it != intron_explain_map.cend(); )
  {
    // Type of it->second is std::vector<std::bitset<gyper::MAX_NUMBER_OF_HAPLOTYPES> >
    if (std::find(longest_uniq_variants.begin(), longest_uniq_variants.end(), it->first) == longest_uniq_variants.end())
    {
      // Remove the variant if it is not found in the longest unique variants vector
      BOOST_LOG_TRIVIAL(info) << "[graphtyper::segment_calling] Removing from intron " << it->first;
      intron_explain_map.erase(it++);
    }
    else
    {
      ++it;
    }
  }
}


void
add_end_on_explain_map(std::map<uint32_t, std::vector<std::bitset<gyper::MAX_NUMBER_OF_HAPLOTYPES> > > & explain_map)
{
  std::vector<uint8_t> has_ended(explain_map.begin()->second.size(), 0);

  for (auto it = explain_map.rbegin(); it != explain_map.rend(); ++it)
  {
    // Type of it->second is std::vector<std::bitset<gyper::MAX_NUMBER_OF_HAPLOTYPES> >
    assert(it->second.size() == explain_map.begin()->second.size());

    for (unsigned i = 0; i < it->second.size(); ++i)
    {
      if (has_ended[i])
      {
        continue;
      }
      else if (it->second[i].any())
      {
        has_ended[i] = 1;
      }
      else
      {
        // Sets all as true
        it->second[i].set();
      }
    }
  }
}


std::size_t
determine_reference_index(std::map<uint32_t, std::vector<std::bitset<gyper::MAX_NUMBER_OF_HAPLOTYPES> > > const & exon_explain_map,
                          std::map<uint32_t, std::vector<std::bitset<gyper::MAX_NUMBER_OF_HAPLOTYPES> > > const & intron_explain_map
                          )
{
  BOOST_LOG_TRIVIAL(info) << "[graphtyper::segment_calling] Determining reference index from explain maps.";
  std::vector<uint32_t> ref_counts(intron_explain_map.begin()->second.size(), 0u);

  for (auto it = exon_explain_map.cbegin(); it != exon_explain_map.cend(); ++it)
  {
    // Type of it->second is std::vector<std::bitset<gyper::MAX_NUMBER_OF_HAPLOTYPES> >
    for (unsigned i = 0; i < it->second.size(); ++i)
    {
      if (it->second[i].test(0))
      {
        assert(i < ref_counts.size());
        ++ref_counts[i];
      }
    }
  }

  for (auto it = intron_explain_map.cbegin(); it != intron_explain_map.cend(); ++it)
  {
    // Type of it->second is std::vector<std::bitset<gyper::MAX_NUMBER_OF_HAPLOTYPES> >
    for (unsigned i = 0; i < it->second.size(); ++i)
    {
      if (it->second[i].test(0))
      {
        assert(i < ref_counts.size());
        ++ref_counts[i];
      }
    }
  }

  int64_t max_ref_counts = -1;
  std::size_t max_ref_counts_index = 0;

  for (std::size_t i = 0; i < ref_counts.size(); ++i)
  {
    if (ref_counts[i] > max_ref_counts)
    {
      max_ref_counts_index = i;
      max_ref_counts = ref_counts[i];
    }
  }

  if (max_ref_counts < static_cast<int64_t>(exon_explain_map.size()) + static_cast<int64_t>(intron_explain_map.size()))
  {
    BOOST_LOG_TRIVIAL(warning) << "[graphtyper::segment_calling] No path is purely reference. " << max_ref_counts << " out of "
                               << exon_explain_map.size() + intron_explain_map.size();
  }

  return max_ref_counts_index;
}


void
put_reference_in_front(std::map<uint32_t, std::vector<std::bitset<gyper::MAX_NUMBER_OF_HAPLOTYPES> > > & explain_map,
                       std::vector<std::string> & hap_ids,
                       std::size_t const ref_index,
                       bool const change_hap_ids
                       )
{
  assert(hap_ids.size() > ref_index);

  if (change_hap_ids)
  {
    std::string ref_str(hap_ids[ref_index]);

    hap_ids.erase(hap_ids.begin() + ref_index);
    hap_ids.insert(hap_ids.begin(), ref_str);
  }

  if (ref_index == 0)
  {
    return;
  }

  for (auto map_it = explain_map.begin(); map_it != explain_map.end(); ++map_it)
  {
    std::vector<std::bitset<gyper::MAX_NUMBER_OF_HAPLOTYPES> > explain_vec(map_it->second);
    std::bitset<gyper::MAX_NUMBER_OF_HAPLOTYPES> ref_bitset(explain_vec[ref_index]);
    explain_vec.erase(explain_vec.begin() + ref_index);
    explain_vec.insert(explain_vec.begin(), ref_bitset);
    map_it->second = std::move(explain_vec);
  }
}


} // anon namespace


namespace gyper
{


void
segment_calling(std::vector<std::string> const & segment_fasta_files,
                VcfWriter & writer,
                std::string const & segment_path,
                std::vector<std::string> const & samples
                )
{
  assert (samples.size() > 0);
  BOOST_LOG_TRIVIAL(info) << "[graphtyper::segment_calling] Segment VCF is at " << segment_path;
  Vcf segment_vcf(WRITE_BGZF_MODE, segment_path);

  for (auto const & sample : samples)
    segment_vcf.sample_names.push_back(sample);

  // Update all maximum log scores
  for (auto & haplototype : writer.haplotypes)
    haplototype.update_max_log_score();

  BOOST_LOG_TRIVIAL(info) << "[graphtyper::segment_calling] Gathering segments from " << segment_fasta_files.size() << " segments.";
  std::vector<Segment> segments;
  using THapPaths = std::vector<GenotypePaths>;
  std::vector<std::map<std::string, THapPaths> > all_haplotype_paths; // haplotype ID to a all its genotype paths results
  std::vector<std::vector<uint8_t> > has_long_exon; // haplotype ID to a all its genotype paths results

  {
    for (auto seg_it = segment_fasta_files.cbegin(); seg_it != segment_fasta_files.cend(); ++seg_it)
    {
      // Type of *seg_it is std::string, it is the fasta filename of the current segment
      std::map<std::string, std::vector<seqan::Dna5String> > mhc_hap = read_haplotypes_from_fasta(*seg_it);
      std::map<std::string, THapPaths> haplotype_paths;

      for (auto hap_it = mhc_hap.cbegin(); hap_it != mhc_hap.cend(); ++hap_it)
      {
        std::cout << "ID " << hap_it->first << std::endl;
        haplotype_paths[hap_it->first] = find_haplotype_paths(hap_it->second);
      }

      std::vector<uint8_t> gene_has_long_exons;

      for (unsigned i = 0; i < mhc_hap.begin()->second.size(); ++i)
      {
        // if (i % 2 == 1 && seqan::length(mhc_hap.begin()->second[i]) >= 2 * K)
        if (i % 2 == 1 && i < 10) // Only check exons 1-4
          gene_has_long_exons.push_back(1);
        else
          gene_has_long_exons.push_back(0);
      }

      has_long_exon.push_back(std::move(gene_has_long_exons));
      all_haplotype_paths.push_back(std::move(haplotype_paths));
    }
  }

  BOOST_LOG_TRIVIAL(info) << "[graphtyper::segment_calling] Iterating paths of segments. ";

  for (auto haplotype_paths_it = all_haplotype_paths.cbegin(); haplotype_paths_it != all_haplotype_paths.cend(); ++haplotype_paths_it)
  {
    // Type of haplotype_paths_it is std::map<std::string, THapPaths>::iterator
    std::vector<std::string> hap_ids;
    using TExplainMap = std::map<uint32_t, std::vector<std::bitset<MAX_NUMBER_OF_HAPLOTYPES> > >;
    TExplainMap exon_explain_map;
    TExplainMap intron_explain_map;

    {
      int i = 0;

      for (auto it = haplotype_paths_it->begin(); it != haplotype_paths_it->end(); ++i, ++it)
      {
        // Type of it->first is std::string
        // Type of it->second is std::vector<std::vector<GenotypePaths> >
        // if (it->second.size() > 0)
        //   std::cout << "[graphtyper::segment_calling] Name = " << it->first << std::endl;
        // std::cout << "[graphtyper::segment_calling] INFO: Number of genotype paths = " << it->second.size() << std::endl;

        std::size_t const & k = std::distance(all_haplotype_paths.cbegin(), haplotype_paths_it);
        std::vector<std::vector<std::pair<uint32_t, std::bitset<MAX_NUMBER_OF_HAPLOTYPES> > > > path_explanations; // Previous path explanation

        for (unsigned j = 0; j < it->second.size(); ++j)
        {
          std::vector<std::pair<uint32_t, std::bitset<MAX_NUMBER_OF_HAPLOTYPES> > > path_explanation;
          GenotypePaths const & path = it->second[j];

          // std::cout << "[graphtyper::segment_calling] INFO: Index " << j << " gt path size = " << it->second[0].paths.size() << std::endl;
          //
          // if (it->second[j].paths.size() == 0)
          //   continue;

          // assert (it->second[j].paths.size() > 0);

          std::cout << "[graphtyper::segment_calling] INFO: " << it->first << ", index " << j << std::endl;

          if (it->second[j].paths.size() == 0)
          {
            std::cout << "NO PATHS" << std::endl;
          }
          else if (it->second[j].paths.size() == 1)
          {
            std::cout << "UNIQUE PATH "
                      << absolute_pos.get_contig_position(it->second[j].paths[0].start_ref_reach_pos()).second << "-"
                      << absolute_pos.get_contig_position(it->second[j].paths[0].end_ref_reach_pos()).second << std::endl;
          }
          else if (it->second[j].paths.size() > 1)
          {
            for (auto const & dup_path : it->second[j].paths)
            {
              std::cout << "DUPLICATED PATH "
                        << absolute_pos.get_contig_position(dup_path.start_ref_reach_pos()).second << "-"
                        << absolute_pos.get_contig_position(dup_path.end_ref_reach_pos()).second << std::endl;
            }
          }

          writer.find_path_explanation(path, path_explanation);
          path_explanations.push_back(std::move(path_explanation));
        }

        // Add to explain maps
        for (unsigned j = 0; j < path_explanations.size(); ++j)
        {
          assert (it->second.size() == path_explanations.size());

          for (unsigned p = 0; p < path_explanations[j].size(); ++p)
          {
            auto & var_explanation = path_explanations[j][p];

            // Type of var_explanation is std::pair<uint32_t, std::bitset<MAX_NUMBER_OF_HAPLOTYPES> >
            assert(k < has_long_exon.size());
            assert(j < has_long_exon[k].size());
            // std::cout << "[graphtyper::segment_calling] Var explain " << var_explanation.first << " ";
            // std::cout << var_explanation.second.count() << " ";
            //
            // for (unsigned i = 0; i < 100; ++i)
            //   std::cout << var_explanation.second.test(i);
            //
            // std::cout << std::endl;

            if (has_long_exon[k][j])
            {
              insert_into_explain_map(exon_explain_map, var_explanation, i, haplotype_paths_it->size());
            }
            else
            {
              insert_into_explain_map(intron_explain_map, var_explanation, i, haplotype_paths_it->size());
            }
          }
        }

        hap_ids.push_back(it->first);
      }
    }

    // Resize all exon explain maps
    for (auto & e : exon_explain_map)
      e.second.resize(hap_ids.size());

    // Resize all intron explain maps
    for (auto & i : intron_explain_map)
      i.second.resize(hap_ids.size());

    BOOST_LOG_TRIVIAL(info) << "[graphtyper::segment_calling] Done creating explain maps.";

    // DEBUG
    // print_explain_map(exon_explain_map);
    // std::cout << std::endl;
    // print_explain_map(intron_explain_map);
    // DEBUG ENDS HERE

    // Remove variants which only a small portion overlaps
    if (intron_explain_map.size() == 0)
    {
      BOOST_LOG_TRIVIAL(error) << "Could not align any introns to the graph. Did you align to the correct graph?";
      std::exit(1);
    }

    remove_out_of_order_variants(exon_explain_map, intron_explain_map);
    remove_insignificant_variants(exon_explain_map);
    remove_insignificant_variants(intron_explain_map);

    BOOST_LOG_TRIVIAL(info) << "[graphtyper::segment_calling] Done removing out of order and insignificant variants";

    // This condition is required to avoid segfault!
    if (exon_explain_map.size() > 0)
    {
      BOOST_LOG_TRIVIAL(info) << "[graphtyper::segment_calling] exon_explain_map sequences.size() = " << exon_explain_map.begin()->second.size();
      add_start_on_explain_map(intron_explain_map);
      add_end_on_explain_map(intron_explain_map);
    }
    else
    {
      BOOST_LOG_TRIVIAL(info) << "[graphtyper::segment_calling] exon_explain_map sequences is empty";
    }

    BOOST_LOG_TRIVIAL(info) << "[graphtyper::segment_calling] intron_explain_map sequences.size() = " << intron_explain_map.begin()->second.size();

    // DEBUG
    if (exon_explain_map.size() > 0)
    {
      print_explain_map(hap_ids, exon_explain_map, 19, 405);
      std::cout << std::endl;
    }

    // if (intron_explain_map.size() > 0)
    //   print_explain_map(intron_explain_map);
    // DEBUG ENDS HERE

    // std::cout << "[caller] Overall number of haplotypes before removing is " << hap_ids.size() << std::endl;
    // remove_non_existing_alleles(hap_ids, exon_explain_map, intron_explain_map);
    // std::cout << "[caller] Overall number of haplotypes after removing is " << hap_ids.size() << std::endl;

    std::size_t ref_index = determine_reference_index(exon_explain_map, intron_explain_map);
    put_reference_in_front(exon_explain_map, hap_ids, ref_index, false); // Last parameter is change hap_ids
    put_reference_in_front(intron_explain_map, hap_ids, ref_index, true);

    BOOST_LOG_TRIVIAL(info) << "[graphtyper::segment_calling] Generating scores with reference " << hap_ids[0];

    // Create segment for these results
    BOOST_LOG_TRIVIAL(info) << "[graphtyper::segment_calling] Creating a new segment.";
    assert (haplotype_paths_it->size() > 0);
    assert (haplotype_paths_it->begin()->second.size() > 0);
    unsigned s = 0;
    unsigned e = 0;

    while ((haplotype_paths_it->begin()->second.begin() + s)->longest_paths().size() == 0
           )
    {
      ++s;

      if (s == haplotype_paths_it->begin()->second.size())
      {
        std::cerr << "warning: could not find a segment which matched" << std::endl;
        --s;
        break;
      }

      assert (s < haplotype_paths_it->begin()->second.size());
    }

    while ((haplotype_paths_it->begin()->second.rbegin() + e)->longest_paths().size() == 0 &&
           e != haplotype_paths_it->begin()->second.size() - 1
           )
    {
      ++e;
    }

    assert ((haplotype_paths_it->begin()->second.begin() + s)->longest_paths().size() > 0);
    assert ((haplotype_paths_it->begin()->second.begin() + e)->longest_paths().size() > 0);
    Path const longest_path_start = (haplotype_paths_it->begin()->second.begin() + s)->longest_paths().front();
    Path const longest_path_end = (haplotype_paths_it->begin()->second.rbegin() + e)->longest_paths().front();
    int64_t seq_size = static_cast<int64_t>(longest_path_end.end_correct_pos()) - static_cast<int64_t>(longest_path_start.start_correct_pos()) + 1;
    uint32_t segment_start = longest_path_start.start_correct_pos();
    BOOST_LOG_TRIVIAL(info) << "[graphtyper::segment_calling] Segment sequence size = " << seq_size;

    if (seq_size < 0)
    {
      seq_size = static_cast<int64_t>(longest_path_start.end_correct_pos()) - static_cast<int64_t>(longest_path_end.start_correct_pos()) + 1;
      segment_start = longest_path_end.start_correct_pos();
    }

    BOOST_LOG_TRIVIAL(info) << "[graphtyper::segment_calling] Number of hap_ids is " << hap_ids.size();

    Segment seg(segment_start, static_cast<uint32_t>(seq_size), hap_ids);
    std::vector<std::vector<uint32_t> > hap_scores(samples.size());
    // Segment created

    // Check the score of all the exons, if there are any exons (if we only have full sequences without features, we assume all the sequences are introns)
    if (exon_explain_map.size() > 0)
    {
      for (uint32_t s = 0; s < samples.size(); ++s)
      {
        // std::cout << "derp function starting" << std::endl;
        std::vector<uint32_t> hap_score = writer.explain_map_to_haplotype_scores(s, exon_explain_map);
        // std::cout << "derp function done" << std::endl;
        assert(hap_score.size() > 0);
        auto max_score_it = std::max_element(hap_score.begin(), hap_score.end());
        uint32_t const max_score = *max_score_it;
        BOOST_LOG_TRIVIAL(info) << "[graphtyper::segment_calling] Highest exon score is " << max_score;
        // std::cout << "Best index should be " << to_pair(std::distance(hap_score.begin(), max_score_it)).first
        //           << "/" << to_pair(std::distance(hap_score.begin(), max_score_it)).second << std::endl;
        std::vector<std::pair<uint32_t, uint32_t> > best_indexes;

        for (uint32_t i = 0; i < hap_score.size(); ++i)
        {
          if (hap_score[i] >= max_score)
            best_indexes.push_back(to_pair(i));
        }

        assert(best_indexes.size() > 0);
        BOOST_LOG_TRIVIAL(info) << "[graphtyper::segment_calling] Number of best indexes are " << best_indexes.size();

        if (best_indexes.size() <= 100)
        {
          for (auto const & best_index : best_indexes)
            BOOST_LOG_TRIVIAL(info) << "[graphtyper::segment_calling] Best alleles: " << hap_ids[best_index.first] << "/" << hap_ids[best_index.second];

          if (best_indexes.size() > 1)
          {
            BOOST_LOG_TRIVIAL(info) << "[graphtyper::segment_calling] The best exon score is not unique, but there are less than 100 best exon scores";

            // Add intron scores
            std::vector<uint32_t> intron_scores(best_indexes.size(), 0);
            uint32_t max_intron_score = 0;
            uint32_t second_max_intron_score = 0;

            for (unsigned i = 0; i < best_indexes.size(); ++i)
            {
              intron_scores[i] = writer.explain_map_specific_indexes_to_haplotype_scores(s, best_indexes[i], intron_explain_map);

              if (intron_scores[i] > max_intron_score)
              {
                second_max_intron_score = max_intron_score;
                max_intron_score = intron_scores[i];
              }
              else if (intron_scores[i] != max_intron_score && intron_scores[i] > second_max_intron_score)
              {
                // Also check if we found a larger second largest score
                second_max_intron_score = intron_scores[i];
              }
            }

            assert (second_max_intron_score <= max_intron_score);
            BOOST_LOG_TRIVIAL(info) << "[graphtyper::segment_calling] Max and second max intron score is "
                                    << max_intron_score << ", " << second_max_intron_score;

            // Increase scores of alleles with the most likely introns
            if (max_intron_score > 0)
            {
              for (unsigned i = 0; i < intron_scores.size(); ++i)
              {
                assert (intron_scores[i] <= max_intron_score);

                if (intron_scores[i] == max_intron_score)
                {
                  BOOST_LOG_TRIVIAL(info) << "[graphtyper::segment_calling] Increasing scores of " << hap_ids[best_indexes[i].first] << "/"
                                          << hap_ids[best_indexes[i].second];
                  hap_score[to_index(best_indexes[i].first, best_indexes[i].second)] += std::max(10u, (max_intron_score - second_max_intron_score) / 2);
                }
              }
            }
          }
          else if (best_indexes.size() == 1)
          {
            BOOST_LOG_TRIVIAL(info) << "[graphtyper::segment_calling] Unique best exon score";
          }
          else
          {
            BOOST_LOG_TRIVIAL(info) << "[graphtyper::segment_calling] No best exon score";
          }
        }
        else
        {
          BOOST_LOG_TRIVIAL(info) << "[graphtyper::segment_calling] There are more than 100 best exon scores";
        }

        // std::cout << "Inserting scores..." << std::endl;
        seg.insert_score(hap_score);
        // std::cout << "Done inserting scores" << std::endl;
      }
    }
    else
    {
      for (uint32_t s = 0; s < samples.size(); ++s)
      {
        assert(intron_explain_map.size() > 0);
        std::vector<uint32_t> hap_score = writer.explain_map_to_haplotype_scores(s, intron_explain_map);
        assert (hap_score.size() > 0);
        auto max_it = std::max_element(hap_score.begin(), hap_score.end());
        BOOST_LOG_TRIVIAL(info) << "[graphtyper::segment_calling] Highest total score is " << *max_it;

        for (unsigned k = 0; k < hap_score.size(); ++k)
        {
          if (hap_score[k] == *max_it)
          {
            std::pair<uint16_t, uint16_t> calls =  gyper::to_pair(k);
            BOOST_LOG_TRIVIAL(info) << "[graphtyper::segment_calling] Call with highest total score is " << hap_ids[calls.first] << "/" << hap_ids[calls.second];
          }
        }

        seg.insert_score(hap_score);
      }
    }

    segment_vcf.add_segment(std::move(seg));
  }

  if (segment_vcf.segments.size() > 0)
  {
    segment_vcf.open_for_writing();
    segment_vcf.write_header();
    segment_vcf.write_segments();
    segment_vcf.close_vcf_file();
  }
}


} // namespace gyper
