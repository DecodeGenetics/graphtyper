#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <unordered_set> // std::unordered_set

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/log/trivial.hpp>

#include <graphtyper/graph/graph.hpp>
#include <graphtyper/graph/graph_utils.hpp> // count_mismatches, count_mismatches_backward, add_node_dna_to_sequence
#include <graphtyper/graph/label.hpp>
#include <graphtyper/graph/node.hpp>
#include <graphtyper/graph/var_record.hpp>
#include <graphtyper/graph/sv.hpp>
#include <graphtyper/typer/path.hpp>
#include <graphtyper/typer/variant.hpp>
#include <graphtyper/utilities/type_conversions.hpp>
#include <graphtyper/utilities/options.hpp>


namespace gyper
{

Graph::Graph(bool _use_absolute_positions)
  : use_absolute_positions(_use_absolute_positions)
{}


void
Graph::clear()
{
  genomic_region.clear();
  reference.clear();
  ref_nodes.clear();
  var_nodes.clear();
  ref_reach_to_special_pos.clear();
  ref_reach_poses.clear();
  actual_poses.clear();

//  reference_offset = 0;
  use_absolute_positions = true;
}


void
Graph::add_genomic_region(std::vector<char> && reference_sequence,
                          std::vector<VarRecord> && var_records,
                          GenomicRegion && region)
{
//  region.region_to_refnode = static_cast<uint32_t>(ref_nodes.size());
  genomic_region = std::move(region);

  // Ignore alternative alleles with any 'N' or empty alleles
  for (auto & var : var_records)
  {
    var.alts.erase(
      std::remove_if(var.alts.begin(),
                     var.alts.end(),
                     [](std::vector<char> const & alt)
      {
        return std::find(alt.begin(), alt.end(), 'N') != alt.end() ||
        alt.size() == 0;
      }
                     ), var.alts.end());
  }

  // Remove records with no alternative alleles or N/* in reference allele
  var_records.erase(
    std::remove_if(var_records.begin(),
                   var_records.end(),
                   [this](VarRecord const & rec) -> bool
    {
      return std::find(rec.ref.begin(), rec.ref.end(), 'N') != rec.ref.end() ||
      std::find(rec.ref.begin(), rec.ref.end(), '*') != rec.ref.end() ||
      rec.alts.size() == 0 ||
      rec.pos < genomic_region.begin;
    }
                   ),
    var_records.end()
    );

  if (Options::instance()->add_all_variants)
  {
    BOOST_LOG_TRIVIAL(debug) << "[" << __HERE__ << "] Constructing graph of "
                             << var_records.size()
                             << " variants by finding all possible paths";

    for (long i = 0; i < static_cast<long>(var_records.size()); ++i)
    {
      // Check if the variation record overlaps the next one, and merge them if so
      while (i + 1l < static_cast<long>(var_records.size()) &&
             var_records[i + 1].pos < var_records[i].pos + var_records[i].ref.size()
             )
      {
        if (var_records[i].alts.size() * (var_records[i + 1].alts.size() + 1) >=
            (MAX_NUMBER_OF_HAPLOTYPES - 1))
        {
          // Take backup in case we will need to revert
          VarRecord backup_record_i(var_records[i]);
          VarRecord backup_record_i1(var_records[i + 1]);

          var_records[i + 1].merge(std::move(var_records[i]));

          // Check if there will be too many variants in this record
          if (var_records[i + 1].alts.size() >= (MAX_NUMBER_OF_HAPLOTYPES - 1))
          {
            BOOST_LOG_TRIVIAL(warning) << "[" << __HERE__ << "] Not all possible paths can be made "
                                       << "around this variant position: "
                                       << var_records[i].pos
                                       << " (REF = " << to_string(var_records[i].ref) << ")";

            var_records[i] = backup_record_i;
            var_records[i + 1] = backup_record_i1;

            if (var_records[i + 1].alts.size() + var_records[i].alts.size() <
                (MAX_NUMBER_OF_HAPLOTYPES - 1)
                )
            {
              var_records[i + 1].merge_one_path(std::move(var_records[i]));
            }
            else
            {
              BOOST_LOG_TRIVIAL(warning) << "[" << __HERE__ << "] Could not even add a single path "
                                         << "for that variant!"
                                         << var_records[i].pos
                                         << " (REF = "
                                         << to_string(var_records[i].ref) << ")";
            }
          }
        }
        else
        {
          var_records[i + 1].merge(std::move(var_records[i]));
        }

        if (var_records[i + 1].alts.size() >= (MAX_NUMBER_OF_HAPLOTYPES - 1))
        {
          BOOST_LOG_TRIVIAL(error) << "[" << __HERE__ << "] Found a variant with too many alleles! ("
                                   << (var_records[i + 1].alts.size() + 1)
                                   << ">="
                                   << MAX_NUMBER_OF_HAPLOTYPES << ')';
          std::exit(1);
        }

        var_records[i].clear();
        ++i;
      }
    }
  }
  else
  {
    BOOST_LOG_TRIVIAL(debug) << __HERE__ << " Constructing graph of "
                             << var_records.size()
                             << " variants.";

    for (long i = 0; i < static_cast<long>(var_records.size()); ++i)
    {
      // Check if the variation record overlaps the next one, and merge them if so
      while (i + 1 < static_cast<long>(var_records.size()) &&
             var_records[i + 1].pos < var_records[i].pos + var_records[i].ref.size())
      {
        if (var_records[i].is_sv && var_records[i + 1].is_sv)
        {
          // merge two SVs
          var_records[i + 1].merge_one_path(std::move(var_records[i]));
          var_records[i].clear();
          ++i;
          continue;
        }
        else if (var_records[i].is_sv)
        {
          // Delete next variant if it overlaps an SV breakpoint
          var_records[i + 1] = std::move(var_records[i]);
          var_records[i].clear();
          ++i;
          continue;
        }
        else if (var_records[i + 1].is_sv)
        {
          // Delete previous variant if it overlaps an SV breakpoint
          var_records[i].clear();
          ++i;
          continue;
        }

        if (var_records[i].alts.size() > 100 || (var_records[i + 1].pos - var_records[i].pos) < 4)
        {
          var_records[i + 1].merge_one_path(std::move(var_records[i]));
          var_records[i].clear();
          ++i;
          continue;
        }

        {
          // In a few extreme scenarios we cannot merge only one path here.
          //BOOST_LOG_TRIVIAL(fatal) << "i    : " << var_records[i].to_string();
          //BOOST_LOG_TRIVIAL(fatal) << "i+1  : " << var_records[i + 1].to_string();
          var_records[i + 1].merge(std::move(var_records[i]), 4); // last parameter is EXTRA_SUFFIX
          //BOOST_LOG_TRIVIAL(fatal) << "after: " << var_records[i + 1].to_string();
        }

        if (var_records[i + 1].alts.size() >= (MAX_NUMBER_OF_HAPLOTYPES - 1))
        {
          BOOST_LOG_TRIVIAL(warning) << "[" << __HERE__ << "] Found a variant with too many alleles.";
          var_records[i + 1].alts.resize(MAX_NUMBER_OF_HAPLOTYPES - 2);
        }

        var_records[i].clear(); // Clear variants that have been merged into others
        ++i;
      } // while
    }

    BOOST_LOG_TRIVIAL(debug) << __HERE__ << " Done merging overlapping variant records. Now I have "
                             << var_records.size()
                             << " variants.";
  }

  // Erase alternatives which are identical to the reference sequence
  for (auto & var_record : var_records)
  {
    var_record.alts.erase(std::remove_if(var_record.alts.begin(),
                                         var_record.alts.end(),
                                         [&](std::vector<char> const & alt)
      {
        return alt == var_record.ref;
      }
                                         ), var_record.alts.end());
  }


  // Erase records with no alternatives
  var_records.erase(std::remove_if(var_records.begin(),
                                   var_records.end(),
                                   [](VarRecord const & rec)
    {
      return rec.alts.size() == 0;
    }
                                   ), var_records.end());

  // Remove common suffix
  for (auto & var_record : var_records)
  {
    std::vector<char> const suffix = var_record.get_common_suffix();

    if (suffix.size() > 0)
    {
      BOOST_LOG_TRIVIAL(debug) << __HERE__ << " Removing commong suffix.";
      assert(suffix.size() < var_record.ref.size());
      var_record.ref.erase(var_record.ref.begin() + (var_record.ref.size() - suffix.size()), var_record.ref.end());
      assert(var_record.ref.size() > 0);

      for (auto & alt : var_record.alts)
      {
        alt.erase(alt.begin() + (alt.size() - suffix.size()), alt.end());
        assert(alt.size() > 0);
      }
    }
  }

  BOOST_LOG_TRIVIAL(debug) << __HERE__ << " Adding reference and variants to the graph.";
  bool is_continue{true};

  for (long i = 0; i < static_cast<long>(var_records.size()); ++i)
  {
    BOOST_LOG_TRIVIAL(debug) << __HERE__ << " i = " << i << ", pos = " << var_records[i].pos;
    unsigned const num_var = static_cast<unsigned>(var_records[i].alts.size()) + 1u;
    is_continue = add_reference(var_records[i].pos,
                                num_var,
                                reference_sequence);    // force in the first iteration

    if (is_continue)
    {
      BOOST_LOG_TRIVIAL(debug) << __HERE__ << " adding variant";
      add_variants(std::move(var_records[i]));
    }
    else
    {
      break;
    }
  }

  if (is_continue)
  {
    BOOST_LOG_TRIVIAL(debug) << __HERE__ << " Adding final reference sequence behind the last variant.";
    add_reference(static_cast<uint32_t>(reference_sequence.size()) + genomic_region.begin, 0, reference_sequence);
  }

  // If we chose to use absolute positions we need to change all labels
  if (use_absolute_positions)
  {
    BOOST_LOG_TRIVIAL(debug) << __HERE__ << " Updating labels.";
    uint32_t const offset = genomic_region.get_absolute_position(1);
    unsigned r = 0;
    assert(r < ref_nodes.size());

    while (ref_nodes[r].out_degree() != 0)
    {
      ref_nodes[r].change_label_order(offset);

      for (auto v : ref_nodes[r].get_vars())
      {
        var_nodes[v].change_label_order(offset);
      }

      ++r;
    }

    ref_nodes[r].change_label_order(offset);
  }

  // Keep the reference_sequence
  reference = std::move(reference_sequence);

  // Set offset
//  reference_offset = genomic_region.begin;
}


uint16_t
Graph::get_variant_num(uint32_t v) const
{
  assert(v < var_nodes.size());
  return static_cast<uint16_t>(v - ref_nodes[var_nodes[v].get_out_ref_index() - 1].get_var_index(0));
}


uint32_t
Graph::get_variant_order(long variant_id) const
{
  assert(variant_id < static_cast<long>(var_nodes.size()));
  return var_nodes[variant_id].get_label().order;
}


std::vector<char>
Graph::get_all_ref() const
{
  if (ref_nodes.size() == 0)
    return std::vector<char>(0);

  std::vector<char> ref(0);
  unsigned r = 0;
  unsigned v = 0;

  while (ref_nodes[r].out_degree() != 0)
  {
    // insert reference
    add_node_dna_to_sequence(ref, ref_nodes[r].get_label(), 0, 0xFFFFFFFFULL);

    // insert the variant node which contains the reference sequence
    add_node_dna_to_sequence(ref, var_nodes[v].get_label(), 0, 0xFFFFFFFFULL);

    v += ref_nodes[r].out_degree();
    ++r;
  }

  add_node_dna_to_sequence(ref, ref_nodes[r].get_label(), 0, 0xFFFFFFFFULL);
  return ref;
}


void
Graph::generate_reference_genome()
{
  reference = get_all_ref();
//  reference_offset = genomic_region.begin;
}


void
Graph::create_special_positions()
{
  ref_reach_to_special_pos.clear();
  ref_reach_poses.clear();

  for (unsigned r = 0; r < ref_nodes.size() - 1; ++r)
  {
    if (ref_nodes[r].out_degree() <= 1)
      continue;

    std::vector<TNodeIndex> const out_vars = ref_nodes[r].get_vars();
    assert(out_vars.size() >= 2);
    uint32_t const ref_label_reach = var_nodes[out_vars[0]].get_label().reach();
    uint32_t max_var_reach = var_nodes[out_vars[1]].get_label().reach();

    // Get the maximum variant reach
    for (auto it = out_vars.begin() + 2; it != out_vars.end(); ++it)
      max_var_reach = std::max(max_var_reach, var_nodes[*it].get_label().reach());

    // Create special position for each position further than the reference
    for (uint32_t reach = ref_label_reach + 1; reach <= max_var_reach; ++reach)
      add_special_pos(reach, ref_label_reach);
  }
}


std::vector<char>
Graph::get_generated_reference_genome(uint32_t & from, uint32_t & to) const
{
  // TODO: Handle multiregions
  // std::string const & chrom = genomic_regions[0].chr;
  uint32_t const abs_first_from = genomic_region.get_absolute_position(genomic_region.begin + 1);
  // uint32_t const abs_to = genomic_region.get_contig_position(to).second;
  from = std::max(abs_first_from, from);
  to = std::min(static_cast<uint32_t>(abs_first_from + reference.size()), to);

  if (to < from)
  {
    //BOOST_LOG_TRIVIAL(warning) << "[graphtyper::graph] Read end is before the reference genome "
    //                           << "starts.";
    //BOOST_LOG_TRIVIAL(debug) << "reference_offset="
    //                         << reference_offset
    //                         << ", to="
    //                         << to
    //                         << ", from="
    //                         << from
    //                         << ", abs_first_from="
    //                         << abs_first_from;

    return std::vector<char>(0);
  }

  return std::vector<char>(reference.begin() + from - abs_first_from,
                           reference.begin() + to - abs_first_from
                           );
}


std::vector<char>
Graph::get_ref(uint32_t from, uint32_t to) const
{
  return get_reference_ref(from, to);
}


std::vector<char>
Graph::get_reference_ref(uint32_t & from, uint32_t & to) const
{
  if (ref_nodes.size() == 0 or ref_nodes.front().get_label().order > to)
  {
    return std::vector<char>(0);
  }

  from = std::max(ref_nodes.front().get_label().order, from);
  to = std::min(ref_nodes.back().get_label().order + ref_nodes.back().get_label().dna.size(),
                static_cast<std::size_t>(to)
                );

  std::vector<char> ref(0);
  unsigned r = 0;
  unsigned v = 0;

  while (ref_nodes[r].out_degree() != 0)
  {
    // insert reference
    add_node_dna_to_sequence(ref, ref_nodes[r].get_label(), from, to);

    // insert the variant node which contains the reference sequence
    add_node_dna_to_sequence(ref, var_nodes[v].get_label(), from, to);

    v += ref_nodes[r].out_degree();
    ++r;
  }

  add_node_dna_to_sequence(ref, ref_nodes[r].get_label(), from, to);
  return ref;
}


std::vector<char>
Graph::walk_random_path(uint32_t from, uint32_t to) const
{
  if (ref_nodes.size() == 0 or ref_nodes.front().get_label().order > to)
  {
    return std::vector<char>(0);
  }

  std::vector<char> ref(0);
  unsigned r = 0;
  unsigned v = 0;

  while (ref_nodes[r].out_degree() != 0)
  {
    // insert reference
    add_node_dna_to_sequence(ref, ref_nodes[r].get_label(), from, to);

    // insert the variant node which contains the reference sequence
    uint32_t rand_v = v + rand() % ref_nodes[r].out_degree();
    add_node_dna_to_sequence(ref, var_nodes[rand_v].get_label(), from, to);

    v += ref_nodes[r].out_degree();
    ++r;
  }

  add_node_dna_to_sequence(ref, ref_nodes[r].get_label(), from, to);
  return ref;
}


std::vector<char>
Graph::get_first_var() const
{
  if (ref_nodes.size() == 0)
  {
    return std::vector<char>(0);
  }

  std::vector<char> sequence;
  unsigned r = 0;
  unsigned v = 0;

  while (ref_nodes[r].out_degree() != 0)
  {
    assert(ref_nodes[r].out_degree() != 1);

    // insert reference
    sequence.insert(sequence.end(), ref_nodes[r].get_label().dna.begin(), ref_nodes[r].get_label().dna.end());

    // insert the variant node which contains the reference
    sequence.insert(sequence.end(), var_nodes[v + 1].get_label().dna.begin(), var_nodes[v + 1].get_label().dna.end());

    v += ref_nodes[r].out_degree();
    ++r;
  }

  sequence.insert(sequence.end(), ref_nodes[r].get_label().dna.begin(), ref_nodes[r].get_label().dna.end());
  return sequence;
}


template <typename Archive>
void
Graph::serialize(Archive & ar, const unsigned int /*version*/)
{
  ar & ref_nodes;
  ar & var_nodes;
  ar & use_absolute_positions;
  ar & is_sv_graph;
  ar & genomic_region;
  ar & ref_reach_poses;
  ar & actual_poses;
  ar & ref_reach_to_special_pos;
  ar & SVs;
  ar & contigs;
}


std::size_t
Graph::size() const
{
  return ref_nodes.size() + var_nodes.size();
}


void
Graph::add_variants(VarRecord && record)
{
  // Create a nodes for the reference and the alternative variants
  {
    Label reference_label(record.pos, std::move(record.ref), 0);
    var_nodes.push_back(VarNode(std::move(reference_label), ref_nodes.size()));
  }

  for (long i = 0; i < static_cast<long>(record.alts.size()); ++i)
  {
    Label alternative_label(record.pos, std::move(record.alts[i]), i + 1);
    var_nodes.push_back(VarNode(std::move(alternative_label), ref_nodes.size()));
  }
}


bool
Graph::add_reference(unsigned end_pos,
                     unsigned const & num_var,
                     std::vector<char> const & reference_sequence)
{
  if (end_pos > reference_sequence.size() + genomic_region.begin)
  {
    end_pos = reference_sequence.size() + genomic_region.begin;
  }

  unsigned start_pos = genomic_region.begin;

  if (var_nodes.size() > 0)
  {
    assert(ref_nodes.size() > 0);
    auto const v = ref_nodes[ref_nodes.size() - 1].get_var_index(0);
    assert(v < var_nodes.size());
    Label const & previous_var_label = var_nodes[v].get_label();
    start_pos = previous_var_label.order + previous_var_label.dna.size();
  }

  // Make sure end_pos is larger or equal to start_pos
  bool is_continue = start_pos <= end_pos;
  assert(start_pos >= genomic_region.begin);

  auto start_it = (start_pos - genomic_region.begin) < reference_sequence.size() ?
                  reference_sequence.begin() + (start_pos - genomic_region.begin) :
                  reference_sequence.end();

  auto end_it = (end_pos - genomic_region.begin) < reference_sequence.size() ?
                reference_sequence.begin() + (end_pos - genomic_region.begin) :
                reference_sequence.end();

  std::vector<char> current_dna(start_it, end_it);

  // Create a vector of indexes
  std::vector<TNodeIndex> var_indexes(num_var);

  for (long i = 0; i < static_cast<long>(num_var); ++i)
    var_indexes[i] = i + var_nodes.size();

  ref_nodes.push_back(RefNode(Label(start_pos, std::move(current_dna), 0), std::move(var_indexes)));
  return is_continue;
}


void
Graph::break_apart_haplotypes(std::vector<Genotype> gts,
                              std::vector<Haplotype> & haplotypes,
                              int32_t max_read_length) const
{
  if (gts.size() == 0)
    return;

  if (max_read_length <= 1 || gts.size() == 1)
  {
    for (auto const & gt : gts)
    {
      Haplotype new_hap;

      if (gt.num > MAX_NUMBER_OF_HAPLOTYPES)
      {
        BOOST_LOG_TRIVIAL(warning) << __HERE__ << " There is a single genotype with " << gt.num
                                   << " variant sequences, more than the maximum number of "
                                   << "haplotypes.";
      }

      new_hap.add_genotype(Genotype(gt));
      haplotypes.push_back(std::move(new_hap));
    }

    return;
  }

  Haplotype new_hap;

  for (uint32_t i = 0; i < gts.size(); ++i)
  {
    new_hap.add_genotype(Genotype(gts[i]));
    assert(new_hap.gts.size() >= 1);
    uint32_t const v = gts[i].first_variant_node;
    uint32_t const r = var_nodes[v].get_out_ref_index();

    if (ref_nodes[r].get_label().dna.size() >= static_cast<uint32_t>(max_read_length) || i == gts.size() - 1)
    {
      // Check if we need to lower the max read length further
      if (new_hap.has_too_many_genotypes())
      {
        // keep breaking them apart
        break_apart_haplotypes(new_hap.gts, haplotypes, max_read_length - 1);
      }
      else
      {
        assert(new_hap.get_genotype_ids().size() > 0);
        haplotypes.push_back(std::move(new_hap));
      }

      new_hap.clear();
    }
  }
}


std::vector<Haplotype>
Graph::get_all_haplotypes(uint32_t variant_distance) const
{
  if (var_nodes.size() == 0)
    return std::vector<Haplotype>(0);

  std::vector<Haplotype> haplotypes;
  uint32_t v = 0;

  if (Options::instance()->is_one_genotype_per_haplotype)
  {
    for (uint32_t r = 0; r < ref_nodes.size() - 1; ++r)
    {
      assert(v < var_nodes.size());
      auto const & ref_node = ref_nodes[r];
      auto const & var_node = var_nodes[v];

      Haplotype hap;
      hap.add_genotype(Genotype(var_node.get_label().order, ref_node.out_degree(), uint32_t(v)));
      haplotypes.push_back(std::move(hap));

      v += ref_node.out_degree();
    }
  }
  else
  {
    Haplotype hap;

    // out_degree() == 1 means we have reached a sequence of Ns
    for (uint32_t r = 0; r < ref_nodes.size() - 1; ++r)
    {
      auto const & ref_node = ref_nodes[r];
      auto const & var_node = var_nodes[v];

      hap.add_genotype(Genotype(var_node.get_label().order, ref_node.out_degree(), uint32_t(v)));

      assert(v + ref_node.out_degree() >= var_nodes.size() ||
             var_node.get_label().order != var_nodes[v + ref_node.out_degree()].get_label().order
             );

      assert(var_node.get_out_ref_index() == r + 1);

      // Check if we should move on to the next haplotype, instead of adding more to this one
      if (ref_nodes[r + 1].get_label().dna.size() >= variant_distance)
      {
        if (hap.has_too_many_genotypes())
          break_apart_haplotypes(hap.gts, haplotypes, variant_distance - 1);
        else
          haplotypes.push_back(std::move(hap));

        hap.clear(); // This haplotype has been added, clear this one and move on to the next one
      }

      v += ref_node.out_degree();
    }

    // Add the last haplotype, if it is not empty
    if (!hap.gts.empty())
    {
      if (hap.has_too_many_genotypes())
        break_apart_haplotypes(hap.gts, haplotypes, variant_distance - 1);
      else
        haplotypes.push_back(std::move(hap));
    }
  }

// #ifndef NDEBUG
//   for (auto & hap : haplotypes)
//     hap.check_for_duplicate_haplotypes();
// #endif // NDEBUG

  return haplotypes;
}


std::vector<char>
Graph::get_sequence_of_a_haplotype_call(std::vector<Genotype> const & gts,
                                        uint32_t const haplotype_call) const
{
  assert(gts.size() > 0);
  uint32_t rem = haplotype_call;
  uint32_t num = 1;

  for (auto gt_it = gts.cbegin(); gt_it != gts.cend(); ++gt_it)
    num *= gt_it->num;

  std::vector<uint32_t> var_calls;
  var_calls.reserve(gts.size());

  for (unsigned i = 0; i < gts.size(); ++i)
  {
    num /= gts[i].num;
    var_calls.push_back(rem / num);
    rem %= num;
  }

  assert(var_calls.size() > 0);
  uint32_t const first_v = gts[0].first_variant_node + var_calls[0];
  uint32_t const first_r = var_nodes[first_v].get_out_ref_index() - 1;
  assert(first_r < ref_nodes.size());
  std::vector<char> sequence(0);

  // Get previous base
  if (ref_nodes[first_r].get_label().dna.size() > 0)
  {
    // If the reference before the variant is non-empty, we can get the common base from there
    sequence.push_back(ref_nodes[first_r].get_label().dna[ref_nodes[first_r].get_label().dna.size() - 1]);
  }
  else if (gts[0].first_variant_node != 0)
  {
    // Otherwise we try to get it from the previous variant (but we can only do that if there is such variant)
    uint32_t const prev_r = var_nodes[gts[0].first_variant_node - 1ul].get_out_ref_index() - 1;
    assert(prev_r < ref_nodes.size());
    sequence.push_back(var_nodes[ref_nodes[prev_r].get_var_index(0)].get_label().dna.back());
  }
  else
  {
    // If both previous methods fail, we use 'N'
    sequence.push_back('N');
  }

  assert(sequence.size() == 1);
  sequence.insert(sequence.end(), var_nodes[first_v].get_label().dna.begin(), var_nodes[first_v].get_label().dna.end());

  for (unsigned i = 1; i < var_calls.size(); ++i)
  {
    uint32_t const v = gts[i].first_variant_node + var_calls[i];
    uint32_t const r = var_nodes[v].get_out_ref_index() - 1;

    // Add reference
    sequence.insert(sequence.end(), ref_nodes[r].get_label().dna.begin(), ref_nodes[r].get_label().dna.end());

    // Add variant
    sequence.insert(sequence.end(), var_nodes[v].get_label().dna.begin(), var_nodes[v].get_label().dna.end());
  }

  return sequence;
}


std::vector<std::vector<char> >
Graph::get_all_sequences_of_a_genotype(Genotype const & gt) const
{
  assert(gt.first_variant_node < var_nodes.size());
  uint32_t const r = var_nodes[gt.first_variant_node].get_out_ref_index() - 1;
  std::vector<char> first_base;

  if (ref_nodes[r].get_label().dna.size() > 0)
    first_base.push_back(ref_nodes[r].get_label().dna.back());
  else if (r != 0)
    first_base.push_back(var_nodes[ref_nodes[r - 1].get_vars()[0]].get_label().dna.back());
  else
    first_base.push_back('N');

  std::vector<TNodeIndex> const & vars = ref_nodes[r].get_vars();
  std::vector<std::vector<char> > seqs(vars.size(), first_base);

  for (unsigned i = 0; i < vars.size(); ++i)
    seqs[i].insert(seqs[i].end(), var_nodes[vars[i]].get_label().dna.begin(), var_nodes[vars[i]].get_label().dna.end());

  return seqs;
}


bool
Graph::is_variant_in_graph(Variant const & var) const
{
  if (var_nodes.size() == 0)
    return false;

  uint32_t r = 1;

  while (r < ref_nodes.size())
  {
    if (ref_nodes[r].get_label().order >= var.abs_pos)
      break;

    ++r;
  }

  --r; // Move one back

  if (ref_nodes[r].out_degree() <= 1)
    return false;

  uint32_t v = ref_nodes[r].get_var_index(0);
  Variant new_var(Genotype(var_nodes[v].get_label().order, ref_nodes[r].out_degree(), v));

  // Test if the variant is the same
  {
    Variant new_var2(new_var);
    new_var2.normalize();

    if (new_var2 == var)
      return true;
  }

  //// Try to break down the variant and see if it is the SamReader
  //std::size_t const THRESHOLD = 1;
  //std::vector<Variant> broken_vars = break_down_variant(std::move(new_var), THRESHOLD);
  //
  //for (auto & broken_var : broken_vars)
  //{
  //  broken_var.normalize();
  //
  //  if (broken_var == var)
  //    return true;
  //}

  return false;
}


uint8_t
Graph::get_10log10_num_paths(TNodeIndex const v, uint32_t const MAX_DISTANCE)
{
  auto to_log10 = [](uint32_t degree){
                    return 10.0 * log10(static_cast<double>(degree));
                  };

  assert(v < var_nodes.size());
  assert(var_nodes[v].get_label().reach() >= var_nodes[v].get_label().order);

  // Get the center of the variant, this position will be used to determine if the other variants are too far away or not.
  uint32_t const CENTER = var_nodes[v].get_label().order +
                          (var_nodes[v].get_label().reach() - var_nodes[v].get_label().order) / 2;
  uint32_t r = var_nodes[v].get_out_ref_index() - 1;
  double num_paths = 0.0;

  // Go forward
  while (r < ref_nodes.size())
  {
    if (ref_nodes[r].get_label().reach() > (CENTER + MAX_DISTANCE) || ref_nodes[r].out_degree() == 0)
      break;

    num_paths += to_log10(ref_nodes[r].out_degree());
    ++r;
  }

  r = var_nodes[v].get_out_ref_index() - 1; // Reset r

  // Go backwards
  while (r > 0)
  {
    if (ref_nodes[r - 1].get_label().reach() < (CENTER - MAX_DISTANCE))
      break;

    num_paths += to_log10(ref_nodes[r - 1].out_degree());
    --r;
  }

  if (num_paths >= 254.5)
    return 255u;
  else
    return static_cast<uint8_t>(std::lround(num_paths));
}


std::vector<Location>
Graph::get_locations_of_an_actual_position(uint32_t pos, Path const & path, bool const is_special) const
{
  assert(ref_nodes.size() != 0);
  std::vector<Location> locs(0);

  if ((pos < ref_nodes[0].get_label().order) ||
      (pos > (ref_nodes.back().get_label().reach()))
      )
  {
    return locs;
  }

  for (uint32_t r = 1; r <= ref_nodes.size(); ++r)
  {
    if (r < ref_nodes.size() && ref_nodes[r].get_label().order <= pos)
      continue;

    int rr = r - 1;

    if (pos < (ref_nodes[rr].get_label().order + ref_nodes[rr].get_label().dna.size()))
    {
      // Ref covers this location
      if (!is_special)
      {
        locs.push_back(Location('R' /*type*/,
                                static_cast<uint32_t>(rr) /*node_id*/,
                                ref_nodes[rr].get_label().order /*node_order*/,
                                pos - ref_nodes[rr].get_label().order /*offset*/
                                )
                       );
        break; // There is no way there are also variants at this location if the position is not special
      }

      assert(rr > 0);
      --rr; // Variants behind the reference can only have this location
    }

    long const PADDING = is_sv_graph ? 1000000 : 1000;

    // Check variants behind this reference
    while (rr >= 0 && ref_nodes[rr].get_label().reach() + PADDING > pos)
    {
      // Assume there is no variants larger than 5000 bp
      for (int i = 0; i < static_cast<int>(ref_nodes[rr].out_degree()); ++i)
      {
        uint32_t const v = static_cast<uint32_t>(ref_nodes[rr].get_var_index(i));

        if (pos >= var_nodes[v].get_label().order and pos <= var_nodes[v].get_label().reach())
        {
          // Only add this node if the path has it
          auto find_it = std::find(path.var_order.cbegin(),
                                   path.var_order.cend(),
                                   var_nodes[v].get_label().order
                                   );

          long const j = std::distance(path.var_order.cbegin(), find_it);
          assert(j >= 0);
          assert(i == this->get_variant_num(v));

          if (path.is_empty() || (j < static_cast<long>(path.nums.size()) && path.nums[j].test(i)))
          {
            locs.push_back(
              {'V' /*type*/,
               v /*node_id*/,
               var_nodes[v].get_label().order /*node_order*/,
               pos - var_nodes[v].get_label().order  /*offset*/
              }
              );
          }
        }
      }

      --rr;
    }

    break;
  }

  return locs;
}


std::vector<std::vector<char> >
Graph::get_sequence_from_location(Location const & loc,
                                  uint32_t const length,
                                  std::vector<char> const & prefix
                                  ) const
{
  std::vector<std::vector<char> > seqs(1, prefix);
  long r = ref_nodes.size();

  if (loc.node_type == 'V')
  {
    auto const & l = var_nodes[loc.node_index].get_label();
    assert(loc.offset < l.dna.size());
    std::copy(l.dna.begin() + loc.offset, l.dna.end(), std::back_inserter(seqs[0]));

    if (seqs[0].size() < length)
    {
      r = var_nodes[loc.node_index].get_out_ref_index();
      long const dist = std::min(static_cast<long>(ref_nodes[r].get_label().dna.size()),
                                 static_cast<long>(length - seqs[0].size())
                                 );

      if (dist > 0)
      {
        std::copy(ref_nodes[r].get_label().dna.begin(),
                  ref_nodes[r].get_label().dna.begin() + dist,
                  std::back_inserter(seqs[0])
                  );
      }
    }
  }
  else
  {
    assert(loc.node_type == 'R');
    auto const & l = ref_nodes[loc.node_index].get_label();
    assert(loc.offset < l.dna.size());
    std::copy(l.dna.begin() + loc.offset, l.dna.end(), std::back_inserter(seqs[0]));
    r = loc.node_index;
  }

  while (r < static_cast<long>(ref_nodes.size() - 1))
  {
    // Add reference sequence on alternative node
    long v = ref_nodes[r].get_var_index(0);

    // Label of the reference variant node
    auto const & l = var_nodes[v].get_label();

    // Add alternative sequence alternative node
    for (v = 1; v < static_cast<int>(ref_nodes[r].out_degree()); ++v)
    {
      auto const & lv = var_nodes[v].get_label();

      // Skip if the same size as the reference allele
      if (lv.dna.size() == l.dna.size())
        continue;

      seqs.push_back(seqs[0]);
      std::copy(lv.dna.begin(),
                lv.dna.end(),
                std::back_inserter(seqs[seqs.size() - 1])
                );
    }

    // Extend the reference sequence with the reference variant node
    std::copy(l.dna.begin(), l.dna.end(), std::back_inserter(seqs[0]));

    ++r;
    auto const & lr = ref_nodes[r].get_label();
    bool is_any_remaining = false;

    for (auto & seq : seqs)
    {
      long const dist = std::min(static_cast<long>(lr.dna.size()),
                                 static_cast<long>(length) - static_cast<long>(seq.size())
                                 );

      if (dist > 0)
      {
        std::copy(lr.dna.begin(), lr.dna.begin() + dist, std::back_inserter(seq));

        if (seq.size() < length)
          is_any_remaining = true;
      }
    }

    if (!is_any_remaining)
      break;
  }

  return seqs;
}


std::unordered_set<long>
Graph::reference_distance_between_locations(std::vector<Location> const & ll1,
                                            std::vector<Location> const & ll2
                                            ) const
{
  std::unordered_set<long> distance_map;

  auto get_ref_pos_lambda = [&](Location const & l)
                            {
                              uint32_t pos;

                              if (l.node_type == 'R')
                              {
                                pos = l.node_order + l.offset;
                              }
                              else
                              {
                                assert(l.node_type == 'V');
                                long const node_remainder = var_nodes[l.node_index].get_label().dna.size() - l.offset;
                                assert(node_remainder >= 0);
                                pos = this->var_nodes[l.node_index].get_label().reach() + 1 - node_remainder;
                              }

                              return pos;
                            };

  for (auto const & l1 : ll1)
  {
    uint32_t const ref_pos1 = get_ref_pos_lambda(l1);

    for (auto const & l2 : ll2)
    {
      uint32_t const ref_pos2 = get_ref_pos_lambda(l2);
      distance_map.insert(static_cast<int64_t>(ref_pos2) - static_cast<int64_t>(ref_pos1));
    }
  }

  return distance_map;
}


std::vector<Location>
Graph::get_locations_of_a_position(uint32_t pos, Path const & path) const
{
  bool const IS_SPECIAL = is_special_pos(pos);

  if (IS_SPECIAL)
    pos = actual_poses.at(pos - SPECIAL_START);

  std::vector<Location> locs = this->get_locations_of_an_actual_position(pos, path, IS_SPECIAL);

#ifndef NDEBUG
  if (locs.size() == 0)
  {
    BOOST_LOG_TRIVIAL(warning) << "[graphtyper::graph] Found no location for position " << pos << " (start "
                               << ref_nodes[0].get_label().order
                               << ", end "
                               << ref_nodes.back().get_label().reach() << ")";
  }
#endif

  return locs;
}


std::vector<KmerLabel>
Graph::get_labels_forward(Location const & s,
                          std::vector<char> const & read,
                          uint32_t & max_mismatches
                          ) const
{
  std::vector<KmerLabel> labels;

  std::vector<std::vector<char> > var_and_refs(1);
  std::vector<std::vector<uint32_t> > var_ids(1);
  std::vector<uint32_t> end_pos(1, 0u);
  std::vector<TNodeIndex> vars;

  if (s.node_type == 'V')
  {
    assert(s.node_index < var_nodes.size());
    VarNode const & var = var_nodes[s.node_index];
    var_ids[0].push_back(s.node_index);
    var_and_refs[0] = std::vector<char>(var.get_label().dna.begin() + s.offset,
                                        var.get_label().dna.end()
                                        );

    // Check if variant is enough
    if (var_and_refs[0].size() >= read.size())
    {
      // variant is enough
      end_pos[0] =
        static_cast<uint32_t>(var.get_label().reach() - (var_and_refs[0].size() - read.size()));

      uint32_t const ref_reach =
        var_nodes[ref_nodes[var.get_out_ref_index() - 1].get_vars()[0]].get_label().reach();

      if (end_pos[0] > ref_reach)
        end_pos[0] = get_special_pos(end_pos[0], ref_reach);
    }
    else
    {
      // We also need to add a reference
      RefNode const & ref = ref_nodes[var.get_out_ref_index()];
      vars = ref.get_vars();

      var_and_refs[0].insert(var_and_refs[0].end(),
                             ref.get_label().dna.begin(),
                             ref.get_label().dna.end()
                             );

      end_pos[0] =
        static_cast<uint32_t>(ref.get_label().reach() - (var_and_refs[0].size() - read.size()));
    }
  }
  else
  {
    assert(s.node_type == 'R');
    assert(s.node_index < ref_nodes.size());
    RefNode const & ref = ref_nodes[s.node_index];
    vars = ref.get_vars();

    var_and_refs[0] = std::vector<char>(ref.get_label().dna.begin() + s.offset,
                                        ref.get_label().dna.end()
                                        );

    end_pos[0] =
      static_cast<uint32_t>(ref.get_label().reach() - (var_and_refs[0].size() - read.size()));
  }

  // We are starting on a variant node
  if (vars.size() > 0 && var_and_refs[0].size() < read.size())
  {
    // We are the the end of the graph, and the sequence is not long enough, we need to bail
    std::vector<uint32_t> mismatch_scores(vars.size(), 0);
    uint32_t r = var_nodes[vars[0]].get_out_ref_index();
    bool all_sequences_long_enough = false;
    std::size_t const MAX_VAR_AND_REFS = 128;

    while (not all_sequences_long_enough && var_and_refs.size() < MAX_VAR_AND_REFS && vars.size() > 0)
    {
      all_sequences_long_enough = true;
      assert(r < ref_nodes.size());
      RefNode const & ref = ref_nodes[r];
      std::size_t original_size = var_and_refs.size();

      for (unsigned j = 0; j < original_size; ++j)
      {
        assert(j < var_and_refs.size());    // Should always be less than the current size

        if (var_and_refs[j].size() >= read.size())
          continue;   // Sequence is already large enough

        for (unsigned i = 0; i < vars.size() - 1; ++i)
        {
          assert(j < var_and_refs.size());
          assert(vars[i] < var_nodes.size());
          VarNode const & var = var_nodes[vars[i]];
          std::vector<char> new_seq(var_and_refs[j].begin(), var_and_refs[j].end());
          new_seq.insert(new_seq.end(), var.get_label().dna.begin(), var.get_label().dna.end());

          bool const variant_is_enough = new_seq.size() >= read.size();

          if (not variant_is_enough)
            new_seq.insert(new_seq.end(), ref.get_label().dna.begin(), ref.get_label().dna.end());

          // Only add it if it has less or equal than 'max_mismatches' mismatches
          if (count_mismatches(read, 0, new_seq, 0, max_mismatches) <= max_mismatches)
          {
            std::vector<uint32_t> new_var_id(var_ids[j]);
            new_var_id.push_back(vars[i]);
            var_ids.push_back(std::move(new_var_id));

            // Check if we need to continue further
            if (new_seq.size() < read.size())
              all_sequences_long_enough = false;

            // Update end positions
            if (variant_is_enough)
            {
              end_pos.push_back(var.get_label().reach() - (new_seq.size() - read.size()));

              // Check if the end position is further than the reference reach
              uint32_t const ref_reach =
                var_nodes[ref_nodes[var.get_out_ref_index() - 1].get_vars()[0]].get_label().reach();

              if (end_pos.back() > ref_reach)
                end_pos.back() = get_special_pos(end_pos.back(), ref_reach);
            }
            else
            {
              end_pos.push_back(ref.get_label().reach() - (new_seq.size() - read.size()));
            }

            assert(var_nodes[var_ids.back().back()].get_label().order <= end_pos.back());
            var_and_refs.push_back(std::move(new_seq));
          }
        }

        // The last variant replaces the old seq
        VarNode const & var = var_nodes[vars[vars.size() - 1]];
        var_and_refs[j].insert(var_and_refs[j].end(), var.get_label().dna.begin(), var.get_label().dna.end());

        bool const variant_is_enough = var_and_refs[j].size() >= read.size();

        if (!variant_is_enough)
        {
          var_and_refs[j].insert(var_and_refs[j].end(), ref.get_label().dna.begin(), ref.get_label().dna.end());
        }

        if (count_mismatches(read, 0, var_and_refs[j], 0, max_mismatches) <= max_mismatches)
        {
          var_ids[j].push_back(vars[vars.size() - 1]);

          if (all_sequences_long_enough and var_and_refs[j].size() < read.size())
            all_sequences_long_enough = false;

          // Update end positions
          if (variant_is_enough)
          {
            end_pos[j] = var.get_label().reach() - (var_and_refs[j].size() - read.size());

            // Check if the end position is further than the reference reach
            uint32_t const ref_reach =
              var_nodes[ref_nodes[var.get_out_ref_index() - 1].get_vars()[0]].get_label().reach();
            if (end_pos[j] > ref_reach)
              end_pos[j] = get_special_pos(end_pos[j], ref_reach);
          }
          else
          {
            end_pos[j] = ref.get_label().reach() - (var_and_refs[j].size() - read.size());
          }

          assert(var_nodes[var_ids[j].back()].get_label().order <= end_pos[j]);
          assert(var_ids.size() == end_pos.size());
        }
        else
        {
          // Delete the jth element
          var_and_refs.erase(var_and_refs.begin() + j);
          var_ids.erase(var_ids.begin() + j);
          end_pos.erase(end_pos.begin() + j);

          --original_size;
          --j;
        }
      }

      if (not all_sequences_long_enough)
      {
        // Get new reference node and variant nodes
        assert(r < ref_nodes.size());
        vars = ref_nodes[r].get_vars();
        ++r;
      }
      else
      {
        break;
      }
    }
  }

  std::vector<std::vector<uint32_t> > best_var_ids;
  std::vector<uint32_t> best_end_pos;

  // Iterate all possible sequences
  for (unsigned j = 0; j < var_and_refs.size(); ++j)
  {
    if (var_and_refs[j].size() < read.size())
      continue;

    uint32_t mismatches = count_mismatches(read, 0, var_and_refs[j], 0, max_mismatches);

    if (mismatches > max_mismatches)
    {
      continue;
    }
    else if (mismatches < max_mismatches)
    {
      max_mismatches = mismatches; // Found alignment with fewer mismatches
      best_var_ids.clear();
      best_var_ids.push_back(var_ids[j]);
      best_end_pos.clear();
      best_end_pos.push_back(end_pos[j]);
    }
    else
    {
      best_var_ids.push_back(var_ids[j]);
      best_end_pos.push_back(end_pos[j]);
    }
  }

  if (best_var_ids.size() == 0)
    return labels;

  assert(best_var_ids.size() == best_end_pos.size());

  if (best_var_ids.size() > 0)
  {
    for (unsigned j = 0; j < best_var_ids.size(); ++j)
    {
      uint32_t start_pos = s.node_order + s.offset;

      // Check if we need to use a special positions for the end position
      if (s.node_type == 'V')
      {
        uint32_t const ref_reach =
          var_nodes[ref_nodes[var_nodes[s.node_index].get_out_ref_index() - 1].get_vars()[0]].get_label().reach();
        if (start_pos > ref_reach)
          start_pos = get_special_pos(start_pos, ref_reach);
      }

      // Check if we are overlapping any variant node
      if (best_var_ids[j].size() == 0)
      {
        labels.push_back(KmerLabel(start_pos, best_end_pos[j]));
      }
      else
      {
        for (auto const & good_var : best_var_ids[j])
        {
          assert(var_nodes[good_var].get_label().order <= best_end_pos[j]);

          labels.push_back(KmerLabel(start_pos,
                                     best_end_pos[j],
                                     good_var)
                           );
        }
      }
    }
  }

  return labels;
}


std::vector<KmerLabel>
Graph::get_labels_backward(Location const & e,
                           std::vector<char> const & read,
                           uint32_t & max_mismatches
                           ) const
{
  std::vector<KmerLabel> labels;

  std::vector<std::vector<char> > var_and_refs(1);
  std::vector<std::vector<uint32_t> > var_ids(1);
  std::vector<uint32_t> start_pos(1, 0u);
  std::vector<TNodeIndex> vars;

  if (e.node_type == 'V')
  {
    assert(e.node_index < var_nodes.size());
    VarNode const & var = var_nodes[e.node_index];
    var_ids[0].push_back(e.node_index);
    var_and_refs[0] = std::vector<char>(var.get_label().dna.begin(), var.get_label().dna.begin() + e.offset + 1);

    // Check if adding the variant was enough
    if (var_and_refs[0].size() >= read.size())
    {
      start_pos[0] = var.get_label().order + (var_and_refs[0].size() - read.size());

      // Check if we need to use a special positions
      uint32_t const ref_reach = var_nodes[ref_nodes[var.get_out_ref_index() - 1].get_vars()[0]].get_label().reach();
      if (start_pos[0] > ref_reach)
        start_pos[0] = get_special_pos(start_pos[0], ref_reach);
    }
    else
    {
      uint32_t const r = var.get_out_ref_index() - 1;
      RefNode const & ref = ref_nodes[r];
      var_and_refs[0].insert(var_and_refs[0].begin(), ref.get_label().dna.begin(), ref.get_label().dna.end());
      start_pos[0] = ref.get_label().order + (var_and_refs[0].size() - read.size());

      if (r != 0)
        vars = ref_nodes[r - 1].get_vars();
    }
  }
  else
  {
    assert(e.node_type == 'R');
    assert(e.node_index < ref_nodes.size());
    RefNode const & ref = ref_nodes[e.node_index];

    if (e.node_index != 0)
      vars = ref_nodes[e.node_index - 1].get_vars(); // Only if we are not on the first reference node, we can get the vars

    var_and_refs[0] = std::vector<char>(ref.get_label().dna.begin(), ref.get_label().dna.begin() + e.offset + 1);
    start_pos[0] = ref.get_label().order + (var_and_refs[0].size() - read.size());
  }

  // We are starting on a variant node
  if (vars.size() > 0 && var_and_refs[0].size() < read.size())
  {
    std::vector<uint32_t> mismatch_scores(vars.size(), 0);
    uint32_t r = var_nodes[vars[0]].get_out_ref_index() - 1;
    bool all_sequences_long_enough = false;
    std::size_t const MAX_VAR_AND_REFS = 128;

    while (not all_sequences_long_enough and var_and_refs.size() < MAX_VAR_AND_REFS && vars.size() > 0)
    {
      all_sequences_long_enough = true;
      assert(r < ref_nodes.size());
      RefNode const & ref = ref_nodes[r];
      std::size_t original_size = var_and_refs.size();

      for (unsigned j = 0; j < original_size; ++j)
      {
        assert(j < var_and_refs.size());  // Should always be less than the current size

        if (var_and_refs[j].size() >= read.size())
          continue; // Sequence is already large enough

        for (unsigned i = 0; i < vars.size() - 1; ++i)
        {
          assert(j < var_and_refs.size());

          if (var_and_refs[j].size() < read.size())
          {
            assert(i < vars.size());
            assert(vars[i] < var_nodes.size());
            VarNode const & var = var_nodes[vars[i]];
            std::vector<char> new_seq(var.get_label().dna.begin(), var.get_label().dna.end());
            new_seq.insert(new_seq.end(), var_and_refs[j].begin(), var_and_refs[j].end());

            bool const variant_is_enough = new_seq.size() >= read.size();

            if (not variant_is_enough)
            {
              new_seq.insert(new_seq.begin(), ref.get_label().dna.begin(), ref.get_label().dna.end());
            }

            // Only add it if it has less or equal than 'max_mismatches' mismatches
            if (count_mismatches_backward(read, 0, new_seq, 0, max_mismatches) <= max_mismatches)
            {
              std::vector<uint32_t> new_var_id(var_ids[j]);
              new_var_id.push_back(vars[i]);
              var_ids.push_back(std::move(new_var_id));

              // Check if we need to continue further
              if (new_seq.size() < read.size())
              {
                all_sequences_long_enough = false;
              }

              // Update end positions
              if (variant_is_enough)
              {
                start_pos.push_back(var.get_label().order + (new_seq.size() - read.size()));

                // Check if we need to use a special positions
                uint32_t const ref_reach =
                  var_nodes[ref_nodes[var.get_out_ref_index() - 1].get_vars()[0]].get_label().reach();
                if (start_pos.back() > ref_reach)
                  start_pos.back() = get_special_pos(start_pos.back(), ref_reach);
              }
              else
              {
                start_pos.push_back(ref.get_label().order + (new_seq.size() - read.size()));
              }

              var_and_refs.push_back(std::move(new_seq));
            }
          }
        }

        // The last variant replaces the old seq
        VarNode const & var = var_nodes[vars[vars.size() - 1]];
        var_and_refs[j].insert(var_and_refs[j].begin(), var.get_label().dna.begin(), var.get_label().dna.end());

        bool const variant_is_enough = var_and_refs[j].size() >= read.size();

        if (!variant_is_enough)
        {
          var_and_refs[j].insert(var_and_refs[j].begin(), ref.get_label().dna.begin(), ref.get_label().dna.end());
        }

        if (count_mismatches_backward(read, 0, var_and_refs[j], 0, max_mismatches) <= max_mismatches)
        {
          var_ids[j].push_back(vars[vars.size() - 1]);

          if (var_and_refs[j].size() < read.size())
            all_sequences_long_enough = false;

          // Update end positions
          if (variant_is_enough)
          {
            start_pos[j] = var.get_label().order + (var_and_refs[j].size() - read.size());

            // Check if we need to use a special positions
            uint32_t const ref_reach =
              var_nodes[ref_nodes[var.get_out_ref_index() - 1].get_vars()[0]].get_label().reach();
            if (start_pos[j] > ref_reach)
              start_pos[j] = get_special_pos(start_pos[j], ref_reach);
          }
          else
          {
            start_pos[j] = ref.get_label().order + (var_and_refs[j].size() - read.size());
          }

          assert(var_ids.size() == start_pos.size());
        }
        else
        {
          // Delete the jth element
          var_and_refs.erase(var_and_refs.begin() + j);
          var_ids.erase(var_ids.begin() + j);
          start_pos.erase(start_pos.begin() + j);

          --original_size;
          --j;
        }
      }

      if (not all_sequences_long_enough)
      {
        // // Get new reference node and variant nodes
        // r = var_nodes[vars[0]].get_out_ref_index() - 1;
        // assert(r < ref_nodes.size());

        if (r != 0)
        {
          --r;
          assert(ref_nodes[r].get_vars()[0] != vars[0]);
          vars = ref_nodes[r].get_vars();
        }
        else
        {
          vars.clear();
          break;
        }
      }
      else
      {
        break;
      }
    }
  }

  std::vector<std::vector<uint32_t> > best_var_ids;
  std::vector<uint32_t> best_start_pos;

  // Iterate all possible sequences
  for (unsigned j = 0; j < var_and_refs.size(); ++j)
  {
    if (var_and_refs[j].size() < read.size())
      continue;

    uint32_t const mismatches = count_mismatches_backward(read, 0, var_and_refs[j], 0, max_mismatches);

    if (mismatches < max_mismatches)
    {
      max_mismatches = mismatches;
      best_var_ids.clear();
      best_var_ids.push_back(var_ids[j]);
      best_start_pos.clear();
      best_start_pos.push_back(start_pos[j]);
    }
    else if (mismatches == max_mismatches)
    {
      best_var_ids.push_back(var_ids[j]);
      best_start_pos.push_back(start_pos[j]);
    }
  }

  if (best_var_ids.size() == 0)
    return labels;

  assert(best_var_ids.size() == best_start_pos.size());

  for (unsigned j = 0; j < best_var_ids.size(); ++j)
  {
    uint32_t end_pos = e.node_order + e.offset;

    // Check if we need to use a special positions for the end position
    if (e.node_type == 'V')
    {
      uint32_t const ref_reach =
        var_nodes[ref_nodes[var_nodes[e.node_index].get_out_ref_index() - 1].get_vars()[0]].get_label().reach();
      if (end_pos > ref_reach)
        end_pos = get_special_pos(end_pos, ref_reach);
    }

    // Check if we are overlapping any variant node
    if (best_var_ids[j].size() == 0)
    {
      labels.push_back(KmerLabel(best_start_pos[j], end_pos));
    }
    else
    {
      for (auto const & good_var : best_var_ids[j])
      {
        labels.push_back(KmerLabel(best_start_pos[j],
                                   end_pos,
                                   good_var)
                         );
      }
    }
  }

  return labels;
}


std::vector<KmerLabel>
Graph::iterative_dfs(std::vector<Location> const & start_locations,
                     std::vector<Location> const & end_locations,
                     std::vector<char> const & subread,
                     uint32_t & max_mismatches
                     ) const
{
  std::vector<KmerLabel> labels;
  assert(start_locations.size() > 0);
  assert(end_locations.size() > 0);
  assert(subread.size() > 0);
  std::size_t const MAX_LOCATIONS = 1024;

  if (start_locations.size() > MAX_LOCATIONS || end_locations.size() > MAX_LOCATIONS)
    return labels;

  auto add_if_better =
    [&labels, &max_mismatches](std::vector<KmerLabel> && new_labels, uint32_t const mismatches)
    {
      if (new_labels.size() > 0)
      {
        if (mismatches < max_mismatches)
        {
          max_mismatches = mismatches;
          labels = std::move(new_labels);
        }
        else if (mismatches == max_mismatches)
        {
          std::move(new_labels.begin(), new_labels.end(), std::back_inserter(labels));
        }
      }
    };

  // Check if node type of start location is unavailable ('U'). In this case we need to walk the graph backwards
  if (start_locations.size() == 1 and start_locations[0].is_unavailable())
  {
    for (auto const & e : end_locations)
    {
      uint32_t mismatches = max_mismatches;
      std::vector<KmerLabel> new_labels = get_labels_backward(e, subread, mismatches);
      add_if_better(std::move(new_labels), mismatches);
    }
  }
  else
  {
    for (auto const & s : start_locations)
    {
      uint32_t mismatches = max_mismatches;
      std::vector<KmerLabel> new_labels = get_labels_forward(s, subread, mismatches);
      add_if_better(std::move(new_labels), mismatches);
    }
  }

  return labels;
}


/**
 * SPECIAL POS
 */
void
Graph::add_special_pos(uint32_t const actual_pos, uint32_t const ref_reach)
{
  ref_reach_poses.push_back(ref_reach);
  actual_poses.push_back(actual_pos);
  auto find_it = ref_reach_to_special_pos.find(ref_reach);

  if (find_it != ref_reach_to_special_pos.end())
  {
    find_it->second.push_back(SPECIAL_START + this->ref_reach_poses.size() - 1);
  }
  else
  {
    ref_reach_to_special_pos[ref_reach] = std::vector<uint32_t>(1, SPECIAL_START + ref_reach_poses.size() - 1);
  }
}


uint32_t
Graph::get_special_pos(uint32_t const pos, uint32_t const ref_reach) const
{
  assert(pos > ref_reach);
  assert(std::distance(ref_reach_to_special_pos.begin(), ref_reach_to_special_pos.end()) > 0);
  assert(ref_reach_to_special_pos.count(ref_reach) == 1);
  assert(pos - ref_reach - 1 < ref_reach_to_special_pos.at(ref_reach).size());
  return ref_reach_to_special_pos.at(ref_reach).at(pos - ref_reach - 1);
}


bool
Graph::is_special_pos(uint32_t const pos) const
{
  return pos >= SPECIAL_START && (pos - SPECIAL_START) < ref_reach_poses.size();
}


uint32_t
Graph::get_ref_reach_pos(uint32_t const pos) const
{
  if (is_special_pos(pos))
    return ref_reach_poses.at(pos - SPECIAL_START);
  else
    return pos;
}


uint32_t
Graph::get_actual_pos(uint32_t const pos) const
{
  if (is_special_pos(pos))
    return actual_poses.at(pos - SPECIAL_START);
  else
    return pos;
}


/**
 * ERROR CHECKING
 */

bool
Graph::check() const
{
  return check_ACGTN_only()
         && check_empty_variant_dna()
         && check_increasing_order();
  //&& check_if_order_follows_reference();
}


bool
Graph::check_ACGTN_only() const
{
  bool no_non_ACGTN_errors = true;

  // Reference nodes
  for (std::size_t r = 0; r < ref_nodes.size(); ++r)
  {
    Label const & label = ref_nodes[r].get_label();

    for (long i = 0; i < static_cast<long>(label.dna.size()); ++i)
    {
      char const c = label.dna[i];

      switch (c)
      {
      case 'A': continue;

      case 'C': continue;

      case 'G': continue;

      case 'T': continue;

      case 'N': continue;

      default:
      {
        BOOST_LOG_TRIVIAL(warning) << "[" << __HERE__ << "] Reference node " << r
                                   << " has a " << c << " (int=" << ((int)c) << ") at "
                                   << genomic_region.get_contig_position(label.order, *this).first << ":"
                                   << genomic_region.get_contig_position(label.order, *this).second << " + " << i;
        no_non_ACGTN_errors = false;
      }
      }
    }
  }

  // Variant nodes
  for (std::size_t v = 0; v < var_nodes.size(); ++v)
  {
    Label const & label = var_nodes[v].get_label();

    for (long i = 0; i < static_cast<long>(label.dna.size()); ++i)
    {
      char const c = label.dna[i];

      switch (c)
      {
      case 'A':
      case 'C':
      case 'G':
      case 'T':
      case 'N': continue;

      case '<': // Ignore tag
      {
        ++i;
        while (i < static_cast<long>(label.dna.size()) && (label.dna[i] != '>' || label.dna[i] != '<'))
          ++i;

        --i;
        continue;
      }

      default:
      {
        BOOST_LOG_TRIVIAL(warning) << "[" << __HERE__ << "] Variant node " << v
                                   << " has a " << c << " (int=" << ((int)c) << ") at "
                                   << genomic_region.get_contig_position(label.order, *this).first << ":"
                                   << genomic_region.get_contig_position(label.order, *this).second << " + " << i;
        no_non_ACGTN_errors = false;
      }
      }
    }
  }

  return no_non_ACGTN_errors;
}


bool
Graph::check_empty_variant_dna() const
{
  bool no_empty_variant_dna = true;

  for (long v = 0; v < static_cast<long>(var_nodes.size()); ++v)
  {
    if (var_nodes[v].get_label().dna.size() == 0)
    {
      BOOST_LOG_TRIVIAL(warning) << "[" << __HERE__ << "] Variant node " << v << " has an empty dna sequence.";
      no_empty_variant_dna = false;
    }
  }

  return no_empty_variant_dna;
}


bool
Graph::check_increasing_order() const
{
  bool no_increasing_order_errors = true;

  if (size() == 0)
  {
    assert(var_nodes.size() == 0);
    return true;
  }

  // Reference nodes
  uint32_t old_order = ref_nodes[0].get_label().order;

  for (long r = 1; r < static_cast<long>(ref_nodes.size()); ++r)
  {
    if (ref_nodes[r].get_label().order <= old_order)
    {
      std::cerr << "[graphtyper::graph]: WARNING: Reference node " << r << " has a order smaller or equal than " <<
        old_order << "\n";
      no_increasing_order_errors = false;
    }
  }

  if (var_nodes.size() > 0)
  {
    old_order = var_nodes[0].get_label().order;

    // Variant nodes
    for (long v = 1; v < static_cast<long>(var_nodes.size()); ++v)
    {
      if (var_nodes[v].get_label().order < old_order)
      {
        std::cerr << "[graph]: WARNING: Variant node " << v << " has a order smaller than " << old_order << "\n";
        no_increasing_order_errors = false;
      }
    }
  }

  return no_increasing_order_errors;
}


bool
Graph::check_if_order_follows_reference() const
{
  bool order_follows_reference = true;

  if (ref_nodes.size() <= 1)
  {
    assert(var_nodes.size() == 0);
    return true;
  }

  uint32_t r = 1;
  uint32_t v = 0;

  while (ref_nodes[r].out_degree() > 0)
  {
    assert(r < ref_nodes.size());
    assert(v < var_nodes.size());

    if (ref_nodes[r - 1].get_label().order + ref_nodes[r - 1].get_label().dna.size() != var_nodes[v].get_label().order)
    {
      BOOST_LOG_TRIVIAL(warning) << "[" << __HERE__ << "] Variant orders do not match "
                                 << ref_nodes[r - 1].get_label().order + ref_nodes[r - 1].get_label().dna.size()
                                 << " and "
                                 << std::string(std::begin(var_nodes[v].get_label().dna),
                     std::end(var_nodes[v].get_label().dna));
      order_follows_reference = false;
    }

    if (var_nodes[v].get_label().order + var_nodes[v].get_label().dna.size() !=
        ref_nodes[r].get_label().order
        )
    {
      BOOST_LOG_TRIVIAL(warning) << "[" << __HERE__ << "] Reference orders do not match "
                                 << (var_nodes[v].get_label().order + var_nodes[v].get_label().dna.size())
                                 << " and "
                                 << ref_nodes[r].get_label().order;
      order_follows_reference = false;
    }

    v += ref_nodes[r - 1].out_degree();
    ++r;
  }

  return order_follows_reference;
}


bool
Graph::is_snp(Genotype const & gt) const
{
  long v = gt.first_variant_node;
  assert(v < static_cast<long>(var_nodes.size()));
  auto const & var = var_nodes[v];

  if (var.get_label().dna.size() > 1)
    return false;

  long r = var.get_out_ref_index() - 1; // Go back one ref node
  assert(r >= 0);
  assert(r < static_cast<long>(ref_nodes.size()));
  auto const & ref = ref_nodes[r];

  for (long o = 1; o < static_cast<long>(ref.out_degree()); ++o)
  {
    assert(v + o < static_cast<long>(var_nodes.size()));

    if (var_nodes[v + o].get_label().dna.size() > 1)
      return false;
  }

  return true;
}


bool
Graph::is_homopolymer(Genotype const & gt) const
{
  long v = gt.first_variant_node;
  assert(v < static_cast<long>(var_nodes.size()));
  auto const & var = var_nodes[v];

  // Check if we can find N identical bases
  long constexpr N = 10;
  long same_base = 0;
  std::vector<char> const & seq1 = var.get_label().dna;
  char prev_base = seq1.size() > 0 ? seq1[0] : 'N';

  // Check reference allele
  for (long s = 1; s < static_cast<long>(seq1.size()); ++s)
  {
    if (seq1[s] == prev_base)
      ++same_base;
    else
      same_base = 0;

    if (same_base >= N)
      return true;

    prev_base = seq1[s];
  }

  if (seq1.size() < 50)
  {
    // Check ref nodes behind up to 50 bp total
    long r = var.get_out_ref_index();
    assert(r < static_cast<long>(ref_nodes.size()));
    auto const & ref = ref_nodes[r];
    std::vector<char> const & seq2 = ref.get_label().dna;
    long const LEN = std::min(50l - static_cast<long>(seq1.size()), static_cast<long>(seq2.size()));

    for (long s = 0; s < LEN; ++s)
    {
      if (seq2[s] == prev_base)
        ++same_base;
      else
        same_base = 0;

      if (same_base >= N)
        return true;

      prev_base = seq2[s];
    }
  }


  return false;
}


std::vector<uint32_t>
Graph::get_var_orders(uint32_t const start, uint32_t const end) const
{
  std::vector<uint32_t> var_orders;
  unsigned r = 0;
  unsigned v = 0;

  while (ref_nodes[r].out_degree() != 0)
  {
    if (var_nodes[v].get_label().reach() >= start)
    {
      if (var_nodes[v].get_label().order <= end)
        var_orders.push_back(var_nodes[v].get_label().order);
      else
        return var_orders;
    }

    v += ref_nodes[r].out_degree();
    ++r;
  }

  return var_orders;
}


void
Graph::print() const
{
  long r = 0;
  long v = 0;

  while (ref_nodes[r].out_degree() != 0)
  {
    {
      auto const & label = ref_nodes[r].get_label();
      auto contig_pos = genomic_region.get_contig_position(label.order, *this);
      std::cout << contig_pos.first << ":" << contig_pos.second << " ";

      if (label.dna.size() < 100)
      {
        std::cout << std::string(label.dna.begin(), label.dna.end()) << "\n";
      }
      else
      {
        std::cout << std::string(label.dna.begin(), label.dna.begin() + 40)
                  << "..(size " << label.dna.size() << ").."
                  << std::string(label.dna.end() - 40, label.dna.end())
                  << "\n";
      }
    }

    {
      auto const & var_label = var_nodes[v].get_label();
      auto contig_pos = genomic_region.get_contig_position(var_label.order, *this);
      std::cout << contig_pos.first << ":" << contig_pos.second << " ";
    }

    for (long i = 0; i < static_cast<long>(ref_nodes[r].out_degree()); ++i)
    {
      auto const & var_label = var_nodes[v + i].get_label();
      std::cout << std::string(var_label.dna.begin(), var_label.dna.end()) << "|";
    }

    std::cout << "\n";
    v += ref_nodes[r].out_degree();
    ++r;
  }

  auto const & label = ref_nodes[r].get_label();
  auto contig_pos = genomic_region.get_contig_position(label.order, *this);
  std::cout << contig_pos.first << ":" << contig_pos.second << " ";

  if (label.dna.size() < 100)
  {
    std::cout << std::string(label.dna.begin(), label.dna.end()) << "\n";
  }
  else
  {
    std::cout << std::string(label.dna.begin(), label.dna.begin() + 40)
              << "..(size " << label.dna.size() << ").."
              << std::string(label.dna.end() - 40, label.dna.end())
              << "\n";
  }
  //std::cout << ref_nodes[r].size() << " " r << "\n";
}


/***************************
 * EXPLICIT INSTANTIATIONS *
 ***************************/

template void Graph::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive &, const unsigned int);
template void Graph::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive &, const unsigned int);

/**
 * GLOBAL INSTANCE
 */

Graph graph;

} // namespace gyper

//BOOST_CLASS_VERSION(gyper::Graph, 2)
