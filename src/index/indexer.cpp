#include <algorithm>
#include <fstream>
#include <iostream>

#include <boost/log/trivial.hpp>

#include <graphtyper/graph/graph.hpp>
#include <graphtyper/graph/graph_serialization.hpp>
#include <graphtyper/index/indexer.hpp>
#include <graphtyper/index/ph_index.hpp>


/*
namespace
{

bool
entry_has_too_many_nonrefs(gyper::IndexEntry const & entry)
{
  uint32_t constexpr MAX_TOTAL_VAR_NUM = 401u;
  return entry.total_var_count > 1 && entry.total_var_num > MAX_TOTAL_VAR_NUM;
}


} // anon namespace
*/

namespace gyper
{


void
index_reference_label(PHIndex & ph_index, TEntryList & mers, Label const & label)
{
  for (int d{0}; d < static_cast<int>(label.dna.size()); ++d)
  {
    char const dna_base = label.dna[d];

    if (dna_base != 'A' && dna_base != 'C' && dna_base != 'G' && dna_base != 'T')
    {
      mers.clear();
      continue;
    }

    for (auto list_it = mers.begin(); list_it != mers.end(); ++list_it)
    {
      for (auto sublist_it = list_it->begin(); sublist_it != list_it->end(); ++sublist_it)
      {
        sublist_it->add_to_dna(dna_base);
      }
    }

    // Add a new element with the new DNA base
    {
      IndexEntry index_entry(label.order + d);
      index_entry.add_to_dna(dna_base);
      mers.push_front(TEntrySublist(1, index_entry));
    }

    if (mers.size() >= K)
    {
      for (auto q_it = mers.back().begin(); q_it != mers.back().end(); ++q_it)
      {
        // Skip invalid labels (e.g. labels with '*')
        if (q_it->valid > 0)
          continue;

        if (q_it->variant_id.size() == 0)
        {
          // KmerLabel has implicit var_id = INVALID_ID
          ph_index.put(q_it->dna, KmerLabel(q_it->start_index, label.order + d));
        }
        else
        {
          std::vector<KmerLabel> new_labels;
          new_labels.reserve(q_it->variant_id.size());

          for (auto var_id : q_it->variant_id)
            new_labels.push_back(KmerLabel(q_it->start_index, label.order + d, var_id));

          ph_index.put(q_it->dna, std::vector<KmerLabel>(new_labels));
        }
      }

      mers.pop_back();
    }
  }
}


void
insert_variant_label(PHIndex & ph_index,
                     TEntryList & mers,
                     std::vector<VarNode> const & var_nodes,
                     TNodeIndex const v,
                     std::size_t const ref_reach)
{
  assert(v < var_nodes.size());
  VarNode const & var_node = var_nodes[v];
  Label const & label = var_node.get_label();

  for (long d{0}; d < static_cast<int>(label.dna.size()); ++d)
  {
    char const dna_base = label.dna[d];
    assert(dna_base != 'N'); // shouldn't happen

    if (dna_base != 'A' && dna_base != 'C' && dna_base != 'G' && dna_base != 'T')
    {
      mers.clear();
      continue;
    }

    for (auto sublist_it = mers.begin(); sublist_it != mers.end(); ++sublist_it)
    {
      for (auto entry_it = sublist_it->begin(); entry_it != sublist_it->end();) // no increment of ++entry_it)
      {
        // Check if we should add the base
        IndexEntry & entry = *entry_it;
        bool is_ok{true};

        //if (var_node.events.size() > 0 || var_node.anti_events.size() > 0)
        //{
        //  BOOST_LOG_TRIVIAL(info) << __HERE__ << " " << var_node.events.size() << " "
        //                          << var_node.anti_events.size() << " "
        //                          << entry.events.size() << " "
        //                          << entry.anti_events.size();
        //}

        for (auto const & anti_event : entry.anti_events)
        {
          if (var_node.events.count(anti_event) == 1)
          {
            //BOOST_LOG_TRIVIAL(info) << __HERE__ << " bad event " << anti_event;
            is_ok = false;
            break;
          }
        }

        if (is_ok)
        {
          entry.add_to_dna(dna_base);

          std::copy(var_node.events.begin(),
                    var_node.events.end(),
                    std::inserter(entry.events, entry.events.begin()));

          std::copy(var_node.anti_events.begin(),
                    var_node.anti_events.end(),
                    std::inserter(entry.anti_events, entry.anti_events.begin()));

          entry.variant_id.insert(v);
          ++entry_it;
        }
        else
        {
          entry_it = sublist_it->erase(entry_it);
        }

        /*
        for (auto const & event : entry.events)
        {
          if (var_node.anti_events.count(event) == 1)
          {
            BOOST_LOG_TRIVIAL(info) << __HERE__ << " bad event " << event.to_string();
            is_ok = false;
          }
        }*/
      }
    }

    // Add a new element with the new DNA base
    uint32_t pos = label.order + d;

    if (pos > ref_reach)
      pos = graph.get_special_pos(pos, static_cast<uint32_t>(ref_reach));

    IndexEntry new_index_entry(pos, static_cast<uint32_t>(v));
    new_index_entry.add_to_dna(dna_base);
    new_index_entry.events = var_node.events;
    new_index_entry.anti_events = var_node.anti_events;

    // If we are using a list
    mers.push_front(TEntrySublist(1, new_index_entry));

    if (mers.size() >= K)
    {
      // Insert to map
      for (auto q_it = mers.back().begin(); q_it != mers.back().end(); ++q_it)
      {
        // Skip invalid labels (e.g. labels with '*')
        if (q_it->valid > 0)
          continue;

        std::vector<KmerLabel> new_labels;
        new_labels.reserve(q_it->variant_id.size());

        for (auto var_id : q_it->variant_id)
          new_labels.push_back(KmerLabel(q_it->start_index, pos, var_id));

        ph_index.put(q_it->dna, std::vector<KmerLabel>(new_labels));
      }

      mers.pop_back();
    }
  }
}


void
append_list(TEntryList & mers, TEntryList && list)
{
  if (mers.size() < list.size())
    mers.resize(list.size());

  auto mer_it = mers.begin();
  auto list_it = list.begin();

  while (list_it != list.end())
  {
    std::move(list_it->begin(), list_it->end(), std::back_inserter(*mer_it));
    ++list_it;
    ++mer_it;
  }

  assert(list.size() <= mers.size());
}


/*
void
remove_large_variants_from_list(TEntryList & list, unsigned const var_count)
{
  for (auto sublist_it = list.begin(); sublist_it != list.end(); ++sublist_it)
  {
    for (auto entry_it = sublist_it->begin(); entry_it != sublist_it->end(); ++entry_it)
    {
      entry_it->total_var_num *= var_count;
      ++entry_it->total_var_count;
    }

    sublist_it->erase(std::remove_if(sublist_it->begin(),
                                     sublist_it->end(),
                                     entry_has_too_many_nonrefs),
                      sublist_it->end());
  }
}*/


void
index_variant(PHIndex & ph_index,
              std::vector<VarNode> const & var_nodes,
              TEntryList & mers,
              unsigned var_count,
              TNodeIndex v)
{
  TEntryList clean_list(mers); // copies all mers, we find new kmers using the copy.

  // Insert reference label
  assert(v < var_nodes.size());
  std::size_t const ref_label_reach = var_nodes[v].get_label().reach();
  insert_variant_label(ph_index, mers, var_nodes, v, ref_label_reach);

  // Remove all labels with large variants
  //remove_large_variants_from_list(clean_list, var_count);
  //unsigned const var_num = var_count;

  // Loops over variants
  while (var_count > 2)
  {
    --var_count;
    ++v;

    TEntryList new_list(clean_list); // copies all mers, we find new kmers using the new copy
    insert_variant_label(ph_index,
                         new_list,
                         var_nodes,
                         v,
                         ref_label_reach);

    append_list(mers, std::move(new_list));
  }

  // No need to copy clean_list on the last variant
  ++v;
  insert_variant_label(ph_index,
                       clean_list,
                       var_nodes,
                       v,
                       ref_label_reach);

  append_list(mers, std::move(clean_list));
}


PHIndex
index_graph(Graph const & graph)
{
  PHIndex ph_index;

  assert(graph.ref_nodes.back().out_degree() == 0);
  uint32_t const start_order = graph.ref_nodes.front().get_label().order;
  uint32_t const end_order = static_cast<uint32_t>(graph.ref_nodes.back().get_label().order +
                                                   graph.ref_nodes.back().get_label().dna.size());
  uint32_t goal_order = start_order;
  uint32_t goal = 0;

  TNodeIndex r{0}; // Reference node index
  TEntryList mers;
  BOOST_LOG_TRIVIAL(debug) << __HERE__ << " The number of reference nodes are "
                           << graph.ref_nodes.size();

  while (r < graph.ref_nodes.size() - 1)
  {
    if (graph.ref_nodes[r].get_label().order >= goal_order)
    {
      BOOST_LOG_TRIVIAL(debug) << __HERE__ << " Indexing progress: " << goal << '%';
      goal_order += (end_order - start_order) / 5;
      goal += 20;
    }

    index_reference_label(ph_index, mers, graph.ref_nodes[r].get_label());

    if (graph.ref_nodes[r].out_degree() > 0)
    {
      index_variant(ph_index,
                    graph.var_nodes,
                    mers,
                    static_cast<int>(graph.ref_nodes[r].out_degree()),
                    graph.ref_nodes[r].get_var_index(0));
    }

    ++r;
  }

  index_reference_label(ph_index, mers, graph.ref_nodes.back().get_label());
  BOOST_LOG_TRIVIAL(debug) << __HERE__ << " Indexing progress: 100%";
  mers.clear();

  // Commit the rest of the buffer before closing
  BOOST_LOG_TRIVIAL(debug) << __HERE__ << " Writing index to disk...";
  BOOST_LOG_TRIVIAL(debug) << __HERE__ << " Done indexing graph.";
  return ph_index;
}


PHIndex
index_graph(std::string const & graph_path)
{
  load_graph(graph_path);

  if (graph.size() == 0)
  {
    BOOST_LOG_TRIVIAL(warning) << __HERE__ << " Trying to index empty graph.";
    PHIndex ph_index; // Make empty index
    return ph_index;
  }

  return index_graph(graph);
}


} // namespace gyper
