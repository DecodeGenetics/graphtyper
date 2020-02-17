#include <algorithm>
#include <fstream>
#include <iostream>

#include <boost/log/trivial.hpp>

#include <graphtyper/graph/graph.hpp>
#include <graphtyper/graph/graph_serialization.hpp>
#include <graphtyper/index/indexer.hpp>

#include <seqan/stream.h>


namespace
{

bool
entry_has_too_many_nonrefs(gyper::IndexEntry const & entry)
{
  uint32_t const MAX_TOTAL_VAR_NUM = 401u; // 131 is a prime
  uint32_t const MAX_TOTAL_VAR_COUNT = 4u;
  return entry.total_var_count > 1 &&
         (entry.total_var_num > MAX_TOTAL_VAR_NUM || entry.total_var_count > MAX_TOTAL_VAR_COUNT);
}


} // anon namespace


namespace gyper
{


void
index_reference_label(Index<RocksDB> & new_index, TEntryList & mers, Label const & label)
{
  for (unsigned d = 0; d < seqan::length(label.dna); ++d)
  {
    if (label.dna[d] == 'N')
    {
      mers.clear();
      continue;
    }

    for (auto list_it = mers.begin(); list_it != mers.end(); ++list_it)
    {
      for (auto sublist_it = list_it->begin(); sublist_it != list_it->end(); ++sublist_it)
      {
        sublist_it->add_to_dna(label.dna[d]);
      }
    }

    // Add a new element with the new DNA base
    {
      IndexEntry index_entry(label.order + d);
      index_entry.add_to_dna(label.dna[d]);
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
          new_index.put(q_it->dna, KmerLabel(q_it->start_index, label.order + d)); // KmerLabel has implicit var_id = INVALID_ID
        }
        else
        {
          std::vector<KmerLabel> new_labels;
          new_labels.reserve(q_it->variant_id.size());

          for (unsigned i = 0; i < q_it->variant_id.size(); ++i)
            new_labels.push_back(KmerLabel(q_it->start_index, label.order + d, q_it->variant_id[i]));

          new_index.put(q_it->dna, std::move(new_labels));
        }
      }

      mers.pop_back();
    }
  }
}


void
insert_variant_label(Index<RocksDB> & new_index,
                     TEntryList & mers,
                     Label const & label,
                     TNodeIndex const v,
                     bool const is_reference,
                     unsigned const var_count,
                     std::size_t const ref_reach
                     )
{
  for (unsigned d = 0; d < seqan::length(label.dna); ++d)
  {
    for (auto sublist_it = mers.begin(); sublist_it != mers.end(); ++sublist_it)
    {
      for (auto entry_it = sublist_it->begin(); entry_it != sublist_it->end(); ++entry_it)
      {
        entry_it->add_to_dna(label.dna[d]);

        if (std::find(entry_it->variant_id.begin(), entry_it->variant_id.end(), v) == entry_it->variant_id.end())
          entry_it->variant_id.push_back(static_cast<unsigned>(v));
      }
    }

    // Add a new element with the new DNA base
    uint32_t pos = label.order + d;

    if (pos > ref_reach)
      pos = graph.get_special_pos(pos, static_cast<uint32_t>(ref_reach));

    IndexEntry new_index_entry(pos, static_cast<uint32_t>(v), is_reference, var_count);
    new_index_entry.add_to_dna(label.dna[d]);

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

        for (unsigned i = 0; i < q_it->variant_id.size(); ++i)
          new_labels.push_back(KmerLabel(q_it->start_index, pos, q_it->variant_id[i]));

        new_index.put(q_it->dna, std::move(new_labels));
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

    sublist_it->erase(std::remove_if(sublist_it->begin(), sublist_it->end(), entry_has_too_many_nonrefs),
                      sublist_it->end());
  }
}


void
index_variant(Index<RocksDB> & new_index,
              std::vector<VarNode> const & var_nodes,
              TEntryList & mers,
              unsigned var_count,
              TNodeIndex v
              )
{
  TEntryList clean_list(mers); // copies all mers, we find new kmers using the copy.

  // Insert reference label
  std::size_t const ref_label_reach = var_nodes[v].get_label().reach();
  insert_variant_label(new_index, mers, var_nodes[v].get_label(), v, true /*is reference*/, 1, ref_label_reach);

  // Remove all labels with large variants
  remove_large_variants_from_list(clean_list, var_count);
  unsigned const var_num = var_count;

  // Loops over variants
  while (var_count > 2)
  {
    --var_count;
    ++v;

    TEntryList new_list(clean_list); // copies all mers, we find new kmers using the new copy
    insert_variant_label(new_index,
                         new_list,
                         var_nodes[v].get_label(),
                         v,
                         false /*is reference*/,
                         var_num,
                         ref_label_reach);
    append_list(mers, std::move(new_list));
  }

  // No need to copy clean_list on the last variant
  ++v;
  insert_variant_label(new_index,
                       clean_list,
                       var_nodes[v].get_label(),
                       v,
                       false /*is reference*/,
                       var_num,
                       ref_label_reach);
  append_list(mers, std::move(clean_list));
}


void
index_graph(std::string const & index_path)
{
  Index<RocksDB> new_index(index_path, true /*clear_first*/, false /*read_only*/);

  assert(graph.ref_nodes.back().out_degree() == 0);
  uint32_t const start_order = graph.ref_nodes.front().get_label().order;
  uint32_t const end_order = static_cast<uint32_t>(graph.ref_nodes.back().get_label().order +
                                                   graph.ref_nodes.back().get_label().dna.size());
  uint32_t goal_order = start_order;
  uint32_t goal = 0;

  TNodeIndex r = 0; // Reference node index
  // TNodeIndex v = 0; // Variant node index
  TEntryList mers;
  BOOST_LOG_TRIVIAL(debug) << "[graphtyper::indexer] The number of reference nodes are " << graph.ref_nodes.size();

  while (r < graph.ref_nodes.size() - 1)
  {
    if (graph.ref_nodes[r].get_label().order >= goal_order)
    {
      BOOST_LOG_TRIVIAL(debug) << "[graphtyper::indexer] Indexing progress: " << goal << '%';
      goal_order += (end_order - start_order) / 5;
      goal += 20;
    }

    index_reference_label(new_index, mers, graph.ref_nodes[r].get_label());

    if (graph.ref_nodes[r].out_degree() > 0)
    {
      index_variant(new_index,
                    graph.var_nodes,
                    mers,
                    static_cast<int>(graph.ref_nodes[r].out_degree()),
                    graph.ref_nodes[r].get_var_index(0)
                    );
    }

    ++r;
  }

  index_reference_label(new_index, mers, graph.ref_nodes.back().get_label());
  BOOST_LOG_TRIVIAL(debug) << "[graphtyper::indexer] Indexing progress: 100" << '%';
  mers.clear();

  // Commit the rest of the buffer before closing
  BOOST_LOG_TRIVIAL(debug) << "[graphtyper::indexer] Writing index to disk...";
  new_index.commit();
  BOOST_LOG_TRIVIAL(debug) << "[graphtyper::indexer] Done.";
}


void
index_graph(std::string const & graph_path, std::string const & index_path)
{
  load_graph(graph_path);

  if (graph.size() == 0)
  {
    BOOST_LOG_TRIVIAL(warning) << "[graphtyper::indexer] WARNING: Trying to index empty graph.";
    return;
  }

  index_graph(index_path);
}


void
load_index(std::string const & index_path)
{
  index.open(index_path.c_str(), false /*clear_first*/, true /*read_only*/);
}


Index<RocksDB>
load_secondary_index(std::string const & index_path)
{
  Index<RocksDB> secondary_index;
  secondary_index.open(index_path.c_str(), false /*clear_first*/, true /*read_only*/);
  return secondary_index;
}


} // namespace gyper
