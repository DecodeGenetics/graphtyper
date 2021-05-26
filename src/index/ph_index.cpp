#include <string>
#include <vector>

#include <graphtyper/utilities/logging.hpp>

#include <parallel_hashmap/phmap.h>

#include <graphtyper/graph/graph.hpp>
#include <graphtyper/index/ph_index.hpp>
#include <graphtyper/utilities/options.hpp>
#include <graphtyper/utilities/type_conversions.hpp>


namespace gyper
{

void
PHIndex::put(uint64_t const key, KmerLabel && label)
{
  hamming0[key].push_back(std::move(label));
}


void
PHIndex::put(uint64_t const key, std::vector<KmerLabel> && labels)
{
  std::move(labels.begin(), labels.end(), std::back_inserter(hamming0[key]));
}


std::vector<KmerLabel>
PHIndex::get(uint64_t const key) const
{
  auto find_it = hamming0.find(key);

  if (find_it != hamming0.end())
    return find_it->second;
  else
    return std::vector<KmerLabel>();
}


std::vector<KmerLabel>
PHIndex::get(std::vector<uint64_t> const & keys) const
{
  std::vector<KmerLabel> labels;
  std::vector<PHtype::const_iterator> results;
  long num_results{0};
  long const NUM_KEYS = keys.size();

  for (long j = 0; j < NUM_KEYS; ++j)
  {
    auto find_it = hamming0.find(keys[j]);

    if (find_it != hamming0.end())
    {
      num_results += find_it->second.size();

      if (NUM_KEYS > 1 && num_results > Options::const_instance()->max_index_labels)
      {
        // Too many results, give up on this kmer
        results.clear();
        break;
      }

      results.push_back(find_it);
    }
  }

  for (auto const res : results)
    std::copy(res->second.begin(), res->second.end(), std::back_inserter(labels));

  return labels;
}


std::vector<std::vector<KmerLabel> >
PHIndex::multi_get(std::vector<std::vector<uint64_t> > const & keys) const
{
  std::vector<std::vector<KmerLabel> > labels(keys.size());
  std::vector<std::vector<PHtype::const_iterator> > results(keys.size());

  for (long i = 0; i < static_cast<long>(keys.size()); ++i)
  {
    long num_results{0};
    long const NUM_KEYS{static_cast<long>(keys[i].size())};

    for (long j = 0; j < NUM_KEYS; ++j)
    {
      auto find_it = hamming0.find(keys[i][j]);

      if (find_it != hamming0.end())
      {
        num_results += find_it->second.size();

        if (NUM_KEYS > 1 && num_results > Options::const_instance()->max_index_labels)
        {
          // Checking multiple keys and finding too many results, give up on this
          results[i].clear();
          break;
        }

        results[i].push_back(find_it);
      }
    }
  }

  assert(labels.size() == keys.size());
  assert(results.size() == keys.size());

  // If there are not too many results, add them to labels and return them
  for (long i = 0; i < static_cast<long>(results.size()); ++i)
  {
    for (auto const res : results[i])
      std::copy(res->second.begin(), res->second.end(), std::back_inserter(labels[i]));
  }

  return labels;
}


void
PHIndex::print() const
{
  for (auto it = hamming0.begin(); it != hamming0.end(); ++it)
  {
    std::cerr << to_dna_str(it->first) << " count=" << it->second.size() << '\n';
  }
}


/*
std::vector<std::vector<KmerLabel> >
PHIndex::multi_get_hamming1(std::vector<std::vector<uint64_t> > const & keys) const
{
  std::vector<std::vector<KmerLabel> > labels(keys.size());

  for (long i = 0; i < static_cast<long>(keys.size()); ++i)
  {
    long const NUM_KEYS = keys[i].size();

    for (long j = 0; j < NUM_KEYS; ++j)
    {
      auto find_it = hamming1.find(keys[i][j]);

      if (find_it != hamming1.end())
      {
        auto hamming0_find_it = hamming0.find(find_it->second);
        assert(hamming0_find_it != hamming0.end());
        std::copy(hamming0_find_it->second.begin(), hamming0_find_it->second.end(), std::back_inserter(labels[i]));
      }
    }
  }

  assert(keys.size() == labels.size());
  return labels;
}
*/


bool
PHIndex::check() const
{
  bool no_errors = true;

  // Check if reference is in the index
  std::vector<char> ref_seq = graph.get_all_ref();

  if (ref_seq.size() < gyper::K)
  {
    print_log(log_severity::warning, "[", __HERE__, "] The graph is too small for any K-mers to be extracted.");
    return true;
  }

  if (std::all_of(ref_seq.cbegin(), ref_seq.cend(), [](char c){
      return c == 'N';
    }))
  {
    print_log(log_severity::warning, "[", __HERE__, "] The reference contains only Ns.");
    return true;
  }

  auto start_it = ref_seq.begin();
  auto final_it = ref_seq.begin() + K;

  std::vector<std::vector<uint64_t> > keys;
  std::size_t constexpr MAX_KEYS = 100000;
  long constexpr STEP_SIZE = 5;
  keys.reserve(MAX_KEYS);
  std::size_t num_keys_full = 0;

  print_log(log_severity::debug, "[", __HERE__, "] Checking index");

  while (std::distance(final_it, ref_seq.end()) >= STEP_SIZE)
  {
    std::vector<char> k_mer(start_it, final_it);
    assert(k_mer.size() == K);

    // Ignore kmers with non ACGT bases
    if (std::all_of(k_mer.begin(), k_mer.end(), [](char c){
        return c == 'A' || c == 'C' || c == 'G' || c == 'T';
      }))
    {
      keys.push_back(std::vector<uint64_t>(1, to_uint64(std::move(k_mer))));

      if (keys.size() == MAX_KEYS)
      {
        ++num_keys_full;

        std::vector<std::vector<KmerLabel> > labels = multi_get(keys);

        for (int i = 0; i < static_cast<int>(labels.size()); ++i)
        {
          if (labels[i].size() == 0)
          {
            assert(keys[i].size() == 1);
            print_log(log_severity::error, "[", __HERE__, "] Could not find kmer at position "
                                    ,std::string(start_it, final_it), " "
                                    , (i + std::distance(ref_seq.begin(), start_it)));
            no_errors = false;
          }
        }

        keys.clear();
      }
    }

    start_it += STEP_SIZE;
    final_it += STEP_SIZE;
  }

  std::vector<std::vector<KmerLabel> > labels = multi_get(keys);

  for (int i = 0; i < static_cast<int>(labels.size()); ++i)
  {
    if (labels[i].size() == 0)
    {
      assert(keys[i].size() == 1);
      print_log(log_severity::error, "[", __HERE__, "] Could not find kmer at position "
                              , to_dna_str(keys[i][0]), " "
                              , (std::distance(ref_seq.begin(), start_it)), " with i=", i);
      no_errors = false;
    }
  }

  print_log(log_severity::debug, "[", __HERE__, "] DONE Checking index");
  return no_errors;
}


} // namespace gyper
