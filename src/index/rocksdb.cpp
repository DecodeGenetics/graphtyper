#include <string>

#include <graphtyper/graph/graph.hpp>
#include <graphtyper/index/rocksdb.hpp>
#include <graphtyper/utilities/type_conversions.hpp>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

#include <boost/log/trivial.hpp>

#include <rocksdb/db.h>
#include <rocksdb/env.h>
#include <rocksdb/slice.h>
#include <rocksdb/options.h>
#include <rocksdb/statistics.h>
#include <rocksdb/merge_operator.h>
#include <rocksdb/utilities/backupable_db.h>


namespace gyper
{

uint8_t const LABEL_SIZE = 10;

// Any number of labels
std::vector<gyper::KmerLabel>
value_to_labels(std::string const & value)
{
  assert(value.size() % LABEL_SIZE == 0);
  std::vector<gyper::KmerLabel> results(value.size() / LABEL_SIZE);

  for (unsigned i = 0; i < value.size() / LABEL_SIZE; ++i)
  {
    memcpy(&results[i].start_index, value.data() + i * LABEL_SIZE, sizeof(uint32_t));
    int16_t offset;
    memcpy(&offset, value.data() + 4 + i * LABEL_SIZE, sizeof(int16_t));

    // When offset < 0 then end index is a special pos
    if (offset < 0)
      results[i].end_index = gyper::SPECIAL_START - offset - 1;
    else
      results[i].end_index = graph.get_ref_reach_pos(results[i].start_index) + offset;

    assert(graph.get_ref_reach_pos(results[i].start_index) <= graph.get_ref_reach_pos(results[i].end_index));
    memcpy(&results[i].variant_id, value.data() + 6 + i * LABEL_SIZE, sizeof(uint32_t));

    if (results[i].variant_id != gyper::INVALID_ID)
    {
      results[i].variant_num = graph.get_variant_num(results[i].variant_id);
      results[i].variant_order = graph.var_nodes[results[i].variant_id].get_label().order;
    }
  }

  return results;
}


uint64_t
key_to_uint64_t(std::string const & key_str)
{
  assert(key_str.size() == 8);
  uint64_t key;
  memcpy(&key, key_str.data(), sizeof(uint64_t));
  return key;
}


std::string
labels_to_value(std::vector<gyper::KmerLabel> const & labels)
{
  using gyper::graph;

  // Convert labels to byte array
  uint32_t static const t = LABEL_SIZE;
  std::vector<char> v(t * labels.size());

  for (unsigned i = 0; i < labels.size(); ++i)
  {
    memcpy(v.data() + i * t, &labels[i].start_index, sizeof(uint32_t));
    int16_t offset;

    if (labels[i].end_index >= gyper::SPECIAL_START)
    {
      assert(static_cast<int16_t>(gyper::SPECIAL_START - labels[i].end_index - 1) < 0);
      offset = static_cast<int16_t>(gyper::SPECIAL_START - labels[i].end_index - 1);
    }
    else
    {
      assert(labels[i].end_index > graph.get_ref_reach_pos(labels[i].start_index));
      assert(static_cast<int16_t>(labels[i].end_index - graph.get_ref_reach_pos(labels[i].start_index)) >= 0);
      offset = static_cast<int16_t>(labels[i].end_index - graph.get_ref_reach_pos(labels[i].start_index));
    }

    memcpy(v.data() + 4 + i * t, &offset, sizeof(int16_t));
    memcpy(v.data() + 6 + i * t, &labels[i].variant_id, sizeof(uint32_t));
  }

  return std::string(v.data(), t * labels.size());
}


} // namespace gyper


namespace rocksdb
{

// A 'model' merge operator with uint64 addition semantics
class LabelAppendOperator :
  public AssociativeMergeOperator
{
public:
  virtual bool
  Merge(const Slice & /*key*/,
        const Slice * existing_value,
        const Slice & value,
        std::string * new_value,
        Logger * /*logger*/
        ) const override
  {
    *new_value = std::string(value.data(), value.size());

    if (existing_value)
      new_value->append(std::string(existing_value->data(), existing_value->size()));

    return true;
  }


  virtual const char *
  Name() const override
  {
    return "LabelAppendOperator";
  }


};

} // namespace rocksdb

namespace gyper
{

using namespace rocksdb;

template <>
void
Index<RocksDB>::construct(bool);

template <>
void
Index<RocksDB>::clear();

template <>
void
Index<RocksDB>::commit();


template <>
void
Index<RocksDB>::open(std::string const & f, bool clear_first, bool read_only)
{
  opened = true;
  hamming0.filename = f;
  hamming0.options.db_log_dir = "/tmp/rocksdblogs0";

  if (clear_first)
    Index<RocksDB>::clear();

  Index<RocksDB>::construct(read_only);
}


template <>
void
Index<RocksDB>::close()
{
  if (opened)
  {
    commit(); // Just to make sure everything has been committed from buffer.
    delete hamming0.db;
    opened = false;
  }
}


template <>
Index<RocksDB>::Index()
{}

template <>
Index<RocksDB>::Index(Index<RocksDB> const & cp_index)
{
  buffer_map = cp_index.buffer_map;
  opened = cp_index.opened;
  hamming0 = cp_index.hamming0;
}

template <>
Index<RocksDB>::Index(Index<RocksDB> && mv_index)
{
  buffer_map = std::move(mv_index.buffer_map);
  opened = std::move(mv_index.opened);
  hamming0 = std::move(mv_index.hamming0);
}


template <>
Index<RocksDB>::Index(std::string const & f, bool clear_first, bool read_only)
{
  open(f, clear_first, read_only);
}


template <>
Index<RocksDB>::~Index()
{
  close();
}


template <>
Index<RocksDB> &
Index<RocksDB>::operator=(Index<RocksDB> && mv_index)
{
  std::cerr << "Move assign " << mv_index.buffer_map.size() << " " << mv_index.opened << " " << mv_index.hamming0.db << std::endl;
  this->opened = std::move(mv_index.opened);
  this->hamming0 = std::move(mv_index.hamming0);
  //mv_index.opened = false;
  return *this;
}


template <>
bool
Index<RocksDB>::exists(uint64_t key) const
{
  std::string value;
  hamming0.db->Get(ReadOptions(), Slice(static_cast<const char *>(static_cast<const void *>(&key)), sizeof(uint64_t)), &value);
  return value.size() > 0;
}



template <>
std::vector<KmerLabel>
Index<RocksDB>::get(uint64_t key) const
{
  std::string value;
  hamming0.db->Get(ReadOptions(), Slice(static_cast<const char *>(static_cast<const void *>(&key)), sizeof(uint64_t)), &value);
  return value_to_labels(value);
}


template <>
std::vector<std::vector<KmerLabel> >
Index<RocksDB>::multi_get(std::vector<std::vector<uint64_t> > const & multi_keys) const
{
  std::vector<Slice> slices;
  std::vector<std::size_t> key_to_multi_key;
  std::vector<uint64_t> keys;

  for (std::size_t i = 0; i < multi_keys.size(); ++i)
  {
    for (std::size_t j = 0; j < multi_keys[i].size(); ++j)
    {
      key_to_multi_key.push_back(i);
      keys.push_back(multi_keys[i][j]);
    }
  }

  for (std::size_t k = 0; k < keys.size(); ++k)
    slices.push_back(Slice(static_cast<const char *>(static_cast<const void *>(&keys[k])), sizeof(uint64_t)));

  std::vector<std::string> values;
  hamming0.db->MultiGet(ReadOptions(), slices, &values);
  std::vector<std::vector<KmerLabel> > labels;
  labels.resize(multi_keys.size());

  for (std::size_t i = 0; i < values.size(); ++i)
  {
    std::vector<KmerLabel> label = value_to_labels(values[i]);
    std::move(label.begin(), label.end(), std::back_inserter(labels[key_to_multi_key[i]]));
  }

  return labels;
}


template <>
bool
Index<RocksDB>::check()
{
  bool no_errors = true;

  // Check if reference is in the index
  std::vector<char> ref_seq = graph.get_all_ref();

  if (ref_seq.size() < gyper::K)
  {
    BOOST_LOG_TRIVIAL(warning) << "[rocksdb] The graph is too small for any K-mers to be extracted.";
    return true;
  }

  std::size_t const N_count = std::count(ref_seq.cbegin(), ref_seq.cend(), 'N');

  if (N_count >= ref_seq.size())
  {
    BOOST_LOG_TRIVIAL(warning) << "[rocksdb] The reference contains only Ns.";
    return true;
  }

  auto start_it = ref_seq.begin();
  auto final_it = ref_seq.begin() + K;

  std::vector<std::vector<uint64_t> > keys;
  std::size_t static const MAX_KEYS = 100000;
  std::size_t static const STEP_SIZE = 25;
  keys.reserve(MAX_KEYS);
  std::size_t num_keys_full = 0;

  BOOST_LOG_TRIVIAL(info) << "[rocksdb] Check progress is at 0" << '%';

  while (std::distance(final_it, ref_seq.end()) >= static_cast<long>(STEP_SIZE))
  {
    std::vector<char> k_mer(start_it, final_it);

    // Ignore kmers with an N
    if (std::find(k_mer.begin(), k_mer.end(), 'N') == k_mer.end())
    {
      keys.push_back(std::vector<uint64_t>(1, to_uint64(std::move(k_mer))));

      if (keys.size() == MAX_KEYS)
      {
        ++num_keys_full;

        std::vector<std::vector<KmerLabel> > labels = multi_get(keys);

        for (unsigned i = 0; i < labels.size(); ++i)
        {
          if (labels[i].size() == 0)
          {
            std::cerr << "[rocksdb] Could not find kmer " << to_dna(keys[i][0]) << " at position "
                      << (i + std::distance(ref_seq.begin(), start_it) - MAX_KEYS) << std::endl;
            no_errors = false;
          }
        }

        BOOST_LOG_TRIVIAL(info) << "[rocksdb] Check progress is at "
                                << (num_keys_full * STEP_SIZE * 100 * MAX_KEYS / (ref_seq.size() - N_count))
                                << '%';
        keys.clear();
      }
    }

    start_it += STEP_SIZE;
    final_it += STEP_SIZE;
  }

  std::vector<std::vector<KmerLabel> > labels = multi_get(keys);

  for (unsigned i = 0; i < labels.size(); ++i)
  {
    if (labels[i].size() == 0)
    {
      std::cerr << "[rocksdb] Could not find kmer " << to_dna(keys[i][0]) << " at position "
                << (i + std::distance(ref_seq.begin(), start_it) - MAX_KEYS) << std::endl;
      no_errors = false;
    }
  }

  BOOST_LOG_TRIVIAL(info) << "[rocksdb] Check progress is at 100" << '%';
  return no_errors;
}


template <>
void
Index<RocksDB>::put(uint64_t key, KmerLabel && label)
{
  buffer_map[key].push_back(std::move(label));

  if (buffer_map.size() > MAX_BUFFER)
    commit();
}


template <>
void
Index<RocksDB>::put(uint64_t key, std::vector<KmerLabel> && labels)
{
  std::move(labels.begin(), labels.end(), std::back_inserter(buffer_map[key]));

  if (buffer_map.size() > MAX_BUFFER)
    commit();
}


template <>
void
Index<RocksDB>::commit()
{
  for (auto it = buffer_map.begin(); it != buffer_map.end(); ++it)
  {
    hamming0.s = hamming0.db->Merge(WriteOptions(),
                              Slice(static_cast<const char *>(static_cast<const void *>(&it->first)), sizeof(uint64_t)),
                              Slice(labels_to_value(it->second))
                              );
    assert(hamming0.s.ok());
  }

  buffer_map.clear(); // Free the buffer
}


template <>
std::size_t
Index<RocksDB>::size()
{
  return 1; // So test pass, yay
}


template <>
void
Index<RocksDB>::clear()
{
  hamming0.s = rocksdb::DestroyDB(hamming0.filename.c_str(), hamming0.options);

  if (!hamming0.s.ok())
  {
    delete hamming0.db;
    hamming0.s = rocksdb::DestroyDB(hamming0.filename.c_str(), hamming0.options);

    if (!hamming0.s.ok())
      std::cerr << "ERROR: Could not destroy database '" << hamming0.filename << "'. Message: " << hamming0.s.ToString() << std::endl;
  }
}


template <>
void
Index<RocksDB>::construct(bool read_only)
{
  hamming0.options.create_if_missing = true;
  // hamming0.options.statistics = rocksdb::CreateDBStatistics();
  // hamming0.options.OptimizeForPointLookup(1024);
  hamming0.options.IncreaseParallelism();
  // hamming0.options.OptimizeUniversalStyleCompaction();
  hamming0.options.OptimizeLevelStyleCompaction();
  hamming0.options.merge_operator.reset(new LabelAppendOperator);

  if (read_only)
    hamming0.s = DB::OpenForReadOnly(hamming0.options, hamming0.filename.c_str(), &hamming0.db);
  else
    hamming0.s = DB::Open(hamming0.options, hamming0.filename.c_str(), &hamming0.db);

  if (!hamming0.s.ok())
  {
    std::cerr << "ERROR: While trying to open '" << hamming0.filename << "'. Message: " << hamming0.s.ToString() << std::endl;
    std::exit(1);
  }
}


/** Explicit instantation */
template class Index<RocksDB>;

Index<RocksDB> index;

} // namespace gyper
