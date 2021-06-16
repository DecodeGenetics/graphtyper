#pragma once

#include <vector>

#include <htslib/sam.h>

namespace gyper
{
class HtsStore
{
  std::vector<bam1_t *> store;

public:
  HtsStore() = default;

  HtsStore(HtsStore const &) = delete;
  HtsStore(HtsStore &&) = delete;
  HtsStore & operator=(HtsStore const &) = delete;
  HtsStore & operator=(HtsStore &&) = delete;

  ~HtsStore()
  {
    clear();
  }

  inline void push(bam1_t * a)
  {
    if (store.size() <= 10000)
      store.push_back(a);
    else
      bam_destroy1(a);
  }

  inline bam1_t * get()
  {
    if (store.size() == 0)
      return bam_init1();

    bam1_t * rec = *(store.end() - 1);
    store.pop_back();
    return rec;
  }

  inline void clear()
  {
    for (auto s : store)
      bam_destroy1(s);
  }
};

} // namespace gyper
