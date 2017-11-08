#pragma once
#include <vector>
#include <iostream>

#include <graphtyper/constants.hpp>

namespace gyper
{

class SegmentCall
{
public:
  std::vector<uint8_t> phred;
  std::pair<uint16_t, uint16_t> call;
};


class Segment
{
public:
  uint32_t id;
  std::size_t ref_size;
  std::vector<std::string> allele_names;
  std::vector<SegmentCall> segment_calls;
  std::string segment_name;
  std::string var_type;
  int32_t extra_id = -1;

  Segment();
  Segment(uint32_t _id, std::size_t ref_size, std::vector<std::string> const & _alts);
  void clear();
  void insert_score(std::vector<uint32_t> const & score);
  void insert_scores(std::vector<std::vector<uint32_t> > const & scores);

  std::string get_ref_string() const;
  std::string get_alt_string() const;
  std::vector<Segment> get_biallelic_segments() const;
};

} // namespace gyper
