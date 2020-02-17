#pragma once

#include <fstream> // std::ifstream
#include <string> // std::string
#include <vector> // std::vector

#include <htslib/sam.h> // bam1_t


namespace gyper
{

using Tqueue = std::pair<bam1_t*, int>;


inline std::vector<std::string> get_lines(std::string const & path)
{
  std::vector<std::string> lines;
  std::string line;
  std::ifstream is(path);

  while (std::getline(is, line))
    lines.push_back(line);

  return lines;
}


inline std::string
get_std_string_seq(const bam1_t* rec)
{
  int32_t const lqseq = rec->core.l_qseq;
  std::string seq;
  seq.resize(lqseq);
  uint8_t* seqptr = bam_get_seq(rec);

  for (int32_t i = 0; i < lqseq; ++i)
  {
    seq[i] = seq_nt16_str[bam_seqi(seqptr, i)];
  }

  return seq;
}



/**
 * Comparison functions
 */
inline bool
gt_pos(bam1_t* a, bam1_t* b)
{
  return a->core.tid > b->core.tid || (a->core.tid == b->core.tid && a->core.pos > b->core.pos);
}


inline bool
gt_pos_seq(bam1_t* a, bam1_t* b)
{
  if (a->core.tid > b->core.tid)
    return true;
  else if (b->core.tid > a->core.tid)
    return false;
  else if (a->core.pos > b->core.pos)
    return true;
  else if (b->core.pos > a->core.pos)
    return false;
  else if (a->core.l_qseq > b->core.l_qseq)
    return true;
  else if (b->core.l_qseq > a->core.l_qseq)
    return false;

  // Check the sequences
  int32_t const lqseq = (a->core.l_qseq + 1) / 2;
  uint8_t* a_seqptr = bam_get_seq(a);
  uint8_t* b_seqptr = bam_get_seq(b);

  for (int32_t i = 0; i < lqseq; ++i)
  {
    auto a_seq = *(a_seqptr + i);
    auto b_seq = *(b_seqptr + i);

    if (a_seq > b_seq)
      return true;
    else if (b_seq > a_seq)
      return false;
    // and continue if they are the same
  }

  return false;
}


inline bool
gt_pos_seq_same_pos(bam1_t* a, bam1_t* b)
{
  if (a->core.l_qseq > b->core.l_qseq)
    return true;
  else if (b->core.l_qseq > a->core.l_qseq)
    return false;

  // Check the sequences
  int32_t const lqseq = (a->core.l_qseq + 1) / 2;
  uint8_t* a_seqptr = bam_get_seq(a);
  uint8_t* b_seqptr = bam_get_seq(b);

  for (int32_t i = 0; i < lqseq; ++i)
  {
    auto a_seq = *(a_seqptr + i);
    auto b_seq = *(b_seqptr + i);

    if (a_seq > b_seq)
      return true;
    else if (b_seq > a_seq)
      return false;
    // and continue if they are the same
  }

  return false;
}


inline bool
equal_pos_seq(bam1_t* a, bam1_t* b)
{
  // Check the positions and the lengdth of the sequences of the records
  if (a->core.tid != b->core.tid || a->core.pos != b->core.pos || a->core.l_qseq != b->core.l_qseq)
    return false;

  // Check the sequences
  long const lqseq = (a->core.l_qseq + 1l) / 2l;
  uint8_t* a_seqptr = bam_get_seq(a);
  uint8_t* b_seqptr = bam_get_seq(b);

  for (auto i = 0; i < lqseq; ++i)
  {
    if (*(a_seqptr + i) != *(b_seqptr + i))
      return false;
  }

  return true;
}


inline bool
equal_pos(bam1_t* a, bam1_t* b)
{
  // Check the positions and the lengdth of the sequences of the records
  return a->core.tid == b->core.tid && a->core.pos == b->core.pos;
}


} // namespace hts
