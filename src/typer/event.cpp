#include <string>
#include <vector>

#include <boost/log/trivial.hpp>

#include <graphtyper/constants.hpp>
#include <graphtyper/typer/event.hpp>
#include <graphtyper/utilities/options.hpp>


namespace gyper
{

long
base2index(char const base)
{
  switch (base)
  {
  case 'A':
  {
    return 0;
  }

  case 'C':
  {
    return 1;
  }

  case 'G':
  {
    return 2;
  }

  case 'T':
  {
    return 3;
  }

  case 'N':
  {
    return -1;
  }

  case '*':
  {
    return -2;
  }

  default:
  {
    BOOST_LOG_TRIVIAL(error) << __HERE__ << " Unexpected base '" << base << "'";
    std::exit(1);
  }
  }   // switch ends
}


long
BaseCount::get_depth_without_deleted() const
{
  return acgt[0] + acgt[1] + acgt[2] + acgt[3];
}


long
BaseCount::get_depth_with_deleted() const
{
  return acgt[0] + acgt[1] + acgt[2] + acgt[3] + deleted;
}


std::string
BaseCount::to_string() const
{
  std::ostringstream ss;
  ss << acgt[0] << "," << acgt[1] << "," << acgt[2] << "," << acgt[3] << ":"
     << acgt_qualsum[0] << "," << acgt_qualsum[1] << "," << acgt_qualsum[2]
     << "," << acgt_qualsum[3] << " " << deleted;
  return ss.str();
}


void
BaseCount::add_base(char seq, char qual)
{
  long i = base2index(seq);

  if (i < 0)
  {
    if (i == -1)
    {
      ++unknown;
      return;
    }
    else
    {
      assert(i == -2);
      // -2 is deleted base
      ++deleted;
      return;
    }
  }

  assert(i < 4);
  ++(acgt[i]);
  acgt_qualsum[i] += static_cast<int64_t>(qual);
}


bool
IndelEvent::operator==(IndelEvent const & b) const
{
  return (pos == b.pos) &&
         (type == b.type) &&
         (sequence == b.sequence);
}


bool
IndelEvent::operator!=(IndelEvent const & b) const
{
  return !(*this == b);
}


bool
IndelEvent::operator<(IndelEvent const & b) const
{
  return pos < b.pos || (pos == b.pos && type < b.type) || (pos == b.pos && type == b.type && sequence < b.sequence);
}


bool
SnpEvent::operator==(SnpEvent const & b) const
{
  return (pos == b.pos) && (base == b.base);
}


bool
SnpEvent::operator!=(SnpEvent const & b) const
{
  return !(*this == b);
}


bool
SnpEvent::operator<(SnpEvent const & b) const
{
  return pos < b.pos || (pos == b.pos && base < b.base);
}


double
SnpEventInfo::corrected_support() const
{
  return static_cast<double>(hq_count) + static_cast<double>(lq_count) / 2.0;
}


bool
SnpEventInfo::has_good_support(long const cov) const
{
  int const _depth = hq_count + lq_count;
  bool const is_promising = uniq_pos3 != -1 && hq_count >= 4 && proper_pairs >= 3;
  gyper::Options const & copts = *(Options::const_instance());

  return (copts.no_filter_on_begin_pos || uniq_pos2 != -1)
         &&
         (!copts.filter_on_proper_pairs || (proper_pairs >= 2))
         &&
         (hq_count >= 3)
         &&
         (!copts.filter_on_read_bias ||
          is_promising ||
          (first_in_pairs > 0 && first_in_pairs < _depth))
         &&
         (!copts.filter_on_strand_bias ||
          (is_promising && sequence_reversed > 0 && sequence_reversed < _depth) ||
          (sequence_reversed > 1 && sequence_reversed < (_depth - 1)))
         &&
         ((clipped + 3) <= _depth)
         &&
         ((is_promising && hq_count >= 7) || max_distance > 10)
         &&
         (corrected_support() >= 3.9)
         &&
         (cov <= 0 || (static_cast<double>(_depth) / static_cast<double>(cov) > 0.26));
}


std::string
SnpEventInfo::to_string() const
{
  std::ostringstream ss;

  ss << "hq=" << hq_count
     << ",lq=" << lq_count
     << ",pp=" << proper_pairs
     << ",first=" << first_in_pairs
     << ",reversed=" << sequence_reversed
     << ",clipped=" << clipped
     << ",max_mapq=" << static_cast<long>(max_mapq)
     << ",max_distance=" << static_cast<long>(max_distance)
     << ",uniq_pos1=" << uniq_pos1
     << ",uniq_pos2=" << uniq_pos2
     << ",uniq_pos3=" << uniq_pos3
     << ",phase_num=" << phase.size();

  return ss.str();
}


std::vector<uint8_t>
get_phred_biallelic(uint32_t count, uint32_t anti_count, uint32_t eps)
{
  std::vector<uint8_t> phred(3);

  uint64_t gt00 = count * eps;
  uint64_t gt01 = count + anti_count;
  uint64_t gt11 = anti_count * eps;
  uint64_t const min_gt = std::min({gt00, gt01, gt11});
  gt00 = (gt00 - min_gt) * 3ul;
  gt01 = (gt01 - min_gt) * 3ul;
  gt11 = (gt11 - min_gt) * 3ul;

  uint64_t constexpr MAX = std::numeric_limits<uint8_t>::max();
  phred[0] = std::min(gt00, MAX);
  phred[1] = std::min(gt01, MAX);
  phred[2] = std::min(gt11, MAX);
  return phred;
}


uint32_t
get_log_qual(uint32_t count, uint32_t anti_count, uint32_t eps)
{
  uint32_t const gt00 = count * eps;
  uint32_t const gt01 = count + anti_count;
  uint32_t const gt11 = anti_count * eps;
  uint32_t const gt_alt = std::min(gt01, gt11);

  if (gt00 > gt_alt)
    return gt00 - gt_alt;
  else
    return 0u;
}


uint32_t
get_log_qual_double(double count, double anti_count, double eps)
{
  double const gt00 = count * eps;
  double const gt01 = count + anti_count;
  double const gt11 = anti_count * eps;
  double const gt_alt = std::min(gt01, gt11);

  if (gt00 > gt_alt)
    return static_cast<uint32_t>(gt00 - gt_alt + 0.5);
  else
    return 0u;
}


uint32_t
EventInfo::log_qual(uint32_t eps) const
{
  return get_log_qual(count, anti_count, eps);
}


void
EventInfo::clean_counts()
{
  count = 0;
  anti_count = 0;
  multi_count = 0;
}


bool
apply_indel_event(std::vector<char> & sequence,
                  std::vector<int32_t> & ref_positions,
                  IndelEvent const & indel_event,
                  long const offset,
                  bool const is_debug)
{
  assert(sequence.size() == ref_positions.size());
  long const ref_pos{indel_event.pos - offset};

  if (ref_pos <= 0)
    return false;

  long pos{ref_pos}; // Start the search for the reference position here
  long const event_size{static_cast<long>(indel_event.sequence.size())};
  long const seq_size{static_cast<long>(sequence.size())};
  assert(pos < seq_size);

  if (ref_positions[pos] != ref_pos)
  {
    while ((pos + 1) < seq_size && ref_positions[pos] < ref_pos)
      ++pos;

    while (pos > 0 && ref_positions[pos] > ref_pos)
      --pos;

    // Check again
    if (ref_positions[pos] != ref_pos)
    {
      if (is_debug)
        BOOST_LOG_TRIVIAL(warning) << __HERE__ << " Couldn't find pos for event " << indel_event.to_string();

      return false;
    }
  }

  // Check purity
  {
    long constexpr PURITY_PAD{3};
    long const begin = std::max(0l, pos - PURITY_PAD);
    long const end = std::min(static_cast<long>(ref_positions.size()), pos + PURITY_PAD);

    assert(begin >= 0);
    assert(begin <= end);
    assert(end <= static_cast<long>(ref_positions.size()));

    long prev_ref_pos = ref_positions[begin];

    for (long p{begin + 1}; p < end; ++p)
    {
      if (ref_positions[p] == prev_ref_pos + 1)
      {
        ++prev_ref_pos;
      }
      else
      {
        if (is_debug)
          BOOST_LOG_TRIVIAL(debug) << __HERE__ << " Not pure " << indel_event.to_string();

        return false; // Not pure, skip indel
      }
    }
  }

  if (indel_event.type == 'D')
  {
    if ((pos + event_size) >= static_cast<long>(ref_positions.size()) ||
        ref_positions[pos + event_size] != (ref_pos + event_size))
    {
      if (is_debug)
        BOOST_LOG_TRIVIAL(debug) << __HERE__ << " Cannot apply deletion";

      return false;
    }

    std::vector<char>::iterator end_it = sequence.end();
    std::vector<int32_t>::iterator end_pos_it = ref_positions.end();

    if ((pos + event_size) < seq_size)
    {
      end_it = sequence.begin() + pos + event_size;
      end_pos_it = ref_positions.begin() + pos + event_size;
    }

    sequence.erase(sequence.begin() + pos, end_it);
    ref_positions.erase(ref_positions.begin() + pos, end_pos_it);
  }
  else if (indel_event.type == 'I')
  {
    sequence.insert(sequence.begin() + pos, indel_event.sequence.begin(), indel_event.sequence.end());

    std::vector<int32_t> new_pos(event_size, pos + 1);
    ref_positions.insert(ref_positions.begin() + pos + 1, new_pos.begin(), new_pos.end());
  }
  else
  {
    BOOST_LOG_TRIVIAL(error) << "Unknown type: " << indel_event.type;
    return false;
  }

  assert(sequence.size() == ref_positions.size());
  return true;
}


IndelEvent
make_deletion_event(std::vector<char> const & reference_sequence, long ref_offset, int32_t pos, long count)
{
  IndelEvent new_event(pos, 'D');

  assert(ref_offset >= 0);
  assert(ref_offset + count < static_cast<long>(reference_sequence.size()));

  new_event.sequence = std::vector<char>(reference_sequence.begin() + ref_offset,
                                         reference_sequence.begin() + ref_offset + count);

  return new_event;
}


IndelEvent
make_insertion_event(int32_t pos, std::vector<char> && event_sequence)
{
  IndelEvent new_event(pos, 'I');
  new_event.sequence = std::move(event_sequence);
  return new_event;
}


} // namespace gyper
