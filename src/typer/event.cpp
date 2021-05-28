#include <string>
#include <vector>

#include <boost/functional/hash.hpp> // boost::hash_range
#include <graphtyper/utilities/logging.hpp>


#include <graphtyper/constants.hpp>
#include <graphtyper/typer/event.hpp>
#include <graphtyper/utilities/options.hpp>


namespace
{
/*
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
*/

} // anon namespace


namespace gyper
{

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
  if (qual == 0)
  {
    ++unknown;
    return;
  }

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
    print_log(log_severity::error, __HERE__, " Unexpected base '", base, "'");
    std::exit(1);
  }
  }   // switch ends
}


bool
Event::operator==(Event const & b) const
{
  return (pos == b.pos) &&
         (type == b.type) &&
         (sequence == b.sequence);
}


bool
Event::operator!=(Event const & b) const
{
  return !(*this == b);
}


bool
Event::operator<(Event const & b) const
{
  // order is indel, deletions, snp
  auto const order_a = (type == 'D') + 2 * (type == 'X');
  auto const order_b = (b.type == 'D') + 2 * (b.type == 'X');

  return pos < b.pos ||
         (pos == b.pos && order_a < order_b) ||
         (pos == b.pos && order_a == order_b && sequence < b.sequence);
}


void
EventSupport::clear()
{
  hq_count = 0;
  lq_count = 0;
  proper_pairs = 0;
  first_in_pairs = 0;
  sequence_reversed = 0;
  clipped = 0;
  max_mapq = 0;
  max_distance = 0;
  uniq_pos1 = -1;
  uniq_pos2 = -1;
  uniq_pos3 = -1;
  //phase.clear();
  // do not clear indel only things
}


int
EventSupport::get_raw_support() const
{
  return hq_count + lq_count;
}


double
EventSupport::corrected_support() const
{
  return static_cast<double>(hq_count) + static_cast<double>(lq_count) / 2.0;
}


bool
EventSupport::has_good_support(long cov) const
{
  gyper::Options const & copts = *(Options::const_instance());

  if (cov < 1)
    cov = 1;

  int const raw_support = get_raw_support();

  bool const is_very_promising =
    uniq_pos3 != -1 &&
    (
      (hq_count >= 8 && ((static_cast<double>(raw_support) / static_cast<double>(cov)) >= 0.35)) ||
      (hq_count >= 7 && ((static_cast<double>(raw_support) / static_cast<double>(cov)) >= 0.40))
    ) && (!copts.filter_on_proper_pairs || proper_pairs >= 6);

  bool const is_promising =
    uniq_pos3 != -1 &&
    (
      (hq_count >= 7 && ((static_cast<double>(raw_support) / static_cast<double>(cov)) >= 0.20)) ||
      (hq_count >= 6 && ((static_cast<double>(raw_support) / static_cast<double>(cov)) >= 0.30)) ||
      (hq_count >= 5 && ((static_cast<double>(raw_support) / static_cast<double>(cov)) >= 0.40))
    ) &&
    (!copts.filter_on_proper_pairs || proper_pairs >= 4);

  return (copts.no_filter_on_begin_pos || uniq_pos2 != -1)
         &&
         (!copts.filter_on_proper_pairs || (proper_pairs >= 2))
         &&
         (hq_count >= 3)
         &&
         (!copts.filter_on_read_bias ||
          is_promising ||
          (first_in_pairs > 0 && first_in_pairs < raw_support))
         &&
         (is_very_promising ||
          !copts.filter_on_strand_bias ||
          (is_promising && sequence_reversed > 0 && sequence_reversed < raw_support) ||
          (sequence_reversed > 1 && sequence_reversed < (raw_support - 1)))
         &&
         (clipped <= 1 || (clipped + 5) <= raw_support)
         &&
         (max_distance >= 10 || (is_promising && hq_count >= 10))
         &&
         (corrected_support() >= 3.9)
         &&
         (((static_cast<double>(raw_support) / static_cast<double>(cov)) > 0.26) || is_promising);
}


std::string
EventSupport::to_string() const
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
     << ",phase_num=" << phase.size()
     << ",anti_count=" << anti_count
     << ",span=" << span
     << ",max_log_qual=" << max_log_qual
     << ",max_log_qual_file_i=" << max_log_qual_file_i;

  return ss.str();
}


uint32_t
EventSupport::log_qual(uint32_t eps) const
{
  return get_log_qual(hq_count + lq_count, anti_count, eps);
}


bool
EventSupport::is_good_indel(uint32_t eps) const
{
  long const depth = hq_count + lq_count + anti_count + multi_count;

  if (hq_count <= 6 || sequence_reversed <= 0 || sequence_reversed >= depth || proper_pairs <= 4 ||
      (hq_count < 10 && max_mapq <= 10)
      )
  {
    return false;
  }

  assert(sequence_reversed <= depth);
  long const qual = 3 * get_log_qual(hq_count + lq_count, anti_count, eps);

  if (qual < 50)
    return false;

  double const qd = static_cast<double>(qual) / static_cast<double>(depth);
  return qd >= 3.5;
}


bool
apply_indel_event(std::vector<char> & sequence,
                  std::vector<int32_t> & ref_positions,
                  Event const & indel_event,
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

  if (pos >= seq_size)
    return false;

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
        print_log(log_severity::warning, __HERE__, " Couldn't find pos for event ", indel_event.to_string());

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
          print_log(log_severity::debug, __HERE__, " Not pure ", indel_event.to_string());

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
        print_log(log_severity::debug, __HERE__, " Cannot apply deletion");

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
    print_log(log_severity::error, "Unknown type: ", indel_event.type);
    return false;
  }

  assert(sequence.size() == ref_positions.size());
  return true;
}


Event
make_deletion_event(std::vector<char> const & reference_sequence, long ref_offset, int32_t pos, long count)
{
  Event new_event(pos, 'D');

  assert(ref_offset >= 0);
  assert(ref_offset + count < static_cast<long>(reference_sequence.size()));

  new_event.sequence = std::vector<char>(reference_sequence.begin() + ref_offset,
                                         reference_sequence.begin() + ref_offset + count);

  return new_event;
}


Event
make_insertion_event(int32_t pos, std::vector<char> && event_sequence)
{
  Event new_event(pos, 'I');
  new_event.sequence = std::move(event_sequence);
  return new_event;
}


std::size_t
EventHash::operator()(gyper::Event const & e) const
{
  std::size_t h1 = std::hash<uint32_t>()(e.pos);
  std::size_t h2 = std::hash<char>()(e.type);
  std::size_t h3 = 42 + boost::hash_range(e.sequence.begin(), e.sequence.end());
  return h1 ^ (h2 << 1) ^ (h3 + 0x9e3779b9);
}


} // namespace gyper
