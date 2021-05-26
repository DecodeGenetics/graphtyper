#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>

#include <graphtyper/utilities/logging.hpp>

#include <graphtyper/typer/segment.hpp>
#include <graphtyper/utilities/graph_help_functions.hpp>


namespace
{

std::vector<uint8_t>
get_segment_phred_from_log2_scores(std::vector<uint32_t> const & log_score)
{
  using namespace gyper;
  assert(log_score.size() > 0);

  // Check if all phred scores are zero by finding the first non-zero phred score
  auto find_it = std::find_if(log_score.begin() + 1,
                              log_score.end(),
                              [&](uint32_t const val)
    {
      return val != log_score[0];
    });

  if (find_it == log_score.end())
  {
    print_log(log_severity::warning, "[graphtyper::segment] All PHRED scores zero");
    return std::vector<uint8_t>(log_score.size(), 0u); // If no non-zero phred score is found, the phred scores of the haplotype should be all zero as well
  }


  std::vector<uint8_t> segment_phred(log_score.size(), 255);

  // First find out what the maximum log score is
  uint32_t const max_log_score = *std::max_element(log_score.begin(), log_score.end());

  for (std::size_t i = 0; i < log_score.size(); ++i)
  {
    assert (log_score[i] <= max_log_score);
    double constexpr LOG10_HALF_times_10 = 3.0102999566398119521373889472449302676818988146210854131;
    uint64_t const phred_score = std::llround(static_cast<double>(max_log_score - log_score[i]) * LOG10_HALF_times_10);

    if (phred_score < 255u)
      segment_phred[i] = static_cast<uint8_t>(phred_score);
  }

  return segment_phred;
}


} // anon namespace


namespace gyper
{

Segment::Segment(){}


Segment::Segment(uint32_t _id, std::size_t _ref_size, std::vector<std::string> const & _allele_names)
  : id(_id)
  , ref_size(_ref_size)
  , allele_names(_allele_names)
  , var_type("XG")
{
  assert(allele_names.size() > 1);

  if (allele_names[0].substr(0, 4) == std::string("HLA-") ||
      allele_names[0].substr(0, 4) == std::string("TAP1") ||
      allele_names[0].substr(0, 4) == std::string("TAP2") ||
      allele_names[0].substr(0, 4) == std::string("MICA") ||
      allele_names[0].substr(0, 4) == std::string("MICB") ||
      allele_names[0].substr(0, 3) == std::string("HFE")
      )
  {
    var_type = "H";
  }
}


void
Segment::clear()
{
  allele_names.clear();
  allele_names.shrink_to_fit();
  segment_calls.clear();
  segment_calls.shrink_to_fit();
  var_type.clear();
}


void
Segment::insert_score(std::vector<uint32_t> const & score)
{
  SegmentCall new_segment_call;
  new_segment_call.phred = get_segment_phred_from_log2_scores(score);
  assert(new_segment_call.phred.size() == score.size());
  assert(new_segment_call.phred.size() > 0);

  auto find_it = std::find(new_segment_call.phred.begin(), new_segment_call.phred.end(), 0);

  if (find_it == new_segment_call.phred.end())
  {
    print_log(log_severity::error, "[graphtyper::segment] No zero PHRED score found in segment!");
    std::exit(1);
  }

  new_segment_call.call = to_pair(std::distance(new_segment_call.phred.begin(), find_it));
  segment_calls.push_back(std::move(new_segment_call));
}


void
Segment::insert_scores(std::vector<std::vector<uint32_t> > const & scores)
{
  segment_calls.reserve(scores.size());

  for (unsigned i = 0; i < scores.size(); ++i)
    insert_score(scores[i]);

  assert(segment_calls.size() == scores.size());
}


std::string
Segment::get_ref_string() const
{
  assert(allele_names.size() > 1);
  std::stringstream hap_ss;
  hap_ss << '<' << allele_names[0] << '>';
  return hap_ss.str();
}


std::string
Segment::get_alt_string() const
{
  assert(allele_names.size() > 1);
  std::stringstream hap_ss;
  hap_ss << '<' << *(allele_names.begin() + 1) << '>';

  for (auto hap_it = allele_names.cbegin() + 2; hap_it != allele_names.cend(); ++hap_it)
  {
    hap_ss << ",<" << *hap_it << '>';
  }

  return hap_ss.str();
}


std::vector<Segment>
Segment::get_biallelic_segments() const
{
  std::vector<Segment> new_segments;
  assert(this->allele_names.size() >= 2);

  for (uint32_t alt_id = 0; alt_id < this->allele_names.size(); ++alt_id)
  {
    Segment new_segment;
    new_segment.id = this->id;
    new_segment.ref_size = this->ref_size;
    new_segment.segment_name = allele_names[alt_id];
    new_segment.var_type = this->var_type;
    new_segment.allele_names.push_back(std::string("!")+this->allele_names[alt_id]);
    new_segment.allele_names.push_back(this->allele_names[alt_id]);
    new_segment.extra_id = alt_id;

    {
      SegmentCall new_segment_call;
      new_segment_call.phred = std::vector<uint8_t>(3, 255u);
      new_segment_call.call = {0xFFFFu, 0xFFFFu};
      new_segment.segment_calls.resize(this->segment_calls.size(), new_segment_call);
    }

    for (unsigned s = 0; s < new_segment.segment_calls.size(); ++s)
    {
      assert (s < segment_calls.size());
      auto const & segment_call = this->segment_calls[s];
      auto & new_calls = new_segment.segment_calls[s];
      uint32_t i = 0;

      for (uint32_t y = 0; y < allele_names.size(); ++y)
      {
        for (uint32_t x = 0; x <= y; ++x, ++i)
        {
          assert (i < segment_call.phred.size());

          if (segment_call.phred[i] == 255u)
            continue;

          // Check for homozygous alt, then heterozygous alt
          if (x == alt_id && y == alt_id)
          {
            if (segment_call.phred[i] == 0 && new_calls.call.first == 0xFFFFu)
              new_calls.call = {1, 1};

            new_calls.phred[2] = std::min(new_calls.phred[2], segment_call.phred[i]); // Homozygous var
          }
          else if (x == alt_id || y == alt_id)
          {
            if (segment_call.phred[i] == 0 && new_calls.call.first == 0xFFFFu)
              new_calls.call = {0, 1};

            new_calls.phred[1] = std::min(new_calls.phred[1], segment_call.phred[i]); // Heterozygous
          }
          else
          {
            if (segment_call.phred[i] == 0 && new_calls.call.first == 0xFFFFu)
              new_calls.call = {0, 0};

            new_calls.phred[0] = std::min(new_calls.phred[0], segment_call.phred[i]); // Homozygous ref
          }
        }
      }

      if (new_calls.call.first == 0xFFFFu)
      {
        print_log(log_severity::error, "[graphtyper::segment] No phred score zero found in a segment!");
        std::exit(1);
      }
    }

    new_segments.push_back(std::move(new_segment));
  }

  return new_segments;
}


} // namespace gyper
