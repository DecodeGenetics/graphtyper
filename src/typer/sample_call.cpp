#include <algorithm>
#include <iostream>

#include <graphtyper/typer/sample_call.hpp>


namespace gyper
{

SampleCall::SampleCall() noexcept
   : ambiguous_depth(0)
{}

SampleCall::SampleCall(std::vector<uint8_t> const & _phred,
                       std::vector<uint16_t> const & _coverage,
                       uint8_t const _ambiguous_depth) noexcept
  : phred(_phred)
  , coverage(_coverage)
  , ambiguous_depth(_ambiguous_depth)
{}


std::pair<uint16_t, uint16_t>
SampleCall::get_gt_call() const
{
  std::pair<uint16_t, uint16_t> gt_call = {0xFFFFu, 0xFFFFu};

  for (uint16_t y = 0; y < coverage.size(); ++y)
  {
    for (uint16_t x = 0; x <= y; ++x)
    {
      if (phred[to_index(x, y)] == 0)
      {
        gt_call.first = x;
        gt_call.second = y;
        y = coverage.size();   // To stop the outer loop
        break;
      }
    }
  }

  if (gt_call.first == 0xFFFFu)
  {
    std::cerr << "ERROR: No phred score zero found on variant" << std::endl;
    std::exit(1);
  }

  assert(gt_call.first < coverage.size());
  assert(gt_call.second < coverage.size());
  return gt_call;
}


uint8_t
SampleCall::get_gq() const
{
  bool seen_zero = false;
  uint8_t next_lowest_phred = 255;

  for (auto const p : phred)
  {
    if (p == 0)
    {
      if (!seen_zero)
        seen_zero = true;
      else
        return 0;
    }
    else
    {
      next_lowest_phred = std::min(next_lowest_phred, p);
    }
  }

  return next_lowest_phred;
}


} // namespce gyper
