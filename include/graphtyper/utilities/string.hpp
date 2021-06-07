#pragma once

#if __has_include(<charconv>)
#  include <charconv>
#endif
#include <string>
#include <string_view>
#include <vector>

#include <graphtyper/utilities/logging.hpp>

namespace gyper
{
//----------------------------------- starts_with / ends_with ----------------------------
// TODO: this can be implemented much more generically and nicer with c++ ranges / c++20

inline bool starts_with(std::string_view input, std::string_view prefix)
{
  auto in_it = begin(input);
  auto in_end = end(input);

  auto pre_it = begin(prefix);
  auto pre_end = end(prefix);

  if (pre_end - pre_it > in_end - in_it)
    return false;

  while (in_it != in_end && pre_it != pre_end)
  {
    if (*in_it != *pre_it)
      return false;

    ++in_it;
    ++pre_it;
  }

  return true;
}

inline bool ends_with(std::string_view input, std::string_view prefix)
{
  auto in_it = rbegin(input);
  auto in_end = rend(input);

  auto pre_it = rbegin(prefix);
  auto pre_end = rend(prefix);

  if (pre_end - pre_it > in_end - in_it)
    return false;

  while (in_it != in_end && pre_it != pre_end)
  {
    if (*in_it != *pre_it)
      return false;

    ++in_it;
    ++pre_it;
  }

  return true;
}

//----------------------------------- split_on_delim ------------------------------------
// TODO replace with views in C++20

// cannot return views into temporary of string
inline std::vector<std::string_view> split_on_delim(std::string && input, char const delim) = delete;

inline std::vector<std::string_view> split_on_delim(std::string_view input, char const delim)
{
  std::vector<std::string_view> ret;

  for (size_t i = 0, s = 1; i <= input.size(); ++i, ++s)
  {
    if (i == input.size() || input[i] == delim)
    {
      //                begin_ptr                , size
      ret.emplace_back(input.data() + i - (s - 1), s - 1);
      s = 0;
    }
  }

  return ret;
}

//----------------------------------- stoll that works with string_view ------------------------
// TODO remove special case when removing support for GCC7

#if __has_include(<charconv>)
inline int64_t stoi64(std::string_view str)
{
  int64_t ret = 0;
  auto tmp = std::from_chars(str.data(), str.data() + str.size(), ret);

  if (tmp.ec != std::errc())
  {
    print_log(log_severity::error, "Converting \"", str, "\" of size ", str.size(), " to int failed.");
    std::exit(1);
  }

  return ret;
}
#else // GCC7
inline int64_t stoi64(std::string_view str)
{
  std::string tmp{str}; // needless copy to get null-terminated string
  return std::stoll(tmp);
}
#endif

// Returns an iterator to the nth occurence of an element (char). Expects n>=1.
template <typename Tit>
inline Tit find_nth_occurence(Tit begin, Tit end, char const element, int n)
{
  assert(n >= 1);
  int count{0};
  return std::find_if(begin, end, [&count, element, n](char x) { return x == element && ++count == n; });
}

} // namespace gyper
