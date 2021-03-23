#pragma once

#include <cstdio> // std::exit
#include <cstring>
#include <iostream> // std::cout
#include <memory>
#include <sstream> // std::stringstream
#include <string> // std::string

#define INCLUDE_SEQAN_STREAM_IOSTREAM_BGZF_H_
#include "bgzf.h" // part of htslib


namespace gyper
{


class BGZF_stream
{
private:
  BGZF * fp = nullptr;

public:
  std::ostringstream ss;

  BGZF_stream() = default;
  ~BGZF_stream(); // custom
  BGZF_stream(BGZF_stream const &) = delete;
  BGZF_stream(BGZF_stream &&) = delete;
  BGZF_stream & operator=(BGZF_stream const &) = delete;
  BGZF_stream & operator=(BGZF_stream &&) = delete;

  template <class T>
  BGZF_stream & operator<<(T const & x);
  void check_cache();
  void flush();
  int write(void const * data, std::size_t length); // call bgzf_write directly
  int write(std::string const & str);
  void open(std::string const & filename, std::string const & filemode, long const n_threads);
  bool is_open() const;
  void close(); // Close BGZF file

  long MAX_CACHE_SIZE{1000000ll};
};


template <class T>
inline
BGZF_stream &
BGZF_stream::operator<<(T const & x)
{
  // Add to stringstream
  ss << x;
  return *this;
}


inline
void
BGZF_stream::check_cache()
{
  if (static_cast<long>(ss.tellp()) > this->MAX_CACHE_SIZE)
    flush();
}


inline
void
BGZF_stream::flush()
{
  // Write stringstream to BGZF file
  if (!fp)
  {
    std::cout << ss.str(); // Write uncompressed to stdout
  }
  else
  {
    std::string str = ss.str();
    int ret = bgzf_write(fp, str.data(), str.size());

    if (ret < 0)
    {
      std::cerr << "[bgzf_stream] ERROR: Writing to BGZF file failed. No free space on device?" << std::endl;
      std::exit(1);
    }
  }

  // Clear stringstream
  ss.str(std::string());
  ss.clear();
}


inline
int
BGZF_stream::write(void const * data, std::size_t length)
{
  if (!fp)
  {
    std::cout << std::string(reinterpret_cast<const char *>(data), length);
    return length;
  }
  else
  {
    return bgzf_write(fp, data, length);
  }
}


inline
int
BGZF_stream::write(std::string const & str)
{
  if (!fp)
  {
    std::cout << str;
    return str.size();
  }
  else
  {
    return bgzf_write(fp, str.data(), str.size());
  }
}


inline
void
BGZF_stream::open(std::string const & filename, std::string const & filemode, long const n_threads)
{
  if (fp)
    close();

  if (filename.size() > 0 && filename != "-")
  {
    fp = bgzf_open(filename.c_str(), filemode.c_str());

    if (n_threads > 1)
      bgzf_mt(fp, n_threads, 256);
  }
}


inline bool
BGZF_stream::is_open() const
{
  return fp;
}


inline
void
BGZF_stream::close()
{
  flush();

  if (fp)
  {
    bgzf_close(fp);
    fp = nullptr;
  }
}


inline
BGZF_stream::~BGZF_stream()
{
  close();
}


} // namespace gyper
