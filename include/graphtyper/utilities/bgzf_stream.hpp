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
  std::ostringstream ss;


public:
  BGZF_stream() = default;
  BGZF_stream(std::string const & filename, std::string const & filemode);
  ~BGZF_stream();

  template <class T>
  BGZF_stream & operator<<(T const & x);
  void flush();
  void open(std::string const & filename, std::string const & filemode);
  bool is_open() const;
  void close(); // Close BGZF file

  long MAX_CACHE_SIZE = 10000000ll;
};


inline
BGZF_stream::BGZF_stream(std::string const & filename, std::string const & filemode)
{
  open(filename, filemode);
}


template <class T>
inline
BGZF_stream &
BGZF_stream::operator<<(T const & x)
{
  // Add to stringstream
  ss << x;

  // Check if we should flush
  if (static_cast<long>(ss.tellp()) > this->MAX_CACHE_SIZE)
    flush();

  return *this;
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
    std::string const & str = ss.str();
    int ret = bgzf_write(fp, str.data(), str.size());

    if (ret < 0)
    {
      std::cerr << "ERROR: Writing to BGZF file failed." << std::endl;
      std::exit(1);
    }
  }

  // Clear stringstream
  ss.str("");
}


inline
void
BGZF_stream::open(std::string const & filename, std::string const & filemode)
{
  if (fp)
    close();

  if (filename.size() > 0 && filename != "-")
  {
    fp = bgzf_open(filename.c_str(), filemode.c_str());
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
