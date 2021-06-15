#pragma once

#include <cstdio> // std::exit
#include <cstring>
#include <iostream> // std::cout
#include <memory>
#include <sstream> // std::stringstream
#include <string>  // std::string

#include <graphtyper/constants.hpp>

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

  std::string filename{};
  std::string filemode{};
  long n_threads{1};
  long static constexpr MAX_CACHE_SIZE{50000ll};
};

template <class T>
inline BGZF_stream & BGZF_stream::operator<<(T const & x)
{
  // Add to stringstream
  ss << x;
  return *this;
}

inline void BGZF_stream::check_cache()
{
  if (static_cast<long>(ss.tellp()) > MAX_CACHE_SIZE)
    flush();
}

inline void BGZF_stream::flush()
{
  // Write stringstream to BGZF file
  if (fp == nullptr)
  {
    std::cout << ss.str(); // Write uncompressed to stdout
  }
  else
  {
    std::string str = ss.str();
    int written_length = bgzf_write(fp, str.data(), str.size());

    if (written_length < 0)
    {
      std::cerr << "[bgzf_stream] ERROR: Writing to BGZF file failed. "
                << "Exit code: " << written_length << " . No space left on device?" << std::endl;
      std::exit(1);
    }
    else if (written_length != static_cast<long>(str.size()))
    {
      std::cerr << "[bgzf] WARNING: Mismatch between size written and expected: " << written_length
                << " != " << str.size();
    }
  }

  // Clear stringstream
  ss.str(std::string());
  ss.clear();
}

inline int BGZF_stream::write(void const * data, std::size_t length)
{
  if (fp == nullptr)
  {
    std::cout << std::string(reinterpret_cast<const char *>(data), length);
    return length;
  }

  int const written_length = bgzf_write(fp, data, length);

  if (written_length < 0)
  {
    std::cerr << "[bgzf_stream] ERROR: Writing to BGZF file failed. "
              << "Exit code: " << written_length << " . No space left on device?" << std::endl;
    std::exit(1);
  }
  else if (written_length != static_cast<long>(length))
  {
    std::cerr << "[bgzf] WARNING: Mismatch between size written and expected: " << written_length << " != " << length;
  }

  return written_length != static_cast<long>(length);
}

inline int BGZF_stream::write(std::string const & str)
{
  if (fp == nullptr)
  {
    std::cout << str;
    return str.size();
  }

  std::size_t const length = str.size();
  int const written_length = bgzf_write(fp, str.data(), length);

  if (written_length < 0)
  {
    std::cerr << "[bgzf_stream] ERROR: Writing to BGZF file failed. "
              << "Exit code: " << written_length << " . No space left on device?" << std::endl;
    std::exit(1);
  }
  else if (written_length != static_cast<long>(length))
  {
    std::cerr << "[bgzf] WARNING: Mismatch between size written and expected: " << written_length << " != " << length;
  }

  return written_length != static_cast<long>(length);
}

inline void BGZF_stream::open(std::string const & _filename, std::string const & _filemode, long const _n_threads)
{
  if (fp != nullptr)
    close();

  if (_filename.size() > 0 && _filename != "-")
  {
    fp = bgzf_open(_filename.c_str(), _filemode.c_str());

    if (_n_threads > 1)
      bgzf_mt(fp, _n_threads, 256);
  }

  // Options are stored for uncompressed sample names writing
  this->filename = _filename;
  this->filemode = _filemode;
  this->n_threads = _n_threads;
}

inline bool BGZF_stream::is_open() const
{
  return fp != nullptr;
}

inline void BGZF_stream::close()
{
  flush();

  if (fp != nullptr)
  {
    int ret = bgzf_close(fp);

    if (ret != 0)
    {
      std::cerr << " " << __HERE__ << "[bgzf_stream] ERROR: Failed closing bgzf file " << filename << "\n";
      std::exit(1);
    }

    fp = nullptr;
  }
}

inline BGZF_stream::~BGZF_stream()
{
  close();
}

} // namespace gyper
