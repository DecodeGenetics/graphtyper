#pragma once

#include <cassert>
#include <cstdio> // std::exit
#include <cstring>
#include <iostream> // std::cout
#include <memory>
#include <sstream> // std::stringstream
#include <string>  // std::string

#include <graphtyper/constants.hpp>

#include <popvcf/encode.hpp>

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
  std::string buffer_in;
  std::string buffer_out;

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

  popvcf::EncodeData ed;
  std::string filename{};
  std::string filemode{};
  long n_threads{1};
  long static constexpr MAX_CACHE_SIZE{200000ll};
};

inline void BGZF_stream::check_cache()
{
  buffer_in.append(ss.str());
  ss.str(std::string()); // clear stringstream
  ss.clear();

  if (buffer_in.size() > MAX_CACHE_SIZE)
    flush();
}

inline void BGZF_stream::flush()
{
  buffer_in.append(ss.str());
  ss.str(std::string()); // clear stringstream
  ss.clear();

  if (Options::const_instance()->is_on_final_output && Options::const_instance()->encoding == 'p')
  {
    assert(buffer_in.size() >= ed.i);
    ed.bytes_read = buffer_in.size();
    popvcf::encode_buffer(buffer_out, buffer_in, ed);

    if (fp == nullptr)
    {
      std::cout << buffer_out; // Write uncompressed to stdout
    }
    else if (buffer_out.size() > 0)
    {
      int written_length = bgzf_write(fp, buffer_out.data(), buffer_out.size());

      if (written_length < 0)
      {
        std::cerr << "[bgzf_stream] ERROR: Writing to BGZF file failed. "
                  << "Exit code: " << written_length << " . No space left on device?" << std::endl;
        std::exit(1);
      }
      else if (written_length != static_cast<long>(buffer_out.size()))
      {
        std::cerr << "[bgzf] WARNING: Mismatch between size written and expected: " << written_length
                  << " != " << buffer_out.size();
      }
    }

    buffer_in.resize(ed.i);
    buffer_out.resize(0);
  }
  else
  {
    if (buffer_in.empty())
      return;

    // Write stringstream to BGZF file
    if (fp == nullptr)
    {
      std::cout << buffer_in; // Write uncompressed to stdout
    }
    else
    {
      int written_length = bgzf_write(fp, buffer_in.data(), buffer_in.size());

      if (written_length < 0)
      {
        std::cerr << "[bgzf_stream] ERROR: Writing to BGZF file failed. "
                  << "Exit code: " << written_length << " . No space left on device?" << std::endl;
        std::exit(1);
      }
      else if (written_length != static_cast<long>(buffer_in.size()))
      {
        std::cerr << "[bgzf] WARNING: Mismatch between size written and expected: " << written_length
                  << " != " << buffer_in.size();
      }
    }

    buffer_in.resize(0);
  }
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
  assert(buffer_in.empty());
  assert(buffer_out.empty());
  assert(ss.str().empty());

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
