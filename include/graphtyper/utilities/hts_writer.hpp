#pragma once

#include <string>
#include <htslib/sam.h>

namespace gyper
{

class HtsParallelReader;


class HtsWriter
{

private:
  htsFile * fp = nullptr; // htslib file pointer

public:
  HtsWriter() = default;

  // opens a hts file at path in a format
  void open(std::string const & path);

  // closes the hts file
  void close();

  // copies a header from a HtsParallelReader
  void copy_header(HtsParallelReader const & hts_preader);

  // write a hts record
  void write(bam1_t * hts_rec);

};

} // namespace gyper
