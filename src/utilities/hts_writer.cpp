#include <iostream>
#include <string>

#include <graphtyper/utilities/hts_parallel_reader.hpp>
#include <graphtyper/utilities/hts_writer.hpp>

#include <htslib/sam.h>

namespace gyper
{
void HtsWriter::open(std::string const & path)
{
  fp = hts_open(path.c_str(), "wb");

  if (!fp)
  {
    std::cerr << "[graphtyper::hts_writer] ERROR: Could not open HTS file  " << path << std::endl;
    std::exit(1);
  }
}

void HtsWriter::close()
{
  if (fp)
  {
    hts_close(fp);
    fp = nullptr;
  }
}

void HtsWriter::copy_header(HtsParallelReader const & hts_preader)
{
  if (!fp)
  {
    std::cerr << "[graphtyper::hts_writer] ERROR: Cannot write header before opening the file" << std::endl;
    std::exit(1);
  }

  // Get new header
  fp->bam_header = hts_preader.get_header();
  int ret = sam_hdr_write(fp, fp->bam_header);

  if (ret < 0)
  {
    std::cerr << "[graphtyper::hts_writer] ERROR: Could not write header (ret = " << ret << ")" << std::endl;
    std::exit(1);
  }
}

void HtsWriter::write(bam1_t * hts_rec)
{
  if (!fp)
  {
    std::cerr << "[graphtyper::hts_writer] ERROR: Cannot write record before opening the file" << std::endl;
    std::exit(1);
  }

  int ret = sam_write1(fp, fp->bam_header, hts_rec);

  if (ret < 0)
  {
    std::cerr << "[graphtyper::hts_writer] ERROR: Could not write record (ret = " << ret << ")" << std::endl;
    std::exit(1);
  }
}

} // namespace gyper
