#include <algorithm> // std::sort
#include <iostream> // std::cout, std::cerr, std::endl

#include <graphtyper/utilities/hts_reader.hpp>

#include <boost/algorithm/string.hpp>
#include <boost/log/trivial.hpp>

#include <htslib/hfile.h>
#include <htslib/hts.h>


namespace gyper
{

HtsReader::HtsReader(HtsStore & _store)
  : store(_store)
{}


void
HtsReader::open(std::string const & path)
{
  fp = hts_open(path.c_str(), "r");

  if (!fp)
  {
    std::cerr << "ERROR: Could not open BAM file  " << path << std::endl;
    std::exit(1);
  }

  fp->bam_header = sam_hdr_read(fp);

  // Read sample from header
  if (!Options::instance()->get_sample_names_from_filename)
  {
    std::string const header_text(fp->bam_header->text, fp->bam_header->l_text);
    std::vector<std::string> header_lines;

    // Split the header text into lines
    boost::split(header_lines, header_text, boost::is_any_of("\n"));

    for (auto & line_it : header_lines)
    {
      if (boost::starts_with(line_it, "@RG"))
      {
        std::size_t const pos_id = line_it.find("\tID:");
        std::size_t const pos_samp = line_it.rfind("\tSM:");

        if (pos_samp == std::string::npos || pos_id == std::string::npos)
        {
          std::cerr << "[graphtyper::utilities::hts_reader] ERROR: Could not parse RG and sample from header line:"
                    << line_it << std::endl;
          std::exit(1);
        }

        std::size_t pos_id_ends = line_it.find("\t", pos_id + 1);

        // Check if this is the last field
        if (pos_id_ends == std::string::npos)
          pos_id_ends = line_it.size();

        std::size_t pos_samp_ends = line_it.find("\t", pos_samp + 1);

        // Check if this is the last field
        if (pos_samp_ends == std::string::npos)
          pos_samp_ends = line_it.size();

        std::string new_id = line_it.substr(pos_id + 4, pos_id_ends - pos_id - 4);
        std::string new_sample = line_it.substr(pos_samp + 4, pos_samp_ends - pos_samp - 4);

        BOOST_LOG_TRIVIAL(debug) << "[graphtyper::utilities::hts_reader] Added RG: '"
                                 << new_id << "' => '" << new_sample << "'";

        rg2index[new_id] = rg2sample_i.size();
        auto find_it = std::find(samples.begin(), samples.end(), new_sample);

        // check if this is a new sample
        if (find_it == samples.end())
        {
          rg2sample_i.push_back(samples.size());
          samples.push_back(new_sample);
        }
        else
        {
          rg2sample_i.push_back(std::distance(samples.begin(), find_it));
        }
      }
    }
  }

  // Parse sample name from path name if we need to
  if (samples.size() == 0)
  {
    std::string sample = path.substr(path.rfind('/') + 1, path.rfind('.'));

    if (std::count(sample.begin(), sample.end(), '.') > 0)
      sample = sample.substr(0, sample.find('.'));

    samples.push_back(std::move(sample));

    // No read group
    //rg2sample_i.push_back(0);
    //rg2index[""] = 0;
  }

  rec = store.get();
  ret = sam_read1(fp, fp->bam_header, rec);
}


void
HtsReader::close()
{
  if (fp)
  {
    hts_close(fp);
    fp = nullptr;
  }
}


int
HtsReader::set_reference(std::string const & reference_path)
{
  int ret2 = hts_set_fai_filename(fp, reference_path.c_str());

  if (ret2 < 0)
  {
    BOOST_LOG_TRIVIAL(error) << "[graphtyper::utilities::hts_reader] ERROR: "
                             << "Could not open reference FASTA file with filename " << reference_path;
    std::exit(1);
  }

  return ret2;
}


void
HtsReader::set_sample_index_offset(int const new_sample_index_offset)
{
  sample_index_offset = new_sample_index_offset;
}


void
HtsReader::set_rg_index_offset(int new_rg_index_offset)
{
  rg_index_offset = new_rg_index_offset;
}


bam1_t *
HtsReader::get_next_read(bam1_t * old_record)
{
  assert(old_record);

  // If we have some records ready in the vector, return those first
  if (records.size() > 0)
  {
    bam1_t * record = *(records.end() - 1);
    records.pop_back();
    store.push(old_record);   // We do not need this memory now
    return record;
  }

  // We do not have any records ready and we have reached the end of the file, stop now
  if (ret < 0)
  {
    // Deallocate memory when we reach the end of the file
    if (rec)
      store.push(rec);

    store.push(old_record);
    return nullptr;
  }

  // Read until a new position is found
  auto const pos = rec->core.pos;
  records.push_back(rec);
  rec = old_record;
  ret = sam_read1(fp, fp->bam_header, rec);

  // Read while the records have the same position
  while (ret >= 0 && rec->core.pos == pos)
  {
    assert(rec);
    records.push_back(rec);
    rec = store.get();
    assert(rec);
    ret = sam_read1(fp, fp->bam_header, rec);
  }

  std::sort(records.begin(), records.end(), gt_pos_seq_same_pos);

  bam1_t * record = *(records.end() - 1);
  records.pop_back();
  return record;
}


bam1_t *
HtsReader::get_next_read()
{
  // If we have some records ready in the vector, return those first
  if (records.size() > 0)
  {
    bam1_t * record = *(records.end() - 1);
    records.pop_back();
    return record;
  }

  // We do not have any records ready and we have reached the end of the file, stop now
  if (ret < 0)
  {
    // Deallocate memory when we reach the end of the file
    if (rec)
      store.push(rec);

    return nullptr;
  }

  // Read until a new position is found
  auto const pos = rec->core.pos;
  records.push_back(rec);
  rec = store.get();
  assert(rec);
  ret = sam_read1(fp, fp->bam_header, rec);

  // Read while the records have the same position
  while (ret >= 0 && rec->core.pos == pos)
  {
    assert(rec);
    records.push_back(rec);
    rec = store.get();
    assert(rec);
    ret = sam_read1(fp, fp->bam_header, rec);
  }

  std::sort(records.begin(), records.end(), gt_pos_seq_same_pos);

  bam1_t * record = *(records.end() - 1);
  records.pop_back();
  return record;
}


bam1_t *
HtsReader::get_next_pos()
{
  if (ret < 0)
  {
    // Deallocate memory when we reach the end of the file
    bam_destroy1(rec);
    return nullptr;
  }

  bam1_t * record = rec;
  rec = store.get();
  ret = sam_read1(fp, fp->bam_header, rec);
  return record;
}


void
HtsReader::get_sample_and_rg_index(long & sample_i, long & rg_i, bam1_t * rec) const
{
  //assert(rg2index.size() > 0);

  if (rg2sample_i.size() <= 1)
  {
    rg_i = rg_index_offset;
    sample_i = sample_index_offset;
  }
  else
  {
    // Check read group tag from SAM
    uint8_t * rg_tag = bam_aux_get(rec, "RG");

    if (!rg_tag)
    {
      std::cerr << "[graphtyper::utilities::hts_reader] ERROR: Unable to find RG tag in read." << std::endl;
      std::exit(1);
    }

    std::string read_group(reinterpret_cast<char *>(rg_tag + 1)); // Skip 'Z'
    auto find_rg_it = rg2index.find(read_group);

    if (find_rg_it == rg2index.end())
    {
      std::cerr << "[graphtyper::utilities::hts_reader] ERROR: Unable to find read group. " << read_group << std::endl;
      std::exit(1);
    }

    rg_i = find_rg_it->second + rg_index_offset;
    assert(find_rg_it->second < static_cast<long>(rg2sample_i.size()));
    sample_i = rg2sample_i[find_rg_it->second] + sample_index_offset;
  }
}


long
HtsReader::get_num_rg() const
{
  assert(rg2index.size() == rg2sample_i.size());
  return std::max(1l, static_cast<long>(rg2sample_i.size()));
}


} // namespace gyper
