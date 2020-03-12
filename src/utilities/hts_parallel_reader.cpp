#include <string> // std::string
#include <unordered_map> // std::unordered_map
#include <vector> // std::vector

#include <boost/log/trivial.hpp>

#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <htslib/sam.h>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/hts_io.h>

#include <graphtyper/graph/haplotype_calls.hpp>
#include <graphtyper/graph/reference_depth.hpp>
#include <graphtyper/typer/alignment.hpp>
#include <graphtyper/typer/genotype_paths.hpp>
#include <graphtyper/typer/variant_map.hpp>
#include <graphtyper/typer/vcf.hpp>
#include <graphtyper/typer/vcf_writer.hpp>
#include <graphtyper/utilities/hash_seqan.hpp>
#include <graphtyper/utilities/hts_parallel_reader.hpp>
#include <graphtyper/utilities/hts_store.hpp>
#include <graphtyper/utilities/hts_writer.hpp>
#include <graphtyper/utilities/io.hpp>
#include <graphtyper/utilities/options.hpp>

namespace
{

using TMapGPaths = std::unordered_map<std::string, std::pair<gyper::GenotypePaths, gyper::GenotypePaths> >;

#ifndef NDEBUG
void
check_if_maps_are_empty(std::vector<TMapGPaths> const & maps)
{
  for (long i = 0; i < static_cast<long>(maps.size()); ++i)
  {
    auto const & map = maps[i];

    if (map.size() > 0)
    {
      BOOST_LOG_TRIVIAL(warning)
        << "[graphtyper::hts_parallel_reader] Found a non-empty read map for with size "
        << map.size()
        << "! This likely means these reads have the BAM_FPAIRED flag set but have no mate read.";

      for (auto it = map.begin(); it != map.end(); ++it)
        BOOST_LOG_TRIVIAL(debug) << "Leftover read name: " << it->first;
    }
  }
}


#endif // NDEBUG

} // anon namespace


namespace gyper
{

HtsParallelReader::~HtsParallelReader()
{
  close();
}


void
HtsParallelReader::open(std::vector<std::string> const & hts_file_paths,
                        std::string const & reference,
                        std::string const & region)
{
  for (auto const & bam : hts_file_paths)
  {
    HtsReader f(store);
    f.open(bam, region);

    if (!reference.empty())
      f.set_reference(reference);

    f.set_sample_index_offset(samples.size());
    std::copy(f.samples.begin(), f.samples.end(), std::back_inserter(samples));
    f.set_rg_index_offset(num_rg);
    num_rg += f.get_num_rg();
    hts_files.push_back(std::move(f));
  }

  // Output
  long const n = hts_files.size();

  // Read the first record of each hts file
  for (long i = 0; i < n; ++i)
  {
    bam1_t * hts_rec = hts_files[i].get_next_read();

    if (hts_rec)
    {
      heap.push_back(HtsRecord(hts_rec, i));
      std::push_heap(heap.begin(), heap.end(), Cmp_gt_pair_bam1_t_fun);
    }
  }
}


void
HtsParallelReader::close()
{
  for (auto & hts_f : hts_files)
    hts_f.close();
}


bool
HtsParallelReader::read_record(HtsRecord & hts_record)
{
  if (heap.empty())
    return false;

  auto const i = heap[0].file_index;
  bam1_t * old_record = hts_record.record;
  hts_record.record = heap[0].record;
  hts_record.file_index = i;
  bam1_t * hts_rec;

  // reuse an old bam1_t record if it is available
  if (old_record)
    hts_rec = hts_files[i].get_next_read(old_record);
  else
    hts_rec = hts_files[i].get_next_read();

  std::pop_heap(heap.begin(), heap.end(), Cmp_gt_pair_bam1_t_fun);
  assert(hts_record.record == heap.back().record);

  if (hts_rec)
  {
    (heap.end() - 1)->record = hts_rec;
    assert((heap.end() - 1)->file_index == i);
    std::push_heap(heap.begin(), heap.end(), Cmp_gt_pair_bam1_t_fun);
  }
  else
  {
    (heap.end() - 1)->record = nullptr;
    heap.pop_back();
  }

  return true;
}


void
HtsParallelReader::move_record(HtsRecord & to, HtsRecord & from)
{
  if (to.record)
    store.push(to.record);

  to.record = from.record;
  to.file_index = from.file_index;
  from.record = nullptr;
}


std::vector<std::string> const &
HtsParallelReader::get_samples() const
{
  return samples;
}


void
HtsParallelReader::get_sample_and_rg_index(long & sample_i, long & rg_i, HtsRecord const & hts_rec) const
{
  assert(hts_rec.file_index < static_cast<long>(hts_files.size()));
  auto const & hts_f = hts_files[hts_rec.file_index];
  hts_f.get_sample_and_rg_index(sample_i, rg_i, hts_rec.record);
}


long
HtsParallelReader::get_num_rg() const
{
  return num_rg;
}


bam_hdr_t *
HtsParallelReader::get_header() const
{
  if (hts_files.size() == 0)
  {
    std::cerr << "[graphtyper::hts_parallel_reader] ERROR: Cannot get header with no opened files.\n";
    std::exit(1);
  }

  bam_hdr_t * hdr = nullptr;
  std::ostringstream ss;

  {
    bam_hdr_t * curr_hdr = hts_files[0].fp->bam_header;
    // Copy all the text
    ss << std::string(curr_hdr->text, curr_hdr->l_text);
  }

  // Add the @RG lines from other files
  for (long i = 1; i < static_cast<long>(hts_files.size()); ++i)
  {
    bam_hdr_t * curr_hdr = hts_files[i].fp->bam_header;
    std::string const rg = "@RG\t";
    const char * t = curr_hdr->text;
    const char * t_end = curr_hdr->text + (curr_hdr->l_text - 4); // 4 is the size of the header tags

    while (t < t_end)
    {
      long line_size = std::distance(t, std::find(t, t + curr_hdr->l_text, '\n'));

      if (std::equal(rg.begin(), rg.begin() + 4, t))
        ss << std::string(t, line_size) << '\n';

      t += line_size + 1;
    }
  }

  std::string new_header = ss.str();
  hdr = sam_hdr_parse(new_header.size(), new_header.c_str());
  hdr->l_text = new_header.size();
  hdr->text = static_cast<char *>(realloc(hdr->text, sizeof(char) * new_header.size()));
  strncpy(hdr->text, new_header.c_str(), new_header.size());
  return hdr;
}


//***
// Free functions

seqan::BamAlignmentRecord
parse(HtsRecord const & hts_rec)
{
  assert(hts_rec.record);

  seqan::BamAlignmentRecord rec;
  seqan::parse(rec, hts_rec.record);
  return rec;
}


void
get_sequence(seqan::IupacString & seq, seqan::IupacString & rseq, bam1_t const * rec)
{
  auto & core = rec->core;
  assert(core.l_qseq > 0);
  seqan::resize(seq, core.l_qseq);
  seqan::resize(rseq, core.l_qseq);

  auto i = bam_get_seq(rec);

  // Get the sequences
  for (int j = 0; j < core.l_qseq; ++j)
  {
    seq[j] = seq_nt16_str[bam_seqi(i, j)];
    rseq[j] = seq[j];
  }

  seqan::reverseComplement(rseq);
}


void
genotype_only(HtsParallelReader const & hts_preader,
              VcfWriter & writer,
              ReferenceDepth & reference_depth,
              std::vector<TMapGPaths> & maps,
              std::pair<GenotypePaths, GenotypePaths> & prev_paths,
              HtsRecord const & hts_rec,
              seqan::IupacString & seq,
              seqan::IupacString & rseq,
              bool update_prev_paths
              )
{
  assert(hts_rec.record);
  assert(hts_rec.file_index >= 0);
  assert(hts_rec.file_index < static_cast<long>(maps.size()));

  long rg_i = 0;
  long sample_i = 0;
  hts_preader.get_sample_and_rg_index(sample_i, rg_i, hts_rec);

  assert(maps.size() > 0);
  assert(sample_i < static_cast<long>(maps.size()));
  assert(rg_i < static_cast<long>(maps.size()));

  TMapGPaths & map_gpaths = maps[rg_i];
  std::string read_name(reinterpret_cast<char *>(hts_rec.record->data));

  if (update_prev_paths)
  {
    get_sequence(seq, rseq, hts_rec.record);
    prev_paths = align_read(hts_rec.record, seq, rseq);
  }

  std::pair<GenotypePaths, GenotypePaths> geno_paths(prev_paths);

  auto find_it = map_gpaths.find(read_name);

  if (find_it == map_gpaths.end())
  {
    if (hts_rec.record->core.flag & IS_PAIRED)
    {
      // Read is in a pair and we need to wait for its mate
#ifndef NDEBUG
      update_paths(geno_paths, seq, rseq, hts_rec.record);
#else
      update_paths(geno_paths, hts_rec.record);
#endif

      map_gpaths[read_name] = std::move(geno_paths); // Did not find the read name
    }
    else
    {
      GenotypePaths * selected = update_unpaired_read_paths(geno_paths, hts_rec.record);

      if (selected)
      {
        if (graph.is_sv_graph)
        {
          // Add reference depth
          reference_depth.add_genotype_paths(*selected, sample_i);
        }

        writer.update_haplotype_scores_geno(*selected, sample_i);
      }
    }
  }
  else
  {
    // Do stuff with paired reads....
#ifndef NDEBUG
    update_paths(geno_paths, seq, rseq, hts_rec.record);
#else
    update_paths(geno_paths, hts_rec.record);
#endif

    if ((geno_paths.first.flags & IS_FIRST_IN_PAIR) == (find_it->second.first.flags & IS_FIRST_IN_PAIR))
    {
      BOOST_LOG_TRIVIAL(error) << "Reads with name '" << read_name << "' both have IS_FIRST_IN_PAIR="
                               << ((geno_paths.first.flags & IS_FIRST_IN_PAIR) != 0);
      std::exit(1);
    }

    // Find the better pair of geno paths
    std::pair<GenotypePaths *, GenotypePaths *> better_paths =
      get_better_paths(find_it->second, geno_paths);

    if (better_paths.first)
    {
      if (graph.is_sv_graph)
      {
        // Add reference depth
        reference_depth.add_genotype_paths(*better_paths.first, sample_i);
        reference_depth.add_genotype_paths(*better_paths.second, sample_i);
      }

      writer.update_haplotype_scores_geno(*better_paths.first, sample_i);
      writer.update_haplotype_scores_geno(*better_paths.second, sample_i);
    }

    map_gpaths.erase(find_it); // Remove the read name afterwards
  }
}


void
genotype_and_discover(HtsParallelReader const & hts_preader,
                      VcfWriter & writer,
                      ReferenceDepth & reference_depth,
                      VariantMap & varmap,
                      std::vector<TMapGPaths> & maps,
                      std::pair<GenotypePaths, GenotypePaths> & prev_paths,
                      HtsRecord const & hts_rec,
                      seqan::IupacString & seq,
                      seqan::IupacString & rseq,
                      bool update_prev_paths
                      )
{
  assert(hts_rec.record);
  assert(hts_rec.file_index >= 0);

  long rg_i = 0;
  long sample_i = 0;
  hts_preader.get_sample_and_rg_index(sample_i, rg_i, hts_rec);

  assert(maps.size() > 0u);
  assert(sample_i < static_cast<long>(maps.size()));
  assert(rg_i < static_cast<long>(maps.size()));

  TMapGPaths & map_gpaths = maps[rg_i];
  std::string read_name(reinterpret_cast<char *>(hts_rec.record->data));

  if (update_prev_paths)
  {
    get_sequence(seq, rseq, hts_rec.record); // Updates seq and rseq
    prev_paths = align_read(hts_rec.record, seq, rseq);
  }

  std::pair<GenotypePaths, GenotypePaths> geno_paths(prev_paths);
  auto find_it = map_gpaths.find(read_name);

  if (find_it == map_gpaths.end())
  {
    if (hts_rec.record->core.flag & IS_PAIRED)
    {
      // Read is in a pair and we need to wait for its mate
#ifndef NDEBUG
      update_paths(geno_paths, seq, rseq, hts_rec.record);
#else
      update_paths(geno_paths, hts_rec.record);
#endif

      further_update_paths_for_discovery(geno_paths, seq, rseq, hts_rec.record);
      map_gpaths[read_name] = std::move(geno_paths); // Did not find the read name
    }
    else
    {
      GenotypePaths * selected = update_unpaired_read_paths(geno_paths, hts_rec.record);

      if (selected)
      {
        further_update_unpaired_read_paths_for_discovery(*selected, seq, rseq, hts_rec.record);

        // Discover new variants
        {
          std::vector<VariantCandidate> new_vars = selected->find_new_variants();

          if (new_vars.size() > 0)
            varmap.add_variants(std::move(new_vars), sample_i);
        }

        // Add reference depth
        reference_depth.add_genotype_paths(*selected, sample_i);

        // Update haplotype likelihood scores
        writer.update_haplotype_scores_geno(*selected, sample_i);
      }
    }
  }
  else
  {
    // Do stuff with paired reads....
#ifndef NDEBUG
    update_paths(geno_paths, seq, rseq, hts_rec.record);
#else
    update_paths(geno_paths, hts_rec.record);
#endif

    if ((geno_paths.first.flags & IS_FIRST_IN_PAIR) == (find_it->second.first.flags & IS_FIRST_IN_PAIR))
    {
      BOOST_LOG_TRIVIAL(error) << "Reads with name '" << read_name << "' both have IS_FIRST_IN_PAIR="
                               << ((geno_paths.first.flags & IS_FIRST_IN_PAIR) != 0);
      std::exit(1);
    }

    further_update_paths_for_discovery(geno_paths, seq, rseq, hts_rec.record);

    // Find the better pair of geno paths
    std::pair<GenotypePaths *, GenotypePaths *> better_paths =
      get_better_paths(find_it->second, geno_paths);

    if (better_paths.first)
    {
      // Discover new variants
      {
        std::vector<VariantCandidate> new_vars = better_paths.first->find_new_variants();

        if (new_vars.size() > 0)
          varmap.add_variants(std::move(new_vars), sample_i);

        new_vars = better_paths.second->find_new_variants();

        if (new_vars.size() > 0)
          varmap.add_variants(std::move(new_vars), sample_i);
      }

      // Add reference depth
      reference_depth.add_genotype_paths(*better_paths.first, sample_i);
      reference_depth.add_genotype_paths(*better_paths.second, sample_i);

      writer.update_haplotype_scores_geno(*better_paths.first, sample_i);
      writer.update_haplotype_scores_geno(*better_paths.second, sample_i);
    }

    map_gpaths.erase(find_it); // Remove the read name afterwards
  }
}


void
parallel_reader_genotype_only(std::string * out_path,
                              std::vector<std::string> const * hts_paths_ptr,
                              std::string const & output_dir,
                              std::string const & reference,
                              std::string const & region,
                              bool const is_writing_calls_vcf,
                              bool const is_writing_hap)
{
  assert(hts_paths_ptr);
  auto const & hts_paths = *hts_paths_ptr;

  // Inititalize the HTS parallel reader
  HtsParallelReader hts_preader;
  hts_preader.open(hts_paths, reference, region); // "" is reference, "." is the region

  // Set up VcfWriter
  VcfWriter writer(SPLIT_VAR_THRESHOLD - 1);
  writer.set_samples(hts_preader.get_samples());

  if (writer.pns.size() == 0)
  {
    BOOST_LOG_TRIVIAL(error) << "No samples were extracted.";
    std::exit(1);
  }

  std::string const & first_sample = writer.pns[0];

  // Set up reference depth tracks if we are SV calling
  ReferenceDepth reference_depth;

  if (graph.is_sv_graph)
    reference_depth.set_depth_sizes(writer.pns.size());

  std::vector<TMapGPaths> maps; // One map for each file. Each map relates read names to their graph alignments
  maps.resize(hts_preader.get_num_rg());

#ifndef NDEBUG
  if (Options::instance()->stats.size() > 0)
  {
    writer.print_statistics_headers();
    writer.print_variant_details();
    writer.print_variant_group_details();
  }
#endif // NDEBUG

  long num_records = 0;
  long num_duplicated_records = 0;
  std::pair<GenotypePaths, GenotypePaths> prev_paths;
  HtsRecord prev;

  // Read the first record
  if (!hts_preader.read_record(prev))
  {
    BOOST_LOG_TRIVIAL(debug) << "[graphtyper::utilities::hts_parallel_reader] No reads read.";
  }
  else
  {
    ++num_records;
    seqan::IupacString seq;
    seqan::IupacString rseq;
    genotype_only(hts_preader, writer, reference_depth, maps, prev_paths, prev, seq, rseq,
                  true /*update prev_geno_paths*/);
    HtsRecord curr;

    while (hts_preader.read_record(curr))
    {
      // Ignore reads with the following bits set
      if ((curr.record->core.flag & (IS_SECONDARY | IS_QC_FAIL | IS_DUPLICATION | IS_SUPPLEMENTARY)) != 0u)
        continue;

      ++num_records;

      if (equal_pos_seq(prev.record, curr.record))
      {
        // The two records are equal
        ++num_duplicated_records;
        genotype_only(hts_preader, writer, reference_depth, maps, prev_paths, curr, seq, rseq,
                      false /*update prev_geno_paths*/);
      }
      else
      {
        genotype_only(hts_preader, writer, reference_depth, maps, prev_paths, curr, seq, rseq,
                      true /*update prev_geno_paths*/);
        hts_preader.move_record(prev, curr); // move curr to prev
      }
    }

    BOOST_LOG_TRIVIAL(debug) << "[graphtyper::hts_parallel_reader] Num of duplicated records: "
                             << num_duplicated_records << " / " << num_records;
  }

#ifndef NDEBUG
  check_if_maps_are_empty(maps);
#endif // NDEBUG

  // Write haplotype calls
  if (is_writing_hap)
  {
    // Write genotype calls. No need to write haplotype calls in SV genotyping
    std::ostringstream hap_calls_path;
    hap_calls_path << output_dir << "/" << first_sample << ".hap";
    HaplotypeCalls hap_calls(writer.get_haplotype_calls());

    BOOST_LOG_TRIVIAL(debug) << "[graphtyper::hts_parallel_reader] Writing haplotype calls to '"
                             << hap_calls_path.str() << "'";
    save_calls(hap_calls, hap_calls_path.str());
  }

  // Optionally we can write _calls.vcf.gz files for each pool of samples
  if (is_writing_calls_vcf)
  {
    BOOST_LOG_TRIVIAL(debug) << "[graphtyper::hts_parallel_reader] Writing calls to '"
                             << output_dir << "/"
                             << first_sample << "_calls.vcf.gz'";

    std::ostringstream vcf_calls_path;
    vcf_calls_path << output_dir << "/" << first_sample << "_calls.vcf.gz";
    Vcf vcf(WRITE_BGZF_MODE, vcf_calls_path.str());

    // Set sample names
    vcf.sample_names = writer.pns;

    for (long ps = 0; ps < static_cast<long>(writer.haplotypes.size()); ++ps)
    {
      vcf.add_haplotype(writer.haplotypes[ps],
                        true /*clear haplotypes*/,
                        static_cast<uint32_t>(ps)
                        );
    }

    for (auto & var : vcf.variants)
      var.trim_sequences(false); // Don't keep one match

    // If the graph is an SV graph, then reformat the SV record
    reformat_sv_vcf_records(vcf.variants, reference_depth);

    // If there are samples then generate info fields
    if (vcf.sample_names.size() > 0)
    {
      for (auto & var : vcf.variants)
        var.generate_infos();
    }

    vcf.write();
  }

  *out_path = output_dir + "/" + first_sample;
}


void
parallel_reader_with_discovery(std::string * out_path,
                               std::vector<std::string> const * hts_paths_ptr,
                               std::string const & output_dir,
                               std::string const & reference,
                               std::string const & region,
                               long const minimum_variant_support,
                               double const minimum_variant_support_ratio,
                               bool const is_writing_calls_vcf,
                               bool const is_writing_hap)
{
  assert(hts_paths_ptr);
  auto const & hts_paths = *hts_paths_ptr;

  // Initialize the HTS parallel reader
  HtsParallelReader hts_preader;
  hts_preader.open(hts_paths, reference, region);

  // Set up VcfWriter
  VcfWriter writer(SPLIT_VAR_THRESHOLD - 1);
  writer.set_samples(hts_preader.get_samples());
  std::string const & first_sample = writer.pns[0];

  ReferenceDepth reference_depth;
  reference_depth.set_depth_sizes(writer.pns.size());

  VariantMap varmap;
  varmap.set_samples(hts_preader.get_samples());
  varmap.minimum_variant_support = minimum_variant_support;
  varmap.minimum_variant_support_ratio = minimum_variant_support_ratio;

  std::vector<TMapGPaths> maps; // One map for each file. Each map relates read names to their graph alignments
  maps.resize(hts_preader.get_num_rg());

#ifndef NDEBUG
  if (Options::instance()->stats.size() > 0)
  {
    writer.print_statistics_headers();
    writer.print_variant_details();
    writer.print_variant_group_details();
  }
#endif // NDEBUG

  long num_records = 0;
  long num_duplicated_records = 0;
  std::pair<GenotypePaths, GenotypePaths> prev_paths;
  HtsRecord prev;

  // Read the first record
  if (!hts_preader.read_record(prev))
  {
    BOOST_LOG_TRIVIAL(debug) << "[graphtyper::hts_parallel_reader] No reads were read.";
  }
  else
  {
    ++num_records;
    seqan::IupacString seq;
    seqan::IupacString rseq;
    genotype_and_discover(hts_preader, writer, reference_depth, varmap, maps, prev_paths, prev,
                          seq, rseq, true /*update prev_geno_paths*/);
    HtsRecord curr;

    while (hts_preader.read_record(curr))
    {
      // Ignore reads with the following bits set
      if ((curr.record->core.flag & (IS_SECONDARY | IS_QC_FAIL | IS_DUPLICATION | IS_SUPPLEMENTARY)) != 0u)
        continue;

      ++num_records;

      if (equal_pos_seq(prev.record, curr.record))
      {
        // The two records are equal
        ++num_duplicated_records;
        genotype_and_discover(hts_preader, writer, reference_depth, varmap, maps, prev_paths, curr,
                              seq, rseq, false /*update prev_geno_paths*/);
      }
      else
      {
        genotype_and_discover(hts_preader, writer, reference_depth, varmap, maps, prev_paths, curr,
                              seq, rseq, true /*update prev_geno_paths*/);
        hts_preader.move_record(prev, curr); // move curr to prev
      }
    }

    BOOST_LOG_TRIVIAL(debug) << "[graphtyper::hts_parallel_reader] Num of duplicated records: "
                             << num_duplicated_records << " / " << num_records;
  }

#ifndef NDEBUG
  check_if_maps_are_empty(maps);
#endif // NDEBUG

  // Output variants
  varmap.create_varmap_for_all(reference_depth);

#ifndef NDEBUG
  if (Options::instance()->stats.size() > 0)
    varmap.write_stats("2");
#endif // NDEBUG

  std::ostringstream variant_map_path;
  variant_map_path << output_dir << "/" << first_sample << "_variant_map";
  BOOST_LOG_TRIVIAL(debug) << "[graphtyper::caller] Writing variant map to '" << variant_map_path.str();
  save_variant_map(variant_map_path.str(), varmap);

  // Write haplotype calls
  if (is_writing_hap)
  {
    // Write genotype calls. No need to write haplotype calls in SV genotyping
    std::ostringstream hap_calls_path;
    hap_calls_path << output_dir << "/" << writer.pns[0] << ".hap";
    HaplotypeCalls hap_calls(writer.get_haplotype_calls());

    BOOST_LOG_TRIVIAL(debug) << "[graphtyper::hts_parallel_reader] Writing haplotype calls to '"
                             << hap_calls_path.str() << "'";
    save_calls(hap_calls, hap_calls_path.str());
  }

  // Optionally we can write _calls.vcf.gz files for each pool of samples
  if (is_writing_calls_vcf)
  {
    BOOST_LOG_TRIVIAL(debug) << "[graphtyper::hts_parallel_reader] Writing calls to '"
                             << output_dir << "/"
                             << writer.pns[0] << "_calls.vcf.gz'";

    std::ostringstream vcf_calls_path;
    vcf_calls_path << output_dir << "/" << first_sample << "_calls.vcf.gz";
    Vcf vcf(WRITE_BGZF_MODE, vcf_calls_path.str());

    // Set sample names
    vcf.sample_names = writer.pns;

    for (long ps = 0; ps < static_cast<long>(writer.haplotypes.size()); ++ps)
    {
      vcf.add_haplotype(writer.haplotypes[ps],
                        true /*clear haplotypes*/,
                        static_cast<uint32_t>(ps)
                        );
    }

    for (auto & var : vcf.variants)
      var.trim_sequences(false); // Don't keep one match

    // If the graph is an SV graph, then reformat the SV record
    reformat_sv_vcf_records(vcf.variants, reference_depth);

    // If there are samples then generate info fields
    if (vcf.sample_names.size() > 0)
    {
      for (auto & var : vcf.variants)
        var.generate_infos();
    }

    vcf.write();
  }

  *out_path = output_dir + "/" + first_sample;
}


void
sam_merge(std::string const & output_sam, std::vector<std::string> const & input_sams)
{
  {
    HtsParallelReader hts_preader;
    hts_preader.open(input_sams);

    HtsRecord hts_rec;
    HtsWriter hts_writer;
    hts_writer.open(output_sam);
    hts_writer.copy_header(hts_preader);

    while (hts_preader.read_record(hts_rec))
    {
      assert(hts_rec.record);
      hts_writer.write(hts_rec.record);
    }

    hts_writer.close();
    hts_preader.close();
  }

  // Remove old files
  for (std::string const & input_sam : input_sams)
  {
    int ret = unlink(input_sam.c_str());

    if (ret < 0)
    {
      BOOST_LOG_TRIVIAL(warning) << "[graphtyper::hts_parallel_reader] WARNING: Unable to remove " << input_sam;
    }
  }
}


} // namespace gyper
