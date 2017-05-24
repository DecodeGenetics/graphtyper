#include <iostream>
#include <string>

#include <boost/log/trivial.hpp>

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/graph.hpp>
#include <graphtyper/graph/var_record.hpp>
#include <graphtyper/graph/constructor.hpp>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>
#include <seqan/vcf_io.h>


namespace gyper
{

void
read_reference_genome(std::vector<char> & reference_sequence,
                      seqan::FaiIndex & fasta_index,
                      GenomicRegion const & genomic_region,
                      std::string const fasta_filename
                      )
{
  if (!seqan::open(fasta_index, fasta_filename.c_str()))
  {
    if (!seqan::build(fasta_index, fasta_filename.c_str()))
    {
      BOOST_LOG_TRIVIAL(error) << "[graphtyper::constructor] Index could not be loaded or built.";
    }
    else if (!seqan::save(fasta_index))
    {
      BOOST_LOG_TRIVIAL(error) << "[graphtyper::constructor] Index could not be saved to disk.";
    }
  }

  seqan::Dna5String ref_seq;
  unsigned idx = 0;

  if (!seqan::getIdByName(idx, fasta_index, genomic_region.chr.c_str()))
  {
    BOOST_LOG_TRIVIAL(error) << "[graphtyper::constructor] FAI index has no entry for this chromosome";
    assert(false);
  }

  seqan::readRegion(
    ref_seq,
    fasta_index,
    idx,
    genomic_region.begin,
    genomic_region.end
    );

  reference_sequence.clear();
  std::move(seqan::begin(ref_seq), seqan::end(ref_seq), std::back_inserter(reference_sequence));
}


void
open_tabix(seqan::Tabix & tabix_file, std::string const tabix_filename, GenomicRegion const & genomic_region)
{
  seqan::open(tabix_file, tabix_filename.c_str());
  seqan::String<char> header;
  seqan::getHeader(header, tabix_file);

  std::string region_str = genomic_region.to_string();
  seqan::setRegion(tabix_file, region_str.c_str());
}


VarRecord
vcf_to_var_record(seqan::VcfRecord const & vcf_record)
{
  assert(vcf_record.beginPos != seqan::VcfRecord::INVALID_POS);

  // std::cout << "VCF record start pos = " << vcf_record.beginPos << std::endl;

  std::vector<char> ref(0);
  std::move(seqan::begin(vcf_record.ref), seqan::end(vcf_record.ref), std::back_inserter(ref));

  std::vector<std::vector<char> > alts(0);
  seqan::StringSet<seqan::String<char> > seqan_alts;
  seqan::strSplit(seqan_alts, vcf_record.alt, seqan::EqualsChar<','>());

  for (unsigned i = 0; i < seqan::length(seqan_alts); ++i)
  {
    std::vector<char> alt(0);
    std::move(seqan::begin(seqan_alts[i]), seqan::end(seqan_alts[i]), std::back_inserter(alt));
    alts.push_back(std::move(alt));
  }

  return VarRecord(std::move(vcf_record.beginPos), std::move(ref), std::move(alts));
}


void
construct_graph(std::string const & reference_filename,
                std::string const & vcf_filename,
                std::vector<std::string> const & regions,
                bool const use_absolute_positions
                )
{
  seqan::FaiIndex fasta_index;
  graph = Graph(use_absolute_positions);

  for (auto region_it = regions.begin(); region_it != regions.end(); ++region_it)
  {
    BOOST_LOG_TRIVIAL(info) << "[graphtyper::constructor] Constructing graph for region " << *region_it;
    GenomicRegion genomic_region(*region_it);

    BOOST_LOG_TRIVIAL(info) << "[graphtyper::constructor] Reading FASTA file located at " << reference_filename;
    std::vector<char> reference_sequence;
    read_reference_genome(reference_sequence, fasta_index, genomic_region, reference_filename);

    BOOST_LOG_TRIVIAL(info) << "[graphtyper::constructor] INFO: Reading VCF file located at " << vcf_filename;
    seqan::Tabix tabix_file;
    open_tabix(tabix_file, vcf_filename, genomic_region);

    seqan::VcfRecord vcf_record;
    bool success = seqan::readRegion(vcf_record, tabix_file);

    // std::cout << "TEST " << vcf_record.beginPos << " " << static_cast<int32_t>(genomic_region.begin) << std::endl;

    // Ignore VCF record that overlap the region, but start earlier
    while (success && vcf_record.beginPos < static_cast<int32_t>(genomic_region.begin))
    {
      success = seqan::readRegion(vcf_record, tabix_file);
    }

    std::vector<VarRecord> var_records(0);

    if (success)
    {
      var_records.push_back(vcf_to_var_record(vcf_record));
    }

    while (seqan::readRegion(vcf_record, tabix_file))
    {
      if (vcf_record.beginPos + seqan::length(vcf_record.ref) < genomic_region.end)
        var_records.push_back(vcf_to_var_record(vcf_record));
    }

    //BOOST_LOG_TRIVIAL(info) << "Read " << var_records.size() << " records from region " << *region_it;
    genomic_region.check_if_var_records_match_reference_genome(var_records, reference_sequence);

    for (auto & var_record : var_records)
      genomic_region.add_reference_to_record_if_they_have_a_matching_prefix(var_record, reference_sequence);

    genomic_region.check_if_var_records_match_reference_genome(var_records, reference_sequence);
    graph.add_genomic_region(std::move(reference_sequence), std::move(var_records), std::move(genomic_region));
    assert(graph.size() > 0);

    if (!graph.check())
    {
      BOOST_LOG_TRIVIAL(error) << "[graphtyper::graph] Problem creating graph.";
      std::exit(1);
    }
  }

  BOOST_LOG_TRIVIAL(info) << "[graphtyper::constructor] Graph was successfully constructed.";

  // Create all specials positions
  graph.create_special_positions();
}


void
construct_graph(std::string const & reference_filename, std::vector<std::string> const & regions, bool const use_absolute_positions)
{
  seqan::FaiIndex fasta_index;
  graph = Graph(use_absolute_positions);

  for (auto region_it = regions.begin(); region_it != regions.end(); ++region_it)
  {
    GenomicRegion genomic_region(*region_it);

    std::vector<char> reference_sequence;
    read_reference_genome(reference_sequence, fasta_index, genomic_region, reference_filename);
    graph.add_genomic_region(std::move(reference_sequence), std::vector<VarRecord>(0), std::move(genomic_region));
  }

  assert(graph.size() > 0);
}


} // namespace gyper
