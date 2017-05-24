#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

#include <graphtyper/graph/graph_serialization.hpp>
#include <graphtyper/index/rocksdb.hpp>
#include <graphtyper/index/indexer.hpp>
#include <graphtyper/utilities/io.hpp>
#include <graphtyper/typer/alignment.hpp>


namespace gyper
{

void
align_reads_from_fasta(std::string graph_path, std::string index_path, std::string fasta_filename)
{
  std::vector<std::pair<seqan::CharString, seqan::Dna5String> > fasta_sequences =
    gyper::read_fasta_sequences(fasta_filename);

  load_graph(graph_path);
  assert(graph.size() > 0);
  load_index(index_path);

  std::cout.width(80);
  std::cout << "ID";
  std::cout.width(10);
  std::cout << "START";
  std::cout.width(10);
  std::cout << "END";
  std::cout.width(13);
  std::cout << "START_INDEX";
  std::cout.width(11);
  std::cout << "END_INDEX";
  std::cout.width(8);
  std::cout << "LENGTH";
  std::cout.width(10);
  std::cout << "NUM_VARS";
  std::cout << std::endl;

  for (auto it = fasta_sequences.cbegin(); it != fasta_sequences.cend(); ++it)
  {
    // it->first is ID, it->second is the sequence
    // std::cout << "[aligner] Aligning sequence of ID " << it->first << std::endl;
    std::vector<seqan::Dna5String> sequences = {it->second};
    std::vector<GenotypePaths> gt_paths = find_haplotype_paths(sequences);
    assert(gt_paths.size() == 1);
    GenotypePaths const & geno = gt_paths[0];

    if (geno.longest_path_size() == 0)
    {
      std::cout << "[aligner] INFO: No alignment found for sequence with ID = " << it->first << std::endl;
    }
    else if (geno.longest_path_size() <= 1)
    {
      std::cout.width(80);
      std::cout << std::string(seqan::begin(it->first), seqan::end(it->first));
      std::cout.width(10);
      std::cout << "NA";
      std::cout.width(10);
      std::cout << "NA";
      std::cout.width(13);
      std::cout << "NA";
      std::cout.width(11);
      std::cout << "NA";
      std::cout.width(8);
      std::cout << "NA";
      std::cout.width(10);
      std::cout << "NA";
      std::cout << std::endl;
    }
    else
    {
      for (auto const & path : geno.paths)
      {
        std::cout.width(80);
        std::cout << std::string(seqan::begin(it->first), seqan::end(it->first));
        std::cout.width(10);
        std::cout << path.start_ref_reach_pos();
        std::cout.width(10);
        std::cout << path.end_ref_reach_pos();
        std::cout.width(13);
        std::cout << path.read_start_index;
        uint32_t const length = std::min(static_cast<std::size_t>(K + (path.size() - 1) * (K - 1)), seqan::length(it->second));
        std::cout.width(11);
        std::cout << path.read_start_index + length;
        std::cout.width(8);
        std::cout << length;
        std::cout.width(10);
        std::cout << path.var_order.size();
        std::cout << std::endl;
      }
    }
  }
}


} // namespace gyper
