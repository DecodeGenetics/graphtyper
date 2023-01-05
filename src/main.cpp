#include <fstream>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <thread>
#include <unistd.h>
#include <vector>

#include <paw/parser.hpp>

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/constructor.hpp>
#include <graphtyper/graph/graph_serialization.hpp>
#include <graphtyper/typer/vcf.hpp>
#include <graphtyper/typer/vcf_operations.hpp>
#include <graphtyper/utilities/bamshrink.hpp>
#include <graphtyper/utilities/genotype.hpp>
#include <graphtyper/utilities/genotype_camou.hpp>
#include <graphtyper/utilities/genotype_hla.hpp>
#include <graphtyper/utilities/genotype_lr.hpp>
#include <graphtyper/utilities/genotype_sv.hpp>
#include <graphtyper/utilities/io.hpp> // gyper::get_contig_to_lengths
#include <graphtyper/utilities/logging.hpp>
#include <graphtyper/utilities/options.hpp>
#include <graphtyper/utilities/system.hpp>

namespace
{
void add_region(std::unordered_map<std::string, long> const & contig2length,
                std::vector<gyper::GenomicRegion> & regions,
                std::string const & region_str,
                long const REGION_SIZE)
{
  gyper::GenomicRegion region(region_str);
  auto find_it = contig2length.find(region.chr);

  if (find_it == contig2length.end())
  {
    print_log(gyper::log_severity::warning, "Unable to find contig: ", region.chr, " in reference.");
    return;
  }

  // Make sure the region is inside the contig
  region.begin = std::min(region.begin, static_cast<long>(find_it->second - 1));
  region.end = std::min(region.end, static_cast<long>(find_it->second));

  // Split region if it is too large (more than 10% of given REGION_SIZE)
  while (static_cast<long>(region.end - region.begin) > (REGION_SIZE + REGION_SIZE / 10l))
  {
    gyper::GenomicRegion new_region(region);
    new_region.end = new_region.begin + REGION_SIZE;
    region.begin = new_region.end;
    regions.push_back(std::move(new_region));
  }

  regions.push_back(std::move(region));
}

std::vector<gyper::GenomicRegion> get_regions(std::string const & ref_fn,
                                              std::string const & opts_region,
                                              std::string const & opts_region_fn,
                                              long const REGION_SIZE)
{
  std::unordered_map<std::string, long> contig2length = gyper::get_contig_to_lengths(ref_fn + ".fai");

  // region arg
  std::vector<gyper::GenomicRegion> regions;

  if (opts_region.size() > 0)
    add_region(contig2length, regions, opts_region, REGION_SIZE);

  if (opts_region_fn.size() > 0)
  {
    std::string line;
    std::ifstream file_in;
    file_in.open(opts_region_fn);

    while (std::getline(file_in, line))
      add_region(contig2length, regions, line, REGION_SIZE);
  }

  if (regions.size() == 0)
  {
    print_log(gyper::log_severity::error,
              "No regions specified. Either use --region or --region_file option to specify regions.");
    std::exit(1);
  }

  return regions;
}

std::vector<std::string> get_sams(std::string const & opts_sam, std::string const & opts_sams_file)
{
  std::vector<std::string> sams;

  if (opts_sam.size() > 0)
  {
    std::ifstream sam_in(opts_sam);

    if (!sam_in.is_open())
    {
      print_log(gyper::log_severity::error, "Could not open sam file '", opts_sam, "'");
      std::exit(1);
    }

    sams.push_back(opts_sam);
  }

  if (opts_sams_file.size() > 0)
  {
    std::ifstream file_in(opts_sams_file);

    if (!file_in.is_open())
    {
      print_log(gyper::log_severity::error, "Could not open file with SAM/BAM/CRAM paths '", opts_sams_file, "'");
      std::exit(1);
    }

    std::string line;

    while (std::getline(file_in, line))
    {
      if (std::find(line.cbegin(), line.cend(), '\0') != line.cend())
      {
        print_log(gyper::log_severity::error,
                  __HERE__,
                  " Unexpectedly found NULL in input file with SAM/BAM/CRAM paths --sams=",
                  opts_sams_file,
                  ". The file should be uncompressed and contain only file paths.");
        std::exit(1);
      }

      sams.push_back(line);
    }
  }

  if (sams.size() == 0)
  {
    print_log(gyper::log_severity::error, __HERE__, " No SAM/BAM/CRAM files where given as input.");
    std::exit(1);
  }

  return sams;
}

std::vector<double> get_avg_cov_by_readlen(std::string const & avg_cov_by_readlen_fn, long const expected_num)
{
  std::vector<double> avg_cov_by_readlen;

  if (avg_cov_by_readlen_fn.size() == 0)
  {
    avg_cov_by_readlen.resize(expected_num, -1.0);
  }
  else
  {
    std::ifstream ifs(avg_cov_by_readlen_fn);

    if (!ifs.is_open())
    {
      print_log(gyper::log_severity::error, "Could not open avgCovByReadlen file ", avg_cov_by_readlen_fn);
      std::exit(1);
    }

    std::string line;

    while (std::getline(ifs, line))
    {
      double const cov_by_readlen = std::stod(line);
      avg_cov_by_readlen.push_back(cov_by_readlen);
    }

    if (static_cast<long>(avg_cov_by_readlen.size()) != expected_num)
    {
      print_log(gyper::log_severity::error,
                "avg_cov_by_readlen file should have the same number of lines as there are SAMs. ",
                avg_cov_by_readlen.size(),
                " != ",
                expected_num);
      std::exit(1);
    }
  }

  return avg_cov_by_readlen;
}

std::vector<std::string> get_sam_index_paths(std::string const & sam_index_fn, std::vector<std::string> const & sams_fn)
{
  std::vector<std::string> sam_index_paths;

  if (sam_index_fn.size() > 0)
  {
    // sam_index argument was given
    std::ifstream sam_index_f(sam_index_fn);

    if (!sam_index_f.is_open())
    {
      print_log(gyper::log_severity::error, "Could not open sam_index_fn ", sam_index_fn);
      std::exit(1);
    }

    std::string line;

    while (std::getline(sam_index_f, line))
      sam_index_paths.push_back(std::move(line));
  }
  else
  {
    // sam_index argument was not given
    sam_index_paths.reserve(sams_fn.size());

    for (auto const & path_in : sams_fn)
    {
      if (path_in.size() > 5 && std::string(path_in.rbegin(), path_in.rbegin() + 5) == "marc.")
        sam_index_paths.push_back(path_in + ".crai");
      else
        sam_index_paths.push_back(path_in + ".bai");
    }
  }

  if (sam_index_paths.size() != sams_fn.size())
  {
    print_log(gyper::log_severity::error,
              "ERROR: Number if sam_index paths doesn't match the number of BAM/CRAMs ",
              sam_index_paths.size(),
              " vs. ",
              sams_fn.size());
    std::exit(1);
  }

  return sam_index_paths;
}

void setup_logger()
{
  gyper::Options & opts = *(gyper::Options::instance());

  gyper::log_severity severity{};

  if (opts.vverbose)
    severity = gyper::log_severity::debug;
  else if (opts.verbose)
    severity = gyper::log_severity::info;
  else
    severity = gyper::log_severity::warning;

  if (opts.log.size() == 0 || opts.log == "-")
    gyper::log_singleton = std::unique_ptr<gyper::log_singleton_t>{new gyper::log_singleton_t{severity, std::clog}};
  else
    gyper::log_singleton = std::unique_ptr<gyper::log_singleton_t>{new gyper::log_singleton_t{severity, opts.log}};
}

int subcmd_bamshrink(paw::Parser & parser)
{
  bamshrink::Options opts;

  parser.parse_positional_argument(opts.bamPathIn, "bamPathIn", "Input BAM file path.");
  parser.parse_option(opts.avgCovByReadLen, 'a', "avg-cov-by-readlen", "Average coverage divided by read length.", "D");

  parser.parse_option(opts.bamIndex, ' ', "index", "Input BAM bai/CRAM crai index file.");
  parser.parse_option(opts.bamPathOut, 'o', "output", "Output BAM file.");
  parser.parse_option(opts.interval,
                      'i',
                      "interval",
                      "Interval/region to filter on in format chrA:N-M, where chrA is the contig name, "
                      "N is the begin position, and M is the end position of the interval.");

  parser.parse_option(opts.intervalFile, 'I', "interval-file", "File with interval(s)/region(s) to filter on.");

  parser.parse_option(opts.maxFragLen, 'f', "max-fragment-length", "Maximum fragment length allowed.", "N");
  parser.parse_option(opts.minNumMatching, 'm', "min-num-matching", "Minumum number of matching bases in read.", "N");

  parser.finalize();
  setup_logger();
  return bamshrink::main(opts);
}

int subcmd_check(paw::Parser & parser)
{
  std::string graph_fn;
  parser.parse_positional_argument(graph_fn, "GRAPH", "Path to graph.");
  parser.finalize();
  setup_logger();

  gyper::load_graph(graph_fn);
  gyper::graph.check();
  gyper::graph.print();

  return 0;
}

int subcmd_construct(paw::Parser & parser)
{
  gyper::Options & opts = *(gyper::Options::instance());

  std::string graph_fn;
  bool is_sv_graph{false};
  bool use_tabix{false};
  std::string ref_fn;
  std::string region;
  std::string vcf_fn;

  // Parse options
  parser.parse_option(is_sv_graph, ' ', "sv_graph", "Set to construct an SV graph.");
  parser.parse_option(opts.add_all_variants,
                      ' ',
                      "add_all_variants",
                      "Set to create a graph with every possible "
                      "haplotype on overlapping variants.");
  parser.parse_option(use_tabix, ' ', "use_tabix", "Set to use tabix index to extract variants of the given region.");
  parser.parse_option(vcf_fn, ' ', "vcf", "VCF variant input.");

  parser.parse_positional_argument(graph_fn, "GRAPH", "Path to graph.");
  parser.parse_positional_argument(ref_fn, "REF.FA", "Reference genome in FASTA format.");
  parser.parse_positional_argument(region, "REGION", "Genomic region to construct graph for.");

  parser.finalize();
  setup_logger();

  region.erase(std::remove(region.begin(), region.end(), ','), region.end());
  print_log(gyper::log_severity::info, "Constructing a graph for region ", region);

  gyper::check_file_exists(ref_fn);
  gyper::check_file_exists_or_empty(vcf_fn);

  gyper::construct_graph(ref_fn, vcf_fn, region, is_sv_graph, use_tabix);
  gyper::save_graph(graph_fn);
  print_log(gyper::log_severity::info, "Graph saved at ", graph_fn);
  return 0;
}

int subcmd_genotype(paw::Parser & parser)
{
  gyper::Options & opts = *(gyper::Options::instance());
  bool see_advanced_options{false};

  std::string avg_cov_by_readlen_fn;
  std::string sam_index_fn;
  std::string output_dir{"results"};
  std::string opts_region{};
  std::string opts_region_file{};
  std::string ref_fn{};
  std::string sam{};
  std::string sams{};
  std::string encoding{"vcf"};
  bool force_copy_reference{false};
  bool force_no_copy_reference{false};
  bool no_filter_on_proper_pairs{false};
  bool no_filter_on_read_bias{false};
  bool no_filter_on_strand_bias{false};

  // Parse options
  parser.parse_option(see_advanced_options,
                      'a',
                      "advanced",
                      "Set to enable advanced options. "
                      "See a list of all options (including advanced) with 'graphtyper genotype --advanced --help'");

  parser.parse_option(
    avg_cov_by_readlen_fn,
    ' ',
    "avg_cov_by_readlen",
    "File with average coverage by read length (one value per line). "
    "The values are used for subsampling regions with extremely high coverage and should be in the same "
    "order as the BAM/CRAM list.");

  parser.parse_option(opts.no_decompose,
                      ' ',
                      "no_decompose",
                      "Set to prohibit decomposing variants in VCF output, which means complex variants wont be "
                      "split/decomposed into smaller variants.");

  parser.parse_option(force_copy_reference,
                      ' ',
                      "force_copy_reference",
                      "Force copy of the reference FASTA to temporary folder.");

  parser.parse_option(
    force_no_copy_reference,
    ' ',
    "force_no_copy_reference",
    "Force that the reference FASTA is NOT copied to temporary folder. "
    "Useful if you have limit storage on your local disk or the reference is already on your local disk.");

  parser.parse_option(output_dir,
                      'O',
                      "output",
                      "Output directory. Results will be written in <output>/<contig>/<region>.vcf.gz");

  parser.parse_option(
    opts.prior_vcf,
    'p',
    "prior_vcf",
    "Input VCF file with prior variants sites. With this option set GraphTyper will be run normally except the given "
    "input variant sites are used in constructing the initial graph. We recommend only using common variants as a "
    "prior (1% allele frequency or higher). Note that the final output may not necessarily include every prior variant "
    "and GraphTyper may discover other variants as well. If you rather want to call only a specific set of variants "
    "then use the --vcf option instead.");

  parser.parse_option(opts_region,
                      'r',
                      "region",
                      "Genomic region to genotype. Use --region_file if you have more than one region.");

  parser.parse_option(opts_region_file,
                      'R',
                      "region_file",
                      "File with a list of genomic regions to genotype (one per line).");

  parser.parse_option(opts.threads,
                      't',
                      "threads",
                      "Max. number of threads to use. Note that it is not possible to utilize more threads than input "
                      "BAM/CRAMs.");

  parser.parse_option(
    sam,
    's',
    "sam",
    "Input BAM/CRAM to analyze. If you have more than one file then create a list and use --sams instead.");

  parser.parse_option(sams, 'S', "sams", "File with BAM/CRAMs paths to analyze (one per line).");

  parser.parse_option(
    sam_index_fn,
    'i',
    "sams_index",
    "File containing a list of BAM/CRAM indices to use when BAM/CRAM files are queried (one per line). "
    "The index files must be given in the same order as the BAMs/CRAMs.");

  parser.parse_option(opts.vcf,
                      ' ',
                      "vcf",
                      "Input VCF file with variant sites. "
                      "Use this option if you want GraphTyper to only genotype variants from this VCF.");

  if (see_advanced_options)
    parser.see_advanced_options(true);

  parser.parse_advanced_option(encoding, ' ', "encoding", "Select output encoding. Available are: vcf, popvcf");

  // set default compression level as 9 when popvcf, since level 9 is already very fast anyway in that encoding mode
  if (encoding == "popvcf")
    opts.bgzf_compression_level = 9;

  parser.parse_advanced_option(opts.no_asterisks, ' ', "no_asterisks", "Set to avoid using asterisk in VCF output.");

  parser.parse_advanced_option(opts.no_bamshrink, ' ', "no_bamshrink", "(advanced) Set to skip bamShrink.");

  parser.parse_advanced_option(
    opts.no_sample_name_reordering,
    ' ',
    "no_sample_name_reordering",
    "(advanced) Set to prohibit graphtyper from reordering input SAM/BAM/CRAM files. With this option set, the output "
    "VCF should contain samples in the same order as in the input.");

  parser.parse_advanced_option(opts.no_cleanup,
                               ' ',
                               "no_cleanup",
                               "(advanced) Set to skip removing temporary files. Useful for debugging.");

  parser.parse_advanced_option(opts.get_sample_names_from_filename,
                               ' ',
                               "get_sample_names_from_filename",
                               "(advanced) Set to ignore sample names from RG tags and attempt to retrive the sample "
                               "names from filenames instead.");

  parser.parse_advanced_option(opts.is_all_biallelic,
                               ' ',
                               "is_all_biallelic",
                               "(advanced) Set to force all output variants to be biallelic in the VCF output.");

  parser.parse_advanced_option(opts.is_sam_merging_allowed,
                               ' ',
                               "is_sam_merging_allowed",
                               "Set to allow SAM files to be merged interally in graphtyper. For this to be possible "
                               "you must make sure they all have read groups and their name is never duplicates.");

  parser.parse_advanced_option(no_filter_on_proper_pairs,
                               ' ',
                               "no_filter_on_proper_pairs",
                               "Set to disable filter on proper pairs when discovering potentially new variants. This "
                               "should be set if you have unpaired reads.");

  parser.parse_advanced_option(no_filter_on_read_bias,
                               ' ',
                               "no_filter_on_read_bias",
                               "Set to disable filter on read bias.");

  parser.parse_advanced_option(no_filter_on_strand_bias,
                               ' ',
                               "no_filter_on_strand_bias",
                               "Set to disable filter on strand bias.");

  parser.parse_advanced_option(opts.no_filter_on_coverage,
                               ' ',
                               "no_filter_on_coverage",
                               "Set to disable filter on coverage.");

  parser.parse_advanced_option(opts.no_filter_on_begin_pos,
                               ' ',
                               "no_filter_on_begin_pos",
                               "Set to disable filter on number of unique begin position of reads.");

  parser.parse_advanced_option(opts.no_variant_overlapping,
                               ' ',
                               "no_variant_overlapping",
                               "Set to avoid that variants overlap in the VCF output");

  parser.parse_advanced_option(opts.force_no_filter_zero_qual,
                               ' ',
                               "force_no_filter_zero_qual",
                               "Set to force variants to be in the final output despite they have zero quality"
                               " (all calls are ref/ref)");

  parser.parse_advanced_option(opts.normal_and_no_variant_overlapping,
                               ' ',
                               "normal_and_no_variant_overlapping",
                               "Set to output two files for each region, both normal (overlapping)"
                               " and non-overlapping variants.");

  parser.parse_advanced_option(opts.max_files_open,
                               ' ',
                               "max_files_open",
                               "Select how many input files are allowed to be open at the same time. "
                               "See your current limit with 'ulimit -n'.");

  parser.parse_advanced_option(opts.bamshrink_max_fraglen,
                               ' ',
                               "bamshrink_max_fraglen",
                               "Maximum fragment length for bamShrink to consider reads in proper pair.");

  parser.parse_advanced_option(opts.bamshrink_min_matching,
                               ' ',
                               "bamshrink_min_matching",
                               "Minimum amount of matches in reads.");

  parser.parse_advanced_option(opts.bamshrink_is_not_filtering_mapq0,
                               ' ',
                               "bamshrink_is_not_filtering_mapq0",
                               "Set to turn of bamShrink's filtering of MAPQ==0 reads.");

  parser.parse_advanced_option(opts.bamshrink_min_readlen,
                               ' ',
                               "bamshrink_min_readlen",
                               "Minimum read length after adapter trimming.");

  parser.parse_advanced_option(opts.bamshrink_min_readlen_low_mapq,
                               ' ',
                               "bamshrink_min_readlen_low_mapq",
                               "Minimum read length after adapter trimming for low MAPQ reads.");

  parser.parse_advanced_option(opts.bamshrink_min_unpair_readlen,
                               ' ',
                               "bamshrink_min_unpair_readlen",
                               "Minimum read length after adapter trimming for unpaired reads.");

  parser.parse_advanced_option(opts.bamshrink_as_filter_threshold,
                               ' ',
                               "bamshrink_as_filter_threshold",
                               "Threshold for alignment score filter. Lower value is stricter.");

  parser.parse_advanced_option(opts.force_use_input_ref_for_cram_reading,
                               ' ',
                               "force_use_input_ref_for_cram_reading",
                               "Force using the input reference FASTA file when reading CRAMs.");

  parser.parse_advanced_option(opts.genotype_aln_min_support,
                               ' ',
                               "genotype_aln_min_support",
                               "Minimum alignment support for a variant so it can be added to the graph.");

  parser.parse_advanced_option(opts.genotype_aln_min_support_ratio,
                               ' ',
                               "genotype_aln_min_support_ratio",
                               "Minimum alignment support ratio for a variant so it can be added to the graph.");

  parser.parse_advanced_option(opts.genotype_dis_min_support,
                               ' ',
                               "genotype_dis_min_support",
                               "Minimum graph discovery support for a variant so it can be added to the graph.");

  parser.parse_advanced_option(opts.genotype_dis_min_support_ratio,
                               ' ',
                               "genotype_dis_min_support_ratio",
                               "Minimum graph discovery support ratio for a variant so it can be added to the graph.");

  parser.parse_advanced_option(opts.is_csi, 'C', "csi", "If set, graphtyper will make csi indices instead of tbi.");

  parser.parse_advanced_option(opts.is_only_cigar_discovery,
                               ' ',
                               "is_only_cigar_discovery",
                               "If set, graphtyper will only discover variants from the aligner via the cigar.");

  parser.parse_advanced_option(opts.is_discovery_only_for_paired_reads,
                               ' ',
                               "is_discovery_only_for_paired_reads",
                               "If set, graphtyper will only discover variants from paired reads via the cigar.");

  parser.parse_advanced_option(
    opts.impurity_threshold,
    ' ',
    "impurity_threshold",
    "Set a threshold of variant impurity. Should be in range ]0.0, 0.25] where a lower number is stricter.");

  parser.parse_advanced_option(opts.minimum_extract_variant_support,
                               ' ',
                               "minimum_extract_variant_support",
                               "Minimum support for a variant for it to be allowed to go into the next iteration.");

  parser.parse_advanced_option(opts.minimum_extract_score_over_homref,
                               ' ',
                               "minimum_extract_score_over_homref",
                               "Minimum likelihood score over homozygous reference needed for a variant for it to be "
                               "allowed to go into the next iteration.");

  parser.parse_advanced_option(
    opts.primer_bedpe,
    ' ',
    "primer_bedpe",
    "In case you have PCR amplicons sequence data you may specify here a BEDPE file containing "
    "coordinates of the primer pairs. The primer sequence will be used in GraphTyper's alignment but are "
    "ignored in the genotyping model.");

  parser.parse_advanced_option(
    opts.sam_flag_filter,
    ' ',
    "sam_flag_filter",
    "ANY of these bits set in reads will be completely ignored by GraphTyper. Use with care.");

  parser.parse_advanced_option(opts.num_alleles_in_batch,
                               ' ',
                               "num_alleles_in_batch",
                               "How many alleles each batch of internal VCFs has. Increasing this number increases "
                               "memory used, while slightly decreases compute time.");

  parser.parse_advanced_option(opts.is_extra_call_only_iteration,
                               ' ',
                               "is_extra_call_only_iteration",
                               "Set to add an extra call only iteration.");

  parser.parse_advanced_option(opts.is_dropping_genotypes,
                               ' ',
                               "is_dropping_genotypes",
                               "If set, graphtyper will not output genotype calls (only sites and INFO).");

  parser.parse_advanced_option(opts.uncompressed_sample_names,
                               ' ',
                               "uncompressed_sample_names",
                               "If set, graphtyper BGZF VCF (.vcf.gz) files will output samples names in uncompressed "
                               "(0-level) in the header. The byte range of these blocks will also be printed in "
                               "${prefix}.samples_byte_range. Requires external utils: awk, truncate, stat.");

  parser.parse_advanced_option(opts.bgzf_compression_level,
                               ' ',
                               "bgzf_compression_level",
                               "Select bgzf compression level, ranges from 0 (for lowest and fastest) to 9 (highest "
                               "and slowest). Set as -1 to use default bgzf value.");

#ifndef NDEBUG
  parser.parse_advanced_option(opts.stats, ' ', "stats", "Directory for statistics files.");
#endif // ifndef NDEBUG

  if (opts.force_no_filter_zero_qual)
    opts.force_no_filter_bad_alts = true;

  parser.parse_positional_argument(ref_fn, "REF.FA", "Reference genome in FASTA format.");
  parser.finalize();
  setup_logger();

  opts.encoding = (encoding == "popvcf") ? 'p' : 'v';
  opts.filter_on_proper_pairs = !no_filter_on_proper_pairs;
  opts.filter_on_read_bias = !no_filter_on_read_bias;
  opts.filter_on_strand_bias = !no_filter_on_strand_bias;
  print_log(gyper::log_severity::info, "Running the 'genotype' subcommand.");

#ifndef NDEBUG
  // Create stats directory if it doesn't exist
  if (opts.stats.size() > 0 && !gyper::is_directory(opts.stats))
    mkdir(opts.stats.c_str(), 0755);
#endif // ifndef NDEBUG

  // Get the genomic regions to process from the --region and --region_file options
  std::vector<gyper::GenomicRegion> regions = get_regions(ref_fn, opts_region, opts_region_file, 50000);

  // Get the SAM/BAM/CRAM file names
  std::vector<std::string> sams_fn = get_sams(sam, sams);
  long const NUM_SAMPLES = sams_fn.size();

  // If neither force copy reference or force no copy reference we determine it from number of SAMs and
  // whether REF_CACHE is set or not
  bool is_copy_reference =
    force_copy_reference || (!force_no_copy_reference && NUM_SAMPLES >= 100 && !gyper::is_defined_in_env("REF_CACHE"));

  // Get the avgCovByReadLen for each of the SAM/BAM/CRAM
  std::vector<double> avg_cov_by_readlen = get_avg_cov_by_readlen(avg_cov_by_readlen_fn, NUM_SAMPLES);
  std::vector<std::string> sam_index_paths = get_sam_index_paths(sam_index_fn, sams_fn);

  gyper::genotype_regions(ref_fn, sams_fn, sam_index_paths, regions, output_dir, avg_cov_by_readlen, is_copy_reference);
  return 0;
}

int subcmd_genotype_hla(paw::Parser & parser)
{
  gyper::Options & opts = *(gyper::Options::instance());

  std::string avg_cov_by_readlen_fn;
  std::string output_dir{"hla_results"};
  std::string opts_region{};
  std::string opts_region_file{};
  std::string ref_fn{};
  std::string sam{};
  std::string sams{};
  std::string sam_index_fn{};
  std::string hla_vcf{};
  bool force_copy_reference{false};
  bool force_no_copy_reference{false};
  bool see_advanced_options{false};
  std::string interval_fn{};

  // Change the default values
  opts.max_files_open = 1024;

  // Parse options
  parser.parse_option(
    see_advanced_options,
    'a',
    "advanced",
    "Set to enable advanced options. "
    "See a list of all options (including advanced) with 'graphtyper genotype_hla --advanced --help'");

  parser.parse_option(
    avg_cov_by_readlen_fn,
    ' ',
    "avg_cov_by_readlen",
    "File with average coverage by read length (one value per line). "
    "The values are used for subsampling regions with extremely high coverage and should be in the same "
    "order as the BAM/CRAM list.");

  parser.parse_option(interval_fn,
                      'i',
                      "intervals",
                      "BED file with intervals to gather reads from. If empty, reads in region will be used.");

  parser.parse_option(output_dir, 'O', "output", "Output directory.");
  parser.parse_option(opts_region, 'r', "region", "Genomic region to genotype.");
  parser.parse_option(opts_region_file, 'R', "region_file", "File with genomic regions to genotype.");
  parser.parse_option(opts.threads, 't', "threads", "Max. number of threads to use.");
  parser.parse_option(sam, 's', "sam", "SAM/BAM/CRAM to analyze.");
  parser.parse_option(sams, 'S', "sams", "File with SAM/BAM/CRAMs to analyze (one per line).");

  parser.parse_option(
    sam_index_fn,
    'i',
    "sams_index",
    "File containing a list of BAM/CRAM indices to use when BAM/CRAM files are queried (one per line). "
    "The index files must be given in the same order as the BAMs/CRAMs.");

  // parser.parse_option(opts.vcf, ' ', "vcf", "Input VCF file with SNP/indel variant sites.");

  parser.parse_positional_argument(ref_fn, "REF.FA", "Reference genome in FASTA format.");
  parser.parse_positional_argument(hla_vcf, "vcf", "Input VCF file with known HLA variants.");

  if (see_advanced_options)
    parser.see_advanced_options(true);

  parser.parse_advanced_option(opts.get_sample_names_from_filename,
                               ' ',
                               "get_sample_names_from_filename",
                               "Set to ignore sample names from RG tags and attempt to retrive the sample names "
                               "from filenames instead.");

  parser.parse_advanced_option(force_copy_reference,
                               ' ',
                               "force_copy_reference",
                               "Force copy of the reference FASTA to temporary folder.");

  parser.parse_advanced_option(force_no_copy_reference,
                               ' ',
                               "force_no_copy_reference",
                               "Force that the reference FASTA is NOT copied to temporary folder.");

  parser.parse_advanced_option(opts.no_filter_on_coverage,
                               ' ',
                               "no_filter_on_coverage",
                               "Set to disable filter on coverage.");

  // parser.parse_advanced_option(opts.force_no_filter_zero_qual,
  //                             ' ',
  //                             "force_no_filter_zero_qual",
  //                             "Set to force variants to be in the final output despite they have zero quality"
  //                             " (all calls are ref/ref)");

  parser.parse_advanced_option(opts.force_use_input_ref_for_cram_reading,
                               ' ',
                               "force_use_input_ref_for_cram_reading",
                               "Force using the input reference FASTA file when reading CRAMs.");

  parser.parse_advanced_option(opts.no_cleanup,
                               ' ',
                               "no_cleanup",
                               "Set to skip removing temporary files. Useful for debugging.");

  parser.parse_advanced_option(opts.is_csi, 'C', "csi", "If set, graphtyper will make csi indices instead of tbi.");

  parser.parse_advanced_option(opts.max_files_open,
                               ' ',
                               "max_files_open",
                               "Select how many files can be open at the same time.");

  parser.parse_advanced_option(
    opts.sam_flag_filter,
    ' ',
    "sam_flag_filter",
    "ANY of these bits set in reads will be completely ignored by GraphTyper. Use with care.");

  parser.parse_advanced_option(
    opts.primer_bedpe,
    ' ',
    "primer_bedpe",
    "In case you have PCR amplicons sequence data you may specify here a BEDPE file containing "
    "coordinates of the primer pairs. The primer sequence will be used in GraphTyper's alignment but are "
    "ignored in the genotyping model.");

  parser.parse_advanced_option(opts.force_ignore_segment,
                               ' ',
                               "force_ignore_segment",
                               "Set to ignore segment calling. Just print graph variant sites instead.");

  parser.parse_advanced_option(opts.uncompressed_sample_names,
                               ' ',
                               "uncompressed_sample_names",
                               "If set, graphtyper BGZF VCF (.vcf.gz) files will output samples names in uncompressed "
                               "(0-level) in the header. The byte range of these blocks will also be printed in "
                               "${prefix}.samples_byte_range. Requires external utils: awk, truncate, stat.");

  opts.no_decompose = true;
  opts.is_segment_calling = true;
  opts.force_no_filter_zero_qual = true;
  opts.force_no_filter_bad_alts = true;
  opts.output_all_variants = true;
  opts.split_var_threshold = 2;
  opts.hq_reads = true;
  opts.is_one_genotype_per_haplotype = true;

  parser.finalize();
  setup_logger();

  print_log(gyper::log_severity::info, "Running the 'genotype_hla' subcommand.");

#ifndef NDEBUG
  // Create stats directory if it doesn't exist
  if (opts.stats.size() > 0 && !gyper::is_directory(opts.stats))
    mkdir(opts.stats.c_str(), 0755);
#endif // NDEBUG

  // Get the genomic regions to process from the --region and --region_file options
  std::vector<gyper::GenomicRegion> regions = get_regions(ref_fn, opts_region, opts_region_file, 1000000);

  // Get the BAM/CRAM file names
  std::vector<std::string> sams_fn = get_sams(sam, sams);
  std::vector<double> avg_cov_by_readlen = get_avg_cov_by_readlen(avg_cov_by_readlen_fn, sams_fn.size());
  std::vector<std::string> sam_index_paths = get_sam_index_paths(sam_index_fn, sams_fn);

  // If neither force copy reference or force no copy reference we determine it from number of SAMs
  bool is_copy_reference = force_copy_reference || (!force_no_copy_reference && sams_fn.size() >= 100);

  gyper::genotype_hla_regions(ref_fn,
                              hla_vcf,
                              interval_fn,
                              sams_fn,
                              sam_index_paths,
                              avg_cov_by_readlen,
                              regions,
                              output_dir,
                              is_copy_reference);

  return 0;
}

int subcmd_genotype_sv(paw::Parser & parser)
{
  gyper::Options & opts = *(gyper::Options::instance());

  std::string avg_cov_by_readlen_fn;
  std::string output_dir{"sv_results"};
  std::string opts_region{};
  std::string opts_region_file{};
  std::string ref_fn{};
  std::string sam{};
  std::string sams{};
  std::string sv_vcf{};
  std::string encoding{"vcf"};
  bool force_copy_reference{false};
  bool force_no_copy_reference{false};
  bool see_advanced_options{false};
  bool force_filter_zero_qual{false};

  // Change the default values
  opts.max_files_open = 128;

  // Parse options
  parser.parse_option(see_advanced_options,
                      'a',
                      "advanced",
                      "Set to enable advanced options. "
                      "See a list of all options (including advanced) with 'graphtyper genotype_sv --advanced --help'");

  parser.parse_option(
    avg_cov_by_readlen_fn,
    ' ',
    "avg_cov_by_readlen",
    "File with average coverage by read length (one value per line). "
    "The values are used for subsampling regions with extremely high coverage and should be in the same "
    "order as the BAM/CRAM list.");

  parser.parse_option(output_dir, 'O', "output", "Output directory.");
  parser.parse_option(opts_region, 'r', "region", "Genomic region to genotype.");
  parser.parse_option(opts_region_file, 'R', "region_file", "File with genomic regions to genotype.");
  parser.parse_option(opts.threads, 't', "threads", "Max. number of threads to use.");
  parser.parse_option(sam, 's', "sam", "SAM/BAM/CRAM to analyze.");
  parser.parse_option(sams, 'S', "sams", "File with SAM/BAM/CRAMs to analyze (one per line).");

  // parser.parse_option(opts.vcf, ' ', "vcf", "Input VCF file with SNP/indel variant sites.");

  parser.parse_positional_argument(ref_fn, "REF.FA", "Reference genome in FASTA format.");
  parser.parse_positional_argument(sv_vcf,
                                   "vcf",
                                   "Input VCF file with structural variant sites and optionally also SNP/indel sites.");

  if (see_advanced_options)
    parser.see_advanced_options(true);

#ifndef NDEBUG
  parser.parse_advanced_option(opts.stats, ' ', "stats", "Directory for statistics files.");
#endif // NDEBUG

  parser.parse_advanced_option(opts.get_sample_names_from_filename,
                               ' ',
                               "get_sample_names_from_filename",
                               "Set to ignore sample names from RG tags and attempt to retrive the sample names "
                               "from filenames instead.");

  parser.parse_advanced_option(force_copy_reference,
                               ' ',
                               "force_copy_reference",
                               "Force copy of the reference FASTA to temporary folder.");

  parser.parse_advanced_option(force_no_copy_reference,
                               ' ',
                               "force_no_copy_reference",
                               "Force that the reference FASTA is NOT copied to temporary folder.");

  parser.parse_advanced_option(opts.no_filter_on_coverage,
                               ' ',
                               "no_filter_on_coverage",
                               "Set to disable filter on coverage.");

  parser.parse_advanced_option(opts.force_no_filter_zero_qual,
                               ' ',
                               "force_no_filter_zero_qual",
                               "Set to force variants to be in the final output despite they have zero quality"
                               " (all calls are ref/ref). NOTE: This is now the default behaviour, use "
                               "'--force_filter_zero_qual' to turn it off.");

  parser.parse_advanced_option(force_filter_zero_qual,
                               ' ',
                               "force_filter_zero_qual",
                               "Set to filter variants to be in the final output despite they have zero quality"
                               " (all calls are ref/ref)");

  parser.parse_advanced_option(opts.force_use_input_ref_for_cram_reading,
                               ' ',
                               "force_use_input_ref_for_cram_reading",
                               "Force using the input reference FASTA file when reading CRAMs.");

  parser.parse_advanced_option(opts.no_cleanup,
                               ' ',
                               "no_cleanup",
                               "Set to skip removing temporary files. Useful for debugging.");

  parser.parse_advanced_option(opts.is_csi, 'C', "csi", "If set, graphtyper will make csi indices instead of tbi.");

  parser.parse_advanced_option(opts.max_files_open,
                               ' ',
                               "max_files_open",
                               "Select how many files can be open at the same time.");

  parser.parse_advanced_option(
    opts.sam_flag_filter,
    ' ',
    "sam_flag_filter",
    "ANY of these bits set in reads will be completely ignored by GraphTyper. Use with care.");

  parser.parse_advanced_option(
    opts.primer_bedpe,
    ' ',
    "primer_bedpe",
    "In case you have PCR amplicons sequence data you may specify here a BEDPE file containing "
    "coordinates of the primer pairs. The primer sequence will be used in GraphTyper's alignment but are "
    "ignored in the genotyping model.");

  parser.parse_advanced_option(opts.uncompressed_sample_names,
                               ' ',
                               "uncompressed_sample_names",
                               "If set, graphtyper BGZF VCF (.vcf.gz) files will output samples names in uncompressed "
                               "(0-level) in the header. The byte range of these blocks will also be printed in "
                               "${prefix}.samples_byte_range.");

  parser.parse_advanced_option(encoding, ' ', "encoding", "Select output encoding. Available are: vcf, popvcf");

  // set default compression level as 9 when popvcf, since level 9 is already very fast anyway in that encoding mode
  if (encoding == "popvcf")
    opts.bgzf_compression_level = 9;

  // Changed behaviour such that zero qual SVs are not filtered out by default
  if (!force_filter_zero_qual)
    opts.force_no_filter_zero_qual = true;

  if (opts.force_no_filter_zero_qual)
    opts.force_no_filter_bad_alts = true;

  parser.finalize();
  setup_logger();

  opts.encoding = (encoding == "popvcf") ? 'p' : 'v';
  print_log(gyper::log_severity::info, "Running the 'genotype_sv' subcommand.");

#ifndef NDEBUG
  // Create stats directory if it doesn't exist
  if (opts.stats.size() > 0 && !gyper::is_directory(opts.stats))
    mkdir(opts.stats.c_str(), 0755);
#endif // NDEBUG

  // Get the genomic regions to process from the --region and --region_file options
  std::vector<gyper::GenomicRegion> regions = get_regions(ref_fn, opts_region, opts_region_file, 1000000);

  // Get the BAM/CRAM file names
  std::vector<std::string> sams_fn = get_sams(sam, sams);
  std::vector<double> avg_cov_by_readlen = get_avg_cov_by_readlen(avg_cov_by_readlen_fn, sams_fn.size());

  // If neither force copy reference or force no copy reference we determine it from number of SAMs
  bool is_copy_reference = force_copy_reference || (!force_no_copy_reference && sams_fn.size() >= 100);

  gyper::genotype_sv_regions(ref_fn, sv_vcf, sams_fn, avg_cov_by_readlen, regions, output_dir, is_copy_reference);

  return 0;
}

int subcmd_genotype_lr(paw::Parser & parser)
{
  gyper::Options & opts = *(gyper::Options::instance());

  std::string output_dir{"sv_results"};
  std::string opts_region{};
  std::string opts_region_file{};
  std::string ref_fn{};
  std::string sam{};
  std::string sams{};
  bool force_copy_reference{false};
  bool force_no_copy_reference{false};
  bool see_advanced_options{false};

  // Change the default values
  opts.max_files_open = 128;
  opts.is_lr_calling = true;

  // Parse options
  parser.parse_option(see_advanced_options,
                      'a',
                      "advanced",
                      "Set to enable advanced options. "
                      "See a list of all options (including advanced) with 'graphtyper genotype_sv --advanced --help'");

  parser.parse_option(output_dir, 'O', "output", "Output directory.");
  parser.parse_option(opts_region, 'r', "region", "Genomic region to genotype.");
  parser.parse_option(opts_region_file, 'R', "region_file", "File with genomic regions to genotype.");
  parser.parse_option(opts.threads, 't', "threads", "Max. number of threads to use.");
  parser.parse_option(sam, 's', "sam", "SAM/BAM/CRAM to analyze.");
  parser.parse_option(sams, 'S', "sams", "File with SAM/BAM/CRAMs to analyze (one per line).");

  parser.parse_positional_argument(ref_fn, "REF.FA", "Reference genome in FASTA format.");

  if (see_advanced_options)
    parser.see_advanced_options(true);

#ifndef NDEBUG
  parser.parse_advanced_option(opts.stats, ' ', "stats", "Directory for statistics files.");
#endif // NDEBUG

  parser.parse_advanced_option(opts.get_sample_names_from_filename,
                               ' ',
                               "get_sample_names_from_filename",
                               "Set to ignore sample names from RG tags and attempt to retrive the sample names "
                               "from filenames instead.");

  parser.parse_advanced_option(force_copy_reference,
                               ' ',
                               "force_copy_reference",
                               "Force copy of the reference FASTA to temporary folder.");

  parser.parse_advanced_option(force_no_copy_reference,
                               ' ',
                               "force_no_copy_reference",
                               "Force that the reference FASTA is NOT copied to temporary folder.");

  parser.parse_advanced_option(opts.no_filter_on_coverage,
                               ' ',
                               "no_filter_on_coverage",
                               "Set to disable filter on coverage.");

  parser.parse_advanced_option(opts.force_no_filter_zero_qual,
                               ' ',
                               "force_no_filter_zero_qual",
                               "Set to force variants to be in the final output despite they have zero quality"
                               " (all calls are ref/ref)");

  parser.parse_advanced_option(opts.force_use_input_ref_for_cram_reading,
                               ' ',
                               "force_use_input_ref_for_cram_reading",
                               "Force using the input reference FASTA file when reading CRAMs.");

  parser.parse_advanced_option(opts.no_cleanup,
                               ' ',
                               "no_cleanup",
                               "Set to skip removing temporary files. Useful for debugging.");

  parser.parse_advanced_option(opts.is_csi, 'C', "csi", "If set, graphtyper will make csi indices instead of tbi.");

  parser.parse_advanced_option(opts.max_files_open,
                               ' ',
                               "max_files_open",
                               "Select how many files can be open at the same time.");

  parser.parse_advanced_option(
    opts.sam_flag_filter,
    ' ',
    "sam_flag_filter",
    "ANY of these bits set in reads will be completely ignored by GraphTyper. Use with care.");

  parser.parse_advanced_option(opts.lr_mapq_filter,
                               ' ',
                               "lr_mapq_filter",
                               "Filter reads with MAPQ below this threshold.");

  parser.parse_advanced_option(opts.lr_coverage_filter,
                               ' ',
                               "lr_coverage_filter",
                               "Filter reads after seeing coverage above this value. Set as -1 to remove the filter.");

  parser.parse_advanced_option(opts.uncompressed_sample_names,
                               ' ',
                               "uncompressed_sample_names",
                               "If set, graphtyper BGZF VCF (.vcf.gz) files will output samples names in uncompressed "
                               "(0-level) in the header. The byte range of these blocks will also be printed in "
                               "${prefix}.samples_byte_range. Requires external utils: awk, truncate, stat.");

  if (opts.force_no_filter_zero_qual)
    opts.force_no_filter_bad_alts = true;

  parser.finalize();
  setup_logger();

  print_log(gyper::log_severity::info, "Running the 'genotype_lr' subcommand.");

#ifndef NDEBUG
  // Create stats directory if it doesn't exist
  if (opts.stats.size() > 0 && !gyper::is_directory(opts.stats))
    mkdir(opts.stats.c_str(), 0755);
#endif // NDEBUG

  // Get the genomic regions to process from the --region and --region_file options
  std::vector<gyper::GenomicRegion> regions = get_regions(ref_fn, opts_region, opts_region_file, 50000);

  // Get the BAM/CRAM file names
  std::vector<std::string> sams_fn = get_sams(sam, sams);

  // If neither force copy reference or force no copy reference we determine it from number of SAMs
  bool is_copy_reference = force_copy_reference || (!force_no_copy_reference && sams_fn.size() >= 100);

  opts.filter_on_read_bias = false;
  opts.filter_on_proper_pairs = false;

  gyper::genotype_lr_regions(ref_fn, sams_fn, regions, output_dir, is_copy_reference);

  return 0;
}

int subcmd_genotype_camou(paw::Parser & parser)
{
  gyper::Options & opts = *(gyper::Options::instance());

  std::string avg_cov_by_readlen_fn;
  std::string interval_fn;
  std::string output_dir = "results";
  std::string ref_fn;
  std::string sam;
  std::string sams;

  parser.parse_option(avg_cov_by_readlen_fn, ' ', "avg_cov_by_readlen", "File with average coverage by read length.");

  parser.parse_option(opts.max_files_open,
                      ' ',
                      "max_files_open",
                      "Select how many files can be open at the same time.");

  parser.parse_option(opts.no_bamshrink, ' ', "no_bamshrink", "Set to skip bamShrink.");

  parser.parse_option(opts.no_cleanup,
                      ' ',
                      "no_cleanup",
                      "Set to skip removing temporary files. Useful for debugging.");

  parser.parse_option(output_dir,
                      'O',
                      "output",
                      "Output directory. Results will be written in <output>/<contig>/<region>.vcf.gz");

  parser.parse_option(sam, 's', "sam", "SAM/BAM/CRAM to analyze.");
  parser.parse_option(sams, 'S', "sams", "File with SAM/BAM/CRAMs to analyze (one per line).");
  parser.parse_option(opts.threads, 't', "threads", "Max. number of threads to use.");
  parser.parse_option(opts.vcf,
                      ' ',
                      "vcf",
                      "Input VCF file with variant sites. "
                      "Use this option if you want GraphTyper to only genotype variants from this VCF.");

  parser.parse_advanced_option(opts.uncompressed_sample_names,
                               ' ',
                               "uncompressed_sample_names",
                               "If set, graphtyper BGZF VCF (.vcf.gz) files will output samples names in uncompressed "
                               "(0-level) in the header. The byte range of these blocks will also be printed in "
                               "${prefix}.samples_byte_range. Requires external utils: awk, truncate, stat.");

  parser.parse_positional_argument(ref_fn, "REF.FA", "Reference genome in FASTA format.");
  parser.parse_positional_argument(interval_fn,
                                   "interval-file",
                                   "3-column BED type file with interval(s)/region(s) to filter on.");

  parser.finalize();
  setup_logger();

  // Force no filter on MAPQ in camou calling and read alignment in both orientations
  opts.bamshrink_is_not_filtering_mapq0 = true;
  opts.filter_on_mapq = false;
  opts.force_align_both_orientations = true;
  opts.force_no_filter_bad_alts = true;

  // Get the SAM/BAM/CRAM file names
  std::vector<std::string> sams_fn = get_sams(sam, sams);

  // Get the avgCovByReadLen for each of the SAM/BAM/CRAM
  std::vector<double> avg_cov_by_readlen = get_avg_cov_by_readlen(avg_cov_by_readlen_fn, sams_fn.size());

  gyper::genotype_camou(interval_fn, ref_fn, sams_fn, output_dir, avg_cov_by_readlen);

  return 0;
}

int subcmd_vcf_break_down(paw::Parser & parser)
{
  std::string graph_fn;
  std::string output_fn = "-";
  std::string region;
  std::string vcf_fn;

  parser.parse_option(output_fn, 'o', "output", "Output VCF file name.");
  parser.parse_option(region, 'r', "region", "Region to print variant in.");

  parser.parse_positional_argument(graph_fn, "GRAPH", "Path to graph.");
  parser.parse_positional_argument(vcf_fn, "VCF", "Path to VCF file to break down.");

  parser.finalize();
  setup_logger();

  // load the graph
  gyper::load_graph(graph_fn);
  gyper::vcf_break_down(vcf_fn, output_fn, region);
  return 0;
}

int subcmd_vcf_concatenate(paw::Parser & parser)
{
  std::vector<std::string> vcfs;
  bool no_sort{false};
  std::string output_fn{"-"};
  bool sites_only{false};
  std::string region{};
  bool write_tbi{false};
  bool is_sv_vcf{false};

  parser.parse_option(no_sort, ' ', "no_sort", "Set to skip sorting the variants.");
  parser.parse_option(output_fn, 'o', "output", "Output VCF file name.");
  parser.parse_option(sites_only, ' ', "sites_only", "Set to write only variant site information.");
  parser.parse_option(is_sv_vcf, ' ', "sv", "Set if the input VCFs were generated from genotype_sv.");
  parser.parse_option(region, 'r', "region", "Region to print variant in.");
  parser.parse_option(write_tbi, 't', "write_tbi", "Set to write TBI index.");

  parser.parse_remaining_positional_arguments(vcfs, "vcfs...", "VCFs to concatenate");

  parser.finalize();
  setup_logger();

  if (is_sv_vcf)
    gyper::graph.is_sv_graph = true;

  gyper::vcf_concatenate(vcfs, output_fn, no_sort, sites_only, write_tbi, region);
  return 0;
}

int subcmd_vcf_merge(paw::Parser & parser)
{
  gyper::Options & opts = *(gyper::Options::instance());

  std::vector<std::string> vcfs;
  std::string output_fn;
  std::string file_list;
  std::string encoding{"vcf"};
  bool is_sv_vcf{false};

  parser.parse_option(output_fn, 'o', "output", "Output VCF file name.");
  parser.parse_option(file_list, ' ', "file_list", "File containing VCFs to merge.");
  parser.parse_option(is_sv_vcf, ' ', "sv", "Set if the input VCFs were generated from genotype_sv.");
  parser.parse_option(encoding, ' ', "encoding", "Select output encoding. Available are: vcf, popvcf");

  // set default compression level as 9 when popvcf, since level 9 is already very fast anyway in that encoding mode
  if (encoding == "popvcf")
    opts.bgzf_compression_level = 9;

  parser.parse_remaining_positional_arguments(vcfs, "vcfs...", "VCFs to merge");

  parser.finalize();
  setup_logger();

  opts.encoding = (encoding == "popvcf") ? 'p' : 'v';

  if (is_sv_vcf)
    gyper::graph.is_sv_graph = true;

  if (file_list.size() > 0)
  {
    std::ifstream files;
    files.open(file_list);

    if (!files.is_open())
    {
      print_log(gyper::log_severity::error, "Could not open file '", file_list, "'");
      return 1;
    }

    for (std::string line; std::getline(files, line);)
      vcfs.push_back(line);
  }

  gyper::vcf_merge(vcfs, output_fn);
  return 0;
}

} // namespace

int main(int argc, char ** argv)
{
  int ret{0};
  paw::Parser parser(argc, argv);
  parser.set_name("GraphTyper");
  parser.set_version(graphtyper_VERSION_MAJOR, graphtyper_VERSION_MINOR, graphtyper_VERSION_PATCH);

  try
  {
    std::string subcmd{};

    parser.add_subcommand("bamshrink", "Run bamShrink.");
    parser.add_subcommand("call", "Call variants of a graph.");
    parser.add_subcommand("check", "Check a GraphTyper graph (useful for debugging).");
    parser.add_subcommand("construct", "Construct a graph.");
    parser.add_subcommand("genotype", "Run the SNP/indel genotyping pipeline.");
    parser.add_subcommand("genotype_camou", "(WIP) Run the camou SNP/indel genotyping pipeline.");
    parser.add_subcommand("genotype_hla", "(WIP) Run the HLA genotyping pipeline.");
    parser.add_subcommand("genotype_lr", "(WIP) Run the camou LR genotyping pipeline.");
    parser.add_subcommand("genotype_sv", "Run the structural variant (SV) genotyping pipeline.");
    parser.add_subcommand("index", "(deprecated) Index a graph.");
    parser.add_subcommand("vcf_break_down", "Break down/decompose a VCF file.");
    parser.add_subcommand("vcf_concatenate", "Concatenate VCF files.");
    parser.add_subcommand("vcf_merge", "Merge VCF files.");
    parser.parse_subcommand(subcmd);

    gyper::Options & opts = *(gyper::Options::instance());

    parser.parse_option(opts.log, 'l', "log", "Set path to log file.");
    parser.parse_option(opts.verbose, 'v', "verbose", "Set to output verbose logging.");
    parser.parse_option(opts.vverbose, ' ', "vverbose", "Set to output very verbose logging.");

    if (subcmd == "bamshrink")
      ret = subcmd_bamshrink(parser);
    else if (subcmd == "check")
      ret = subcmd_check(parser);
    else if (subcmd == "construct")
      ret = subcmd_construct(parser);
    else if (subcmd == "genotype")
      ret = subcmd_genotype(parser);
    else if (subcmd == "genotype_camou")
      ret = subcmd_genotype_camou(parser);
    else if (subcmd == "genotype_lr")
      ret = subcmd_genotype_lr(parser);
    else if (subcmd == "genotype_hla")
      ret = subcmd_genotype_hla(parser);
    else if (subcmd == "genotype_sv")
      ret = subcmd_genotype_sv(parser);
    else if (subcmd == "index")
    {
      std::cerr << "GraphTyper's 'index' subcommand deprecated. "
                << "Graph is now indexed automatically after construction.\n";
      std::exit(1);
    }
    else if (subcmd == "vcf_break_down")
      ret = subcmd_vcf_break_down(parser);
    else if (subcmd == "vcf_concatenate")
      ret = subcmd_vcf_concatenate(parser);
    else if (subcmd == "vcf_merge")
      ret = subcmd_vcf_merge(parser);
    else if (subcmd.size() == 0)
    {
      parser.finalize();
    }
    else
    {
      parser.finalize();
      std::exit(1);
    }
  }
  catch (paw::exception::help const & e)
  {
    std::cout << e.what();
    return 0;
  }
  catch (std::exception const & e)
  {
    std::cerr << e.what();
    return 1;
  }

  return ret;
}
