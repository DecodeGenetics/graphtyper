#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <thread>
#include <unistd.h>

#include "args.hxx"

#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/console.hpp>

#include <graphtyper/constants.hpp>
#include <graphtyper/graph/sequence_extractor.hpp>
#include <graphtyper/graph/haplotype_extractor.hpp>
#include <graphtyper/graph/constructor.hpp>
#include <graphtyper/graph/graph_serialization.hpp>
#include <graphtyper/index/rocksdb.hpp>
#include <graphtyper/index/indexer.hpp>
#include <graphtyper/typer/aligner.hpp>
#include <graphtyper/typer/caller.hpp>
#include <graphtyper/typer/discovery.hpp>
#include <graphtyper/typer/vcf_operations.hpp>
#include <graphtyper/utilities/adapter_remover.hpp>
#include <graphtyper/utilities/options.hpp>


namespace
{

bool
is_file(std::string filename)
{
  struct stat sb;
  return stat(filename.c_str(), &sb) == 0 && (S_ISREG(sb.st_mode) || S_ISLNK(sb.st_mode));
}


bool
is_directory(std::string filename)
{
  struct stat sb;
  return stat(filename.c_str(), &sb) == 0;
}


void
print_default_help()
{
  // Print this help to see this (not so) awesome ASCII art!
  std::cerr << "-. .-.   .-. .-.   .-. .-.   .-. .-.   .-. .-.   .-. .-\n"
            << "||\\|||\\ /|R|\\|||\\ /|P|\\|||\\ /|T|\\|||\\ /|P|\\|||\\ /|R||\\|\n"
            << "|/ \\|G|\\|||/ \\|A|\\|||/ \\|H|\\|||/ \\|Y|\\|||/ \\|E|\\|||/ \\|\n"
            << "~   `-~ `-`   `-~ `-`   `-~ `-~   `-~ `-`   `-~ `-`   `\n\n"
            << "Graphtyper is a program which aligns and genotype calls sequenced reads using an acyclic graph structure.\n"
            << "\n"
            << "Version: " << graphtyper_VERSION_MAJOR << "." << graphtyper_VERSION_MINOR;

  if (std::string(GIT_NUM_DIRTY_FILES) != std::string("0"))
    std::cerr << "-dirty";

  std::cerr << " (" << GIT_COMMIT_SHORT_HASH << ")\n\n"
            << "Usage: graphtyper <COMMAND> [OPTIONS...]\n"
            << '\n'
            << "Available commands are:\n"
            << "  call             Genotype calls sample(s).\n"
            << "  construct        Construct a graph.\n"
            << "  haplotypes       Extracts called haplotypes into a VCF file.\n"
            << "  index            Indexes a graph.\n"
            << "  vcf_break_down   Breaks down variants in Graphtyper VCF files.\n"
            << "  vcf_merge        Merged Graphtyper genotype calls VCF files.\n"
            << "  vcf_concatenate  Concatenates Graphtyper VCF files.\n"
            << "\n"
            << "For more information about each command enter\n"
            << "  graphtyper <COMMAND> --help\n"
            << "\n"
            << "Author: Hannes Pétur Eggertsson (Hannes.Eggertsson@decode.is)" << std::endl;

}


std::unique_ptr<args::HelpFlag>
add_arg_help(args::ArgumentParser & parser)
{
  return std::unique_ptr<args::HelpFlag>(new args::HelpFlag(parser, "help", "Display this help.", {'h', "help"}));
}


/** Command argument */
using TCommand = args::Positional<std::string>;

std::unique_ptr<TCommand>
add_arg_command(args::ArgumentParser & parser, std::string const & metavar = "COMMAND")
{
  return std::unique_ptr<TCommand>(new TCommand(parser, metavar, "Graphtyper command to execute."));
}


/** Graph argument */
using TGraph = args::Positional<std::string>;

std::unique_ptr<TGraph>
add_arg_graph(args::ArgumentParser & parser)
{
  return std::unique_ptr<TGraph>(new TGraph(parser, "GRAPH", "Graph file to use."));
}


/** Index argument */
using TIndex = args::ValueFlag<std::string>;

std::unique_ptr<TIndex>
add_arg_index(args::ArgumentParser & parser)
{
  return std::unique_ptr<TIndex>(new TIndex(parser, "DIR", "Graphtyper index directory.", {"index"}));
}


/** FASTA argument */
using TFasta = args::Positional<std::string>;

std::unique_ptr<TFasta>
add_arg_fasta(args::ArgumentParser & parser)
{
  return std::unique_ptr<TFasta>(new TFasta(parser, "REFERENCE.fa", "Reference FASTA file to use.", {"fasta"}));
}


/** VCF argument */
using TVcfVal = args::ValueFlag<std::string>;

std::unique_ptr<TVcfVal>
add_arg_vcf(args::ArgumentParser & parser)
{
  return std::unique_ptr<TVcfVal>(new TVcfVal(parser, "FILE.vcf.gz", "VCF bgzipped file with variants.", {"vcf"}));
}


using TVcfPos = args::Positional<std::string>;

std::unique_ptr<TVcfPos>
add_arg_vcf_pos(args::ArgumentParser & parser)
{
  return std::unique_ptr<TVcfPos>(new TVcfPos(parser, "FILE.vcf", "VCF file.", {"vcf"}));
}


/** Regions argument */
using TRegions = args::PositionalList<std::string>;

std::unique_ptr<TRegions>
add_arg_regions(args::ArgumentParser & parser)
{
  return std::unique_ptr<TRegions>(new TRegions(parser, "REGIONS", "Regions to use."));
}


using TRegion = args::Positional<std::string>;

std::unique_ptr<TRegion>
add_arg_region(args::ArgumentParser & parser)
{
  return std::unique_ptr<TRegion>(new TRegion(parser, "REGION", "Region to use."));
}


using TRegionVal = args::ValueFlag<std::string>;

std::unique_ptr<TRegionVal>
add_arg_region_val(args::ArgumentParser & parser)
{
  return std::unique_ptr<TRegionVal>(new TRegionVal(parser, "REGION", "Region to use.", {"r", "region"}, "."));
}


/** LOG argument */
using TLog = args::ValueFlag<std::string>;

std::unique_ptr<TLog>
add_arg_log(args::ArgumentParser & parser)
{
  return std::unique_ptr<TLog>(new TLog(parser, "LOG.txt", "Log filename. If none specified logs are written to standard output.", {"log"}));
}


void
parse_log(TLog & log_arg)
{
  if (log_arg)
  {
    boost::log::core::get()->set_filter
    (
      boost::log::trivial::severity >= boost::log::trivial::debug
    );

    const std::string log_filename = args::get(log_arg);

    if (log_filename == "-")
    {
      boost::log::add_console_log(std::clog, boost::log::keywords::auto_flush = false);
      return;
    }

    // Create a sink
    gyper::Options::instance()->sink =
      boost::log::add_file_log
      (
        boost::log::keywords::file_name = log_filename.c_str(),
        boost::log::keywords::auto_flush = false
      );
  }
  else
  {
    boost::log::add_console_log(std::clog, boost::log::keywords::auto_flush = false);
  }
}


/** Output argument */
using TOutput = args::ValueFlag<std::string>;

std::unique_ptr<TOutput>
add_arg_output(args::ArgumentParser & parser)
{
  return std::unique_ptr<TOutput>(new TOutput(parser, "FILE", "Output file.", {"o", "output"}));
}


std::unique_ptr<TOutput>
add_arg_output_dir(args::ArgumentParser & parser)
{
  return std::unique_ptr<TOutput>(new TOutput(parser, "DIR", "Output directory.", {"o", "output"}, "."));
}


/** SAM argument */
using TSam = args::ValueFlag<std::string>;

std::unique_ptr<TSam>
add_arg_sam(args::ArgumentParser & parser)
{
  return std::unique_ptr<TSam>(new TSam(parser, "file.sam", "A SAM/BAM/CRAM file to read reads from.", {"s", "sam"}));
}


std::unique_ptr<TSam>
add_arg_sams(args::ArgumentParser & parser)
{
  return std::unique_ptr<TSam>(new TSam(parser, "SAMs", "A file with a list of SAM/BAM/CRAM files seperated by newlines.", {"S", "sams"}));
}

/** SAM argument */
using TSegment = args::ValueFlag<std::string>;

std::unique_ptr<TSegment>
add_arg_segment(args::ArgumentParser & parser)
{
  return std::unique_ptr<TSegment>(new TSegment(parser, "segment.fa", "FASTA file with segments.", {"e", "segment"}));
}


/** Threads argument */
using TThreads = args::ValueFlag<unsigned>;

std::unique_ptr<TThreads>
add_arg_threads(args::ArgumentParser & parser)
{
  return std::unique_ptr<TThreads>(new TThreads(parser, "N", "Number of threads to use. Use '0' to use automatic detect hardware concurrency.", {"threads"}));
}


void
parse_threads(TThreads & threads_arg)
{
  using gyper::Options;

  if (threads_arg)
  {
    unsigned const threads = args::get(threads_arg);

    if (threads == 0)
      Options::instance()->threads = std::thread::hardware_concurrency();
    else
      Options::instance()->threads = threads;
  }
}


/** VCFs argument */
using TVcfs = args::PositionalList<std::string>;

std::unique_ptr<TVcfs>
add_arg_vcfs(args::ArgumentParser & parser)
{
  return std::unique_ptr<TVcfs>(new TVcfs(parser, "VCFs", "A list of bgzipped VCF files seperated by newlines.", {"V", "vcfs"}));
}

/** max_index_labels argument */
using TMaxIndexLabels = args::ValueFlag<unsigned>;

std::unique_ptr<TMaxIndexLabels>
add_arg_max_index_labels(args::ArgumentParser & parser)
{
  return std::unique_ptr<TMaxIndexLabels>(
    new TMaxIndexLabels(
      parser, "N", "Maximum number labels a single k-mer can be associated with.", {"max_index_labels"}
    )
  );
}

void
parse_max_index_labels(TMaxIndexLabels & max_index_labels)
{
  if (max_index_labels)
    gyper::Options::instance()->max_index_labels = args::get(max_index_labels);
}


/** max_merge_variant_dist argument */
using TMaxMergeVariantDist = args::ValueFlag<unsigned>;

std::unique_ptr<TMaxMergeVariantDist>
add_arg_mmvd(args::ArgumentParser & parser)
{
  return std::unique_ptr<TMaxMergeVariantDist>(new TMaxMergeVariantDist(parser, "N", "Distance to merge variant into haplotypes.", {"max_merge_variant_dist"}));
}


void
parse_mmvd(TMaxMergeVariantDist & mmvd_arg)
{
  if (mmvd_arg)
    gyper::Options::instance()->max_merge_variant_dist = args::get(mmvd_arg);
}


/** stats argument */
using TStats = args::ValueFlag<std::string>;
std::unique_ptr<TStats>
add_arg_stats(args::ArgumentParser & parser)
{
  return std::unique_ptr<TStats>(new TStats(parser, "DIR", "Output statistics to directory.", {"stats"}));
}


void
parse_stats(TStats & stats_arg)
{
  if (stats_arg && args::get(stats_arg).size() > 0)
  {
    gyper::Options::instance()->stats = args::get(stats_arg);

    if (!is_directory(args::get(stats_arg)))
      mkdir(args::get(stats_arg).c_str(), 0755);
  }
}

/** No new variants argument */
std::unique_ptr<args::Flag>
add_arg_no_new_variants(args::ArgumentParser & parser)
{
  return std::unique_ptr<args::Flag>(new args::Flag(parser, "NO_NEW_VARIANTS", "Set to skip check for new variants.", {"no_new_variants"}));
}


void
parse_no_new_variants(args::Flag & no_new_variants_arg)
{
  if (no_new_variants_arg)
    gyper::Options::instance()->no_new_variants = true;
}


/** get_sample_names_from_filename */
std::unique_ptr<args::Flag>
add_arg_get_sample_names_from_filename(args::ArgumentParser & parser)
{
  return std::unique_ptr<args::Flag>(new args::Flag(parser, "", "Set to get sample names from filenames.", {"get_sample_names_from_filename"}));
}

void
parse_get_sample_names_from_filename(args::Flag & get_sample_names_from_filename_arg)
{
  if (get_sample_names_from_filename_arg)
    gyper::Options::instance()->get_sample_names_from_filename = true;
}

/** Only use HQ reads */
std::unique_ptr<args::Flag>
add_arg_hq_reads(args::ArgumentParser & parser)
{
  return std::unique_ptr<args::Flag>(new args::Flag(parser, "HQ_READS", "Set to use only HQ reads.", {"hq_reads"}));
}


void
parse_hq_reads(args::Flag & hq_reads_arg)
{
  if (hq_reads_arg)
    gyper::Options::instance()->hq_reads = true;
}


/** output_all_variants */
std::unique_ptr<args::Flag>
add_arg_output_all_variants(args::ArgumentParser & parser, std::string cmd = "")
{
  if (cmd == std::string("haplotypes"))
    return std::unique_ptr<args::Flag>(new args::Flag(parser, "OUTPUT_ALL_VARIANTS", "Set to output all haplotypes with enough coverage", {"output_all_variants"}));
  else
    return std::unique_ptr<args::Flag>(new args::Flag(parser, "OUTPUT_ALL_VARIANTS", "Set to output all variants, regardless of their allele number or size", {"output_all_variants"}));
}

void
parse_output_all_variants(args::Flag & arg)
{
  if (arg)
    gyper::Options::instance()->output_all_variants = true;
}


/** Skip VCF sorting */
std::unique_ptr<args::Flag>
add_arg_no_sort(args::ArgumentParser & parser)
{
  return std::unique_ptr<args::Flag>(new args::Flag(parser, "NO_SORT", "Set to skip sorting of VCFs.", {"no_sort"}));
}


/** Minimum variant support argument */
using TMinVarSup = args::ValueFlag<unsigned long>;

std::unique_ptr<TMinVarSup>
add_arg_min_var_sup(args::ArgumentParser & parser)
{
  return std::unique_ptr<TMinVarSup>(
    new TMinVarSup(parser,
                   "N",
                   "Minimum number of variants needed to add it to the output VCF file.",
                   {"minimum_variant_support"}
                   )
    );
}


void
parse_minimum_variant_support(TMinVarSup & minimum_variant_support_arg)
{
  if (minimum_variant_support_arg)
    gyper::Options::instance()->minimum_variant_support = args::get(minimum_variant_support_arg);
}


/** soft_cap_of_variants_in_100_bp_window */
using TSoftCapOfVariantsInWindow = args::ValueFlag<unsigned>;

std::unique_ptr<TSoftCapOfVariantsInWindow>
add_arg_soft_cap_of_variants_in_100_bp_window(args::ArgumentParser & parser)
{
  return std::unique_ptr<TSoftCapOfVariantsInWindow>(
    new TSoftCapOfVariantsInWindow(parser,
                                  "N",
                                  "Soft cap of number of variants that can appear in a 100 bp window.",
                                  {"soft_cap_of_variants_in_100_bp_window"}
                                  )
    );
}


void
parse_soft_cap_of_variants_in_100_bp_window(TSoftCapOfVariantsInWindow & soft_cap_of_variants_in_100_bp_window_arg)
{
  if (soft_cap_of_variants_in_100_bp_window_arg)
    gyper::Options::instance()->soft_cap_of_variants_in_100_bp_window = args::get(soft_cap_of_variants_in_100_bp_window_arg);
}


/** hard_cap_of_variants_in_100_bp_window */
using THardCapOfVariantsInWindow = args::ValueFlag<unsigned>;

std::unique_ptr<THardCapOfVariantsInWindow>
add_arg_hard_cap_of_variants_in_100_bp_window(args::ArgumentParser & parser)
{
  return std::unique_ptr<THardCapOfVariantsInWindow>(
    new THardCapOfVariantsInWindow(parser,
                                  "N",
                                  "Hard cap of number of variants that can appear in a 100 bp window.",
                                  {"hard_cap_of_variants_in_100_bp_window"}
                                  )
    );
}


void
parse_hard_cap_of_variants_in_100_bp_window(TSoftCapOfVariantsInWindow & hard_cap_of_variants_in_100_bp_window_arg)
{
  if (hard_cap_of_variants_in_100_bp_window_arg)
    gyper::Options::instance()->hard_cap_of_variants_in_100_bp_window = args::get(hard_cap_of_variants_in_100_bp_window_arg);
}


/** soft_cap_of_non_snps_in_100_bp_window */
using TSoftCapOfNonSnpsInWindow = args::ValueFlag<unsigned>;

std::unique_ptr<TSoftCapOfNonSnpsInWindow>
add_arg_soft_cap_of_non_snps_in_100_bp_window(args::ArgumentParser & parser)
{
  return std::unique_ptr<TSoftCapOfNonSnpsInWindow>(
    new TSoftCapOfNonSnpsInWindow(parser,
                                  "N",
                                  "Soft cap of number of non-SNP variants that can appear in a 100 bp window.",
                                  {"soft_cap_of_non_snps_in_100_bp_window"}
                                  )
    );
}


void
parse_soft_cap_of_non_snps_in_100_bp_window(TSoftCapOfNonSnpsInWindow & soft_cap_of_non_snps_in_100_bp_window_arg)
{
  if (soft_cap_of_non_snps_in_100_bp_window_arg)
    gyper::Options::instance()->soft_cap_of_non_snps_in_100_bp_window = args::get(soft_cap_of_non_snps_in_100_bp_window_arg);
}


/** hard_cap_of_non_snps_in_100_bp_window */
using THardCapOfNonSnpsInWindow = args::ValueFlag<unsigned>;

std::unique_ptr<THardCapOfNonSnpsInWindow>
add_arg_hard_cap_of_non_snps_in_100_bp_window(args::ArgumentParser & parser)
{
  return std::unique_ptr<THardCapOfNonSnpsInWindow>(
    new THardCapOfNonSnpsInWindow(parser,
                                  "N",
                                  "Hard cap of number of non-SNP variants that can appear in a 100 bp window.",
                                  {"hard_cap_of_non_snps_in_100_bp_window"}
                                  )
    );
}


void
parse_hard_cap_of_non_snps_in_100_bp_window(THardCapOfNonSnpsInWindow & hard_cap_of_non_snps_in_100_bp_window_arg)
{
  if (hard_cap_of_non_snps_in_100_bp_window_arg)
    gyper::Options::instance()->hard_cap_of_non_snps_in_100_bp_window = args::get(hard_cap_of_non_snps_in_100_bp_window_arg);
}


/** Minimum variant support ratio argument */
using TMinVarSupRatio = args::ValueFlag<double>;

std::unique_ptr<TMinVarSupRatio>
add_arg_min_var_sup_ratio(args::ArgumentParser & parser)
{
  return std::unique_ptr<TMinVarSupRatio>(
    new TMinVarSupRatio(parser,
                        "D",
                        "Minimum ratio of variant support by variant coverage to add the variant to the output VCF file.",
                        {"minimum_variant_support_ratio"}
                        )
    );
}


void
parse_minimum_variant_support_ratio(TMinVarSupRatio & minimum_variant_support_ratio_arg)
{
  if (minimum_variant_support_ratio_arg)
    gyper::Options::instance()->minimum_variant_support_ratio = args::get(minimum_variant_support_ratio_arg);
}


/** Maximum homozygous allele balance */
using TMaximumHomozygousAlleleBalance = args::ValueFlag<double>;

std::unique_ptr<TMaximumHomozygousAlleleBalance>
add_arg_maximum_homozygous_allele_balance(args::ArgumentParser & parser)
{
  return std::unique_ptr<TMaximumHomozygousAlleleBalance>(
    new TMaximumHomozygousAlleleBalance(parser,
                        "D",
                        "Maxmimum homozygous allele balance allowed when adding new variants (Default is 0.08)",
                        {"maximum_homozygous_allele_balance"}
                        )
    );
}


void
parse_maximum_homozygous_allele_balance(TMaximumHomozygousAlleleBalance & maximum_homozygous_allele_balance_arg)
{
  if (maximum_homozygous_allele_balance_arg)
    gyper::Options::instance()->maximum_homozygous_allele_balance = args::get(maximum_homozygous_allele_balance_arg);
}


/** Epsilon 0 exponent */
using TEpsilonZeroExponent = args::ValueFlag<uint16_t>;

std::unique_ptr<TEpsilonZeroExponent>
add_arg_epsilon_0_exponent(args::ArgumentParser & parser)
{
  return std::unique_ptr<TEpsilonZeroExponent>(
    new TEpsilonZeroExponent(parser,
                        "N",
                        "Epsilon zero exponent (Default is 13)",
                        {"e", "epsilon_0_exponent"}
                        )
    );
}


void
parse_epsilon_0_exponent(TEpsilonZeroExponent & epsilon_0_exponent_arg)
{
  if (epsilon_0_exponent_arg)
    gyper::Options::instance()->epsilon_0_exponent = args::get(epsilon_0_exponent_arg);
}


/** Read chunk size argument */
using TChunkSize = args::ValueFlag<unsigned long>;

std::unique_ptr<TChunkSize>
add_arg_read_chunk_size(args::ArgumentParser & parser)
{
  return std::unique_ptr<TChunkSize>(new TChunkSize(parser, "N", "Reads are read from file in chunks. This parameter tunes how large each chunk is.", {"read_chunk_size"}));
}


void
parse_read_chunk_size(TChunkSize & read_chunk_size_arg)
{
  if (read_chunk_size_arg)
    gyper::Options::instance()->read_chunk_size = args::get(read_chunk_size_arg);
}


/** Haplotypes argument */
using THaplotypes = args::ValueFlag<std::string>;

std::unique_ptr<THaplotypes>
add_arg_haplotypes(args::ArgumentParser & parser)
{
  return std::unique_ptr<THaplotypes>(new THaplotypes(parser, "FILE", "File with a newline separated list of .hap files.", {"haplotypes"}));
}


/** Phased argument */
std::unique_ptr<args::Flag>
add_arg_phased(args::ArgumentParser & parser)
{
  return std::unique_ptr<args::Flag>(new args::Flag(parser, "PHASED", "Set to output phased haplotypes.", {"p", "phased"}));
}


void
parse_phased(args::Flag & use_phased_arg)
{
  if (use_phased_arg)
    gyper::Options::instance()->phased_output = true;
}


/** always_query_hamming_distance_one argument */
std::unique_ptr<args::Flag>
add_arg_always_query_hamming_distance_one(args::ArgumentParser & parser)
{
  return std::unique_ptr<args::Flag>(new args::Flag(parser, "PHASED", "Set to always_query_hamming_distance_one.", {"always_query_hamming_distance_one"}));
}


void
parse_always_query_hamming_distance_one(args::Flag & arg)
{
  if (arg)
    gyper::Options::instance()->always_query_hamming_distance_one = true;
}


/** Max extracted haplotypes argument */
using TMaxExtractH = args::ValueFlag<unsigned long>;

std::unique_ptr<TMaxExtractH>
add_arg_max_extracted_haplotypes(args::ArgumentParser & parser)
{
  return std::unique_ptr<TMaxExtractH>(
    new TMaxExtractH(parser,
                     "N",
                     "Maximum number of haplotypes to be extracted.",
                     {"max_extracted_haplotypes"},
                     42
                     )
    );
}


void
parse_max_extracted_haplotypes(TMaxExtractH & max_extracted_haplotypes_arg)
{
  if (max_extracted_haplotypes_arg)
    gyper::Options::instance()->max_extracted_haplotypes = args::get(max_extracted_haplotypes_arg);
}

/** Skip breaking down extracted haplotypes argument */
std::unique_ptr<args::Flag>
add_arg_skip_breaking_down_extracted_haplotypes(args::ArgumentParser & parser)
{
  return std::unique_ptr<args::Flag>(
    new args::Flag(parser, "SBDEH", "Set to skip breaking down extracted haplotypes.", {"k", "skip_breaking_down_extracted_haplotypes"})
                                    );
}


template <typename T>
bool
check_required_argument(T & arg, std::string const arg_name)
{
  if (!*arg)
  {
    std::cerr << "ERROR: Argument '" << arg_name << "' was not passed but is required." << std::endl;
    return false;
  }

  return true;
}


void
parse_command_line(args::ArgumentParser & parser, int argc, char ** argv)
{
  try
  {
    parser.ParseCLI(argc, argv);
  }
  catch (args::Help)
  {
    std::cout << parser;
    std::exit(0);
  }
  catch (args::Error e)
  {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    std::exit(1);
  }
}


} // namespace


int
main(int argc, char ** argv)
{
  // Check if there was any available command given
  if (argc <= 1)
  {
    print_default_help();
    std::exit(1);
  }

  std::string const given_command = argv[1];
  std::vector<std::string> available_commands = {"count_index_keys", "call", "check", "construct", "haplotypes", "index", "vcf_merge", "vcf_concatenate", "vcf_break_down", "vcf_update_info"};

  if (std::find(available_commands.begin(), available_commands.end(), given_command) == available_commands.end())
  {
    BOOST_LOG_TRIVIAL(error) << "[graphtyper::main] Command '" << given_command << "' is not avaiable.\n";
    print_default_help();
    std::exit(1);
  }

  // Set default logging severity
  boost::log::core::get()->set_filter
  (
    boost::log::trivial::severity >= boost::log::trivial::warning
  );

  //

  if (std::string(argv[1]) == std::string("count_index_keys"))
  {
    args::ArgumentParser count_index_keys_parser("Graphtyper's count index keys");
    auto help_arg = add_arg_help(count_index_keys_parser);
    auto command_arg = add_arg_command(count_index_keys_parser, argv[1]);
    auto graph_arg = add_arg_graph(count_index_keys_parser);
    auto index_arg = add_arg_index(count_index_keys_parser);

    parse_command_line(count_index_keys_parser, argc, argv);

    std::string index_path;

    if (*index_arg)
      index_path = args::get(*index_arg);
    else
      index_path = args::get(*graph_arg) + std::string("_gti");

    // Check if graph exists
    if (!is_file(args::get(*graph_arg)))
    {
      BOOST_LOG_TRIVIAL(error) << "[graphtyper::main] Could not find a graph located at '" << args::get(*graph_arg) << "'.";
      return 1;
    }

    // Check if index exists
    if (!is_directory(index_path))
    {
      BOOST_LOG_TRIVIAL(error) << "[graphtyper::main] Could not find an index located at '" << args::get(*graph_arg) << "'.";
      return 1;
    }

    gyper::load_graph(args::get(*graph_arg));
    gyper::load_index(index_path);

    // Get all reads
    {
      std::size_t num_keys = 0;
      std::size_t num_values = 0; // Count the number of keys

      rocksdb::Iterator* it = gyper::index.hamming0.db->NewIterator(rocksdb::ReadOptions());

      for (it->SeekToFirst(); it->Valid(); it->Next())
      {
        ++num_keys;
        num_values += gyper::value_to_labels(it->value().ToString()).size();
      }

      std::cout << "Read " << num_keys << " keys and " << num_values << " values from index." << std::endl;
      delete it;
    }
  }
  else if (std::string(argv[1]) == std::string("construct"))
  {
    args::ArgumentParser construct_parser("Graphtyper's graph construction tool");
    auto help_arg = add_arg_help(construct_parser);
    auto command_arg = add_arg_command(construct_parser, argv[1]);
    auto graph_arg = add_arg_graph(construct_parser);
    auto fasta_arg = add_arg_fasta(construct_parser);
    auto regions_arg = add_arg_region(construct_parser);
    auto vcf_arg = add_arg_vcf(construct_parser);
    auto log_arg = add_arg_log(construct_parser);

    parse_command_line(construct_parser, argc, argv);

    parse_log(*log_arg);

    bool SUCCESS = true;
    SUCCESS &= check_required_argument(command_arg, "command");
    SUCCESS &= check_required_argument(graph_arg, "graph");
    SUCCESS &= check_required_argument(fasta_arg, "fasta");
    SUCCESS &= check_required_argument(regions_arg, "regions");

    std::vector<std::string> regions;
    regions.push_back(args::get(*regions_arg));

    for (auto & reg : regions)
      reg.erase(std::remove(reg.begin(), reg.end(), ','), reg.end());

    if (!SUCCESS)
    {
      std::cerr << construct_parser;
      return 1;
    }

    // Check if FASTA exists
    if (!is_file(args::get(*fasta_arg)))
    {
      BOOST_LOG_TRIVIAL(error) << "[graphtyper::main] Could not find a FASTA file located at '" << args::get(*fasta_arg) << "'.";
      return 1;
    }

    // Check if VCF exists
    if (*vcf_arg && !is_file(args::get(*vcf_arg)))
    {
      BOOST_LOG_TRIVIAL(error) << "[graphtyper::main] Could not find a VCF file located at '" << args::get(*vcf_arg) << "'.";
      return 1;
    }

    if (*vcf_arg)
      gyper::construct_graph(args::get(*fasta_arg), args::get(*vcf_arg), regions);
    else
      gyper::construct_graph(args::get(*fasta_arg), regions);

    gyper::save_graph(args::get(*graph_arg));
  }
  else if (std::string(argv[1]) == std::string("index"))
  {
    args::ArgumentParser index_parser("Graphtyper's index construction tool.");
    auto help_arg = add_arg_help(index_parser);
    auto command_arg = add_arg_command(index_parser, argv[1]);
    auto graph_arg = add_arg_graph(index_parser);
    auto index_arg = add_arg_index(index_parser);
    auto log_arg = add_arg_log(index_parser);

    parse_command_line(index_parser, argc, argv);

    parse_log(*log_arg);

    bool SUCCESS = true;
    SUCCESS &= check_required_argument(command_arg, "command");
    SUCCESS &= check_required_argument(graph_arg, "graph");

    // Exit if it failed to get all requred arguments
    if (!SUCCESS)
    {
      std::cerr << index_parser;
      return 1;
    }

    // If the index argument was not provided, try to get it from the graph argument.
    std::string index_path;

    if (*index_arg)
      index_path = args::get(*index_arg);
    else
      index_path = args::get(*graph_arg) + std::string("_gti");

    gyper::index_graph(args::get(*graph_arg), index_path);
  }
  else if (std::string(argv[1]) == std::string("call"))
  {
    args::ArgumentParser call_parser("Graphtyper's genotype calling tool.");
    auto help_arg = add_arg_help(call_parser);
    auto command_arg = add_arg_command(call_parser, argv[1]);
    auto graph_arg = add_arg_graph(call_parser);
    auto regions_arg = add_arg_regions(call_parser);
    auto index_arg = add_arg_index(call_parser);
    auto sam_arg = add_arg_sam(call_parser);
    auto sams_arg = add_arg_sams(call_parser);
    auto segment_arg = add_arg_segment(call_parser);
    auto stats_arg = add_arg_stats(call_parser);
    auto max_index_labels_arg = add_arg_max_index_labels(call_parser);
    auto mmvd_arg = add_arg_mmvd(call_parser);
    // auto gather_unmapped_arg = add_arg_gather_unmapped(call_parser);
    auto log_arg = add_arg_log(call_parser);
    auto minimum_variant_support_arg = add_arg_min_var_sup(call_parser);
    auto minimum_variant_support_ratio_arg = add_arg_min_var_sup_ratio(call_parser);
    auto soft_cap_of_non_snps_in_100_bp_window_arg = add_arg_soft_cap_of_non_snps_in_100_bp_window(call_parser);
    auto soft_cap_of_variants_in_100_bp_window_arg = add_arg_soft_cap_of_variants_in_100_bp_window(call_parser);
    auto hard_cap_of_non_snps_in_100_bp_window_arg = add_arg_hard_cap_of_non_snps_in_100_bp_window(call_parser);
    auto hard_cap_of_variants_in_100_bp_window_arg = add_arg_hard_cap_of_variants_in_100_bp_window(call_parser);
    auto maximum_homozygous_allele_balance_arg = add_arg_maximum_homozygous_allele_balance(call_parser);
    auto epsilon_0_exponent_arg = add_arg_epsilon_0_exponent(call_parser);
    auto no_new_variants_arg = add_arg_no_new_variants(call_parser);
    auto hq_reads_arg = add_arg_hq_reads(call_parser);
    auto output_all_variants_arg = add_arg_output_all_variants(call_parser);
    auto read_chunk_size_arg = add_arg_read_chunk_size(call_parser);
    auto output_arg = add_arg_output_dir(call_parser);
    auto threads_arg = add_arg_threads(call_parser);
    auto get_sample_names_from_filename_arg = add_arg_get_sample_names_from_filename(call_parser);
    auto phased_arg = add_arg_phased(call_parser);
    auto always_query_hamming_distance_one_arg = add_arg_always_query_hamming_distance_one(call_parser);
    // auto use_read_cache_arg = add_arg_use_read_cache(call_parser);

    parse_command_line(call_parser, argc, argv);

    parse_log(*log_arg);
    parse_threads(*threads_arg);
    parse_read_chunk_size(*read_chunk_size_arg);
    // parse_use_read_cache(*use_read_cache_arg);
    parse_minimum_variant_support(*minimum_variant_support_arg);
    parse_minimum_variant_support_ratio(*minimum_variant_support_ratio_arg);
    parse_soft_cap_of_non_snps_in_100_bp_window(*soft_cap_of_non_snps_in_100_bp_window_arg);
    parse_soft_cap_of_variants_in_100_bp_window(*soft_cap_of_variants_in_100_bp_window_arg);
    parse_hard_cap_of_non_snps_in_100_bp_window(*hard_cap_of_non_snps_in_100_bp_window_arg);
    parse_hard_cap_of_variants_in_100_bp_window(*hard_cap_of_variants_in_100_bp_window_arg);
    parse_maximum_homozygous_allele_balance(*maximum_homozygous_allele_balance_arg);
    parse_epsilon_0_exponent(*epsilon_0_exponent_arg);
    parse_stats(*stats_arg);
    parse_max_index_labels(*max_index_labels_arg);
    parse_mmvd(*mmvd_arg);
    parse_get_sample_names_from_filename(*get_sample_names_from_filename_arg);
    // parse_gather_unmapped(*gather_unmapped_arg);
    parse_no_new_variants(*no_new_variants_arg);
    parse_hq_reads(*hq_reads_arg);
    parse_output_all_variants(*output_all_variants_arg);
    parse_phased(*phased_arg);
    parse_always_query_hamming_distance_one(*always_query_hamming_distance_one_arg);

    bool SUCCESS = true;
    SUCCESS &= check_required_argument(command_arg, "command");
    SUCCESS &= check_required_argument(graph_arg, "graph");
    SUCCESS &= check_required_argument(regions_arg, "regions");

    std::vector<std::string> regions = args::get(*regions_arg);

    for (auto & reg : regions)
      reg.erase(std::remove(reg.begin(), reg.end(), ','), reg.end());

    // Either sam or sams is required
    if (!*sam_arg && !*sams_arg)
    {
      BOOST_LOG_TRIVIAL(error) << "[graphtyper::main] Either the 'sam' or 'sams' argument must be passed.";
      SUCCESS = false;
    }

    if (!SUCCESS)
    {
      std::cerr << call_parser;
      return 1; // Exit if it failed to get all requred arguments
    }

    if (!is_directory(args::get(*output_arg)))
      mkdir(args::get(*output_arg).c_str(), 0755);

    // If the index_path argument was not provided, try to get it from the graph argument.
    std::string index_path;

    if (*index_arg)
      index_path = args::get(*index_arg);
    else
      index_path = args::get(*graph_arg) + std::string("_gti");

    std::vector<std::string> fasta_segments;

    if (*segment_arg)
    {
      fasta_segments.push_back(args::get(*segment_arg));
      gyper::Options::instance()->is_segment_calling = true;
    }

    // Parse SAM files
    std::vector<std::string> sams;

    if (*sam_arg)
      sams.push_back(args::get(*sam_arg));

    if (*sams_arg)
    {
      std::ifstream file_in(args::get(*sams_arg));
      std::string line;

      while (std::getline(file_in, line))
        sams.push_back(line);
    }

    // Check if graph exists
    if (!is_file(args::get(*graph_arg)))
    {
      BOOST_LOG_TRIVIAL(error) << "[graphtyper::main] Could not find a graph located at '" << args::get(*graph_arg) << "'.";
      return 1;
    }

    // Check if index exists
    if (!is_directory(index_path))
    {
      BOOST_LOG_TRIVIAL(error) << "[graphtyper::main] Could not find an index located at '" << args::get(*graph_arg) << "'.";
      return 1;
    }

    gyper::call(sams,
                args::get(*graph_arg),
                index_path,
                regions,
                std::move(fasta_segments),
                args::get(*output_arg)
                );
  }
  else if (std::string(argv[1]) == std::string("haplotypes"))
  {
    args::ArgumentParser haplotypes_parser("Graphtyper's haplotype extraction tool.");
    auto help_arg = add_arg_help(haplotypes_parser);
    auto command_arg = add_arg_command(haplotypes_parser, argv[1]);
    auto graph_arg = add_arg_graph(haplotypes_parser);
    auto haplotypes_arg = add_arg_haplotypes(haplotypes_parser);
    auto log_arg = add_arg_log(haplotypes_parser);
    auto max_extracted_haplotypes_arg = add_arg_max_extracted_haplotypes(haplotypes_parser);
    auto skip_breaking_down_extracted_haplotypes = add_arg_skip_breaking_down_extracted_haplotypes(haplotypes_parser);
    auto output_arg = add_arg_output(haplotypes_parser);
    auto region_arg = add_arg_region_val(haplotypes_parser);
    auto output_all_variants_arg = add_arg_output_all_variants(haplotypes_parser, "haplotypes");
    auto minimum_variant_support_arg = add_arg_min_var_sup(haplotypes_parser);

    parse_command_line(haplotypes_parser, argc, argv);

    parse_log(*log_arg);
    parse_max_extracted_haplotypes(*max_extracted_haplotypes_arg);
    parse_output_all_variants(*output_all_variants_arg);
    parse_minimum_variant_support(*minimum_variant_support_arg);

    bool SUCCESS = true;
    SUCCESS &= check_required_argument(command_arg, "command");
    SUCCESS &= check_required_argument(graph_arg, "graph");
    SUCCESS &= check_required_argument(haplotypes_arg, "haplotypes");

    if (!SUCCESS)
    {
      std::cerr << haplotypes_parser;
      return 1; // Exit if it failed to get all requred arguments
    }

    if (not is_file(args::get(*haplotypes_arg)))
    {
      BOOST_LOG_TRIVIAL(error) << "[graphtyper::main] Haplotype file '" << args::get(*haplotypes_arg) << "' was not found.";
      return 1;
    }

    std::string output_file("haps.vcf.gz");

    if (*skip_breaking_down_extracted_haplotypes)
      gyper::Options::instance()->skip_breaking_down_extracted_haplotypes = true;

    if (*output_arg)
      output_file = args::get(*output_arg);

    // Check if graph exists
    if (!is_file(args::get(*graph_arg)))
    {
      BOOST_LOG_TRIVIAL(error) << "[graphtyper::main] Could not find a graph located at '" << args::get(*graph_arg) << "'.";
      return 1;
    }

    gyper::extract_to_vcf(args::get(*graph_arg),
                          args::get(*haplotypes_arg),
                          output_file,
                          args::get(*region_arg)
                          );
  }
  else if (std::string(argv[1]) == std::string("vcf_merge"))
  {
    args::ArgumentParser vcf_merge_parser("Graphtyper's VCF merging tool. This tool is used to merge Graphtyper's genotype calls of multiple sample pools.");
    auto help_arg = add_arg_help(vcf_merge_parser);
    auto command_arg = add_arg_command(vcf_merge_parser, argv[1]);
    auto log_arg = add_arg_log(vcf_merge_parser);
    auto output_arg = add_arg_output(vcf_merge_parser);
    auto vcfs_arg = add_arg_vcfs(vcf_merge_parser);

    parse_command_line(vcf_merge_parser, argc, argv);

    parse_log(*log_arg);

    bool SUCCESS = true;
    SUCCESS &= check_required_argument(command_arg, "command");
    SUCCESS &= check_required_argument(vcfs_arg, "vcfs");

    if (!SUCCESS)
    {
      std::cerr << vcf_merge_parser;
      return 1; // Exit if it failed to get all requred arguments
    }

    std::string output = "-";

    if (*output_arg)
      output = args::get(*output_arg);

    gyper::vcf_merge(args::get(*vcfs_arg), output);
  }
  else if (std::string(argv[1]) == std::string("vcf_concatenate"))
  {
    args::ArgumentParser vcf_concatenate_parser("Graphtyper's VCF concatenation tool. This tool is used to concatenate one or more Graphtyper VCF files. This tool is useful when appending haplotype VCF to the newly discovered variants VCF.");
    auto help_arg = add_arg_help(vcf_concatenate_parser);
    auto command_arg = add_arg_command(vcf_concatenate_parser, argv[1]);
    auto log_arg = add_arg_log(vcf_concatenate_parser);
    auto output_arg = add_arg_output(vcf_concatenate_parser);
    auto vcfs_arg = add_arg_vcfs(vcf_concatenate_parser);
    auto no_sort_arg = add_arg_no_sort(vcf_concatenate_parser);
    auto region_arg = add_arg_region_val(vcf_concatenate_parser);

    parse_command_line(vcf_concatenate_parser, argc, argv);

    parse_log(*log_arg);

    bool SUCCESS = true;
    SUCCESS &= check_required_argument(command_arg, "command");
    SUCCESS &= check_required_argument(vcfs_arg, "vcfs");

    if (!SUCCESS)
    {
      std::cerr << vcf_concatenate_parser;
      return 1; // Exit if it failed to get all requred arguments
    }

    std::string output = "-";

    if (*output_arg)
      output = args::get(*output_arg);

    std::string region = ".";

    if (*region_arg)
      region = args::get(*region_arg);

    gyper::vcf_concatenate(args::get(*vcfs_arg), output, *no_sort_arg, region);
  }
  else if (std::string(argv[1]) == std::string("vcf_break_down"))
  {
    args::ArgumentParser vcf_break_down_parser("Graphtyper's VCF breaking down tool. Useful if you have a VCF file with a lot of multiallelic sites and would like to find all SNPs and indels which are part of more complex variants.");
    auto help_arg = add_arg_help(vcf_break_down_parser);
    auto command_arg = add_arg_command(vcf_break_down_parser, argv[1]);
    auto graph_arg = add_arg_graph(vcf_break_down_parser);
    auto log_arg = add_arg_log(vcf_break_down_parser);
    auto output_arg = add_arg_output(vcf_break_down_parser);
    auto vcf_pos_arg = add_arg_vcf_pos(vcf_break_down_parser);
    auto region_arg = add_arg_region_val(vcf_break_down_parser);

    parse_command_line(vcf_break_down_parser, argc, argv);

    parse_log(*log_arg);

    bool SUCCESS = true;
    SUCCESS &= check_required_argument(command_arg, "command");
    SUCCESS &= check_required_argument(graph_arg, "graph");
    SUCCESS &= check_required_argument(vcf_pos_arg, "vcf");

    if (!SUCCESS)
    {
      std::cerr << vcf_break_down_parser;
      return 1; // Exit if it failed to get all requred arguments
    }

    std::string output = "-";

    if (*output_arg)
      output = args::get(*output_arg);

    // Check if graph exists
    if (!is_file(args::get(*graph_arg)))
    {
      BOOST_LOG_TRIVIAL(error) << "[graphtyper::main] Could not find a graph located at '" << args::get(*graph_arg) << "'.";
      return 1;
    }

    std::string region = ".";

    if (*region_arg)
      region = args::get(*region_arg);

    // load the graph
    gyper::load_graph(args::get(*graph_arg));
    gyper::vcf_break_down(args::get(*vcf_pos_arg), output, region);
  }
  else if (std::string(argv[1]) == std::string("vcf_update_info"))
  {
    args::ArgumentParser vcf_update_info_parser("Graphtyper's VCF update info tool. Useful if you have a Graphtyper VCF file with errors in the INFO fields you want fixed.");
    auto help_arg = add_arg_help(vcf_update_info_parser);
    auto command_arg = add_arg_command(vcf_update_info_parser, argv[1]);
    auto log_arg = add_arg_log(vcf_update_info_parser);
    auto output_arg = add_arg_output(vcf_update_info_parser);
    auto vcf_pos_arg = add_arg_vcf_pos(vcf_update_info_parser);

    parse_command_line(vcf_update_info_parser, argc, argv);

    parse_log(*log_arg);

    bool SUCCESS = true;
    SUCCESS &= check_required_argument(command_arg, "command");
    SUCCESS &= check_required_argument(vcf_pos_arg, "vcf");

    if (!SUCCESS)
    {
      std::cerr << vcf_update_info_parser;
      return 1; // Exit if it failed to get all requred arguments
    }

    std::string output = "-";

    if (*output_arg)
      output = args::get(*output_arg);

    gyper::vcf_update_info(args::get(*vcf_pos_arg), output);
  }

  if (gyper::Options::instance()->sink)
    gyper::Options::instance()->sink->flush(); // Flush sink if there is one

  return 0;
}
