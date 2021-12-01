#include <fstream>
#include <numeric> // std::iota
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include <graphtyper/graph/absolute_position.hpp>
#include <graphtyper/graph/constructor.hpp>
#include <graphtyper/graph/genomic_region.hpp>
#include <graphtyper/graph/graph_serialization.hpp>
#include <graphtyper/graph/haplotype_extractor.hpp>
#include <graphtyper/index/indexer.hpp>
#include <graphtyper/index/ph_index.hpp>
#include <graphtyper/typer/caller.hpp>
#include <graphtyper/typer/variant_map.hpp>
#include <graphtyper/typer/vcf.hpp>
#include <graphtyper/typer/vcf_operations.hpp>
#include <graphtyper/utilities/filesystem.hpp>
#include <graphtyper/utilities/genotype.hpp>
#include <graphtyper/utilities/graph_help_functions.hpp>
#include <graphtyper/utilities/hts_parallel_reader.hpp>
#include <graphtyper/utilities/logging.hpp>
#include <graphtyper/utilities/options.hpp>
#include <graphtyper/utilities/string.hpp>
#include <graphtyper/utilities/system.hpp>

namespace
{
struct AlleleNode
{
  int ac{0};
  int digit{0};

  bool operator<(AlleleNode const & o) const
  {
    return ac < o.ac;
  }

  bool operator==(AlleleNode const & o) const
  {
    return ac == o.ac;
  }
};

} // namespace

namespace gyper
{
void genotype_hla(std::string ref_path,
                  std::string const & hla_vcf_fn,
                  std::string const & interval_fn,
                  std::vector<std::string> const & sams,
                  std::vector<std::string> const & sam_index_paths,
                  std::vector<double> const & avg_cov_by_readlen,
                  GenomicRegion const & genomic_region,
                  std::string const & output_path,
                  bool const is_copy_reference)
{
  // TODO: If the reference is only Ns then output an empty vcf with the sample names
  // TODO: Extract the reference sequence and use that to discover directly from BAM
  gyper::Options const & copts = *(Options::const_instance());
  bool constexpr is_writing_calls_vcf{true};
  bool constexpr is_writing_hap{false};

  long const NUM_SAMPLES = sams.size();

  std::string const region = genomic_region.to_string();
  print_log(log_severity::info, "HLA genotyping region ", region);
  print_log(log_severity::info, "Path to genome is '", ref_path, "'");
  print_log(log_severity::info, "Running with up to ", copts.threads, " threads.");
  print_log(log_severity::info, "Copying data from ", NUM_SAMPLES, " input SAM/BAM/CRAMs to local disk.");
  std::string tmp = create_temp_dir(genomic_region);

  print_log(log_severity::info, "Temporary folder is ", tmp);

  // Create directories
  mkdir(output_path.c_str(), 0755);
  mkdir((output_path + "/" + genomic_region.chr).c_str(), 0755);

  // Copy reference genome to temporary directory
  if (is_copy_reference)
  {
    print_log(log_severity::info, "Copying reference genome FASTA and its index to temporary folder.");

    filesystem::copy_file(ref_path, tmp + "/genome.fa");
    filesystem::copy_file(ref_path + ".fai", tmp + "/genome.fa.fai");

    ref_path = tmp + "/genome.fa";
  }

  std::vector<std::string> shrinked_sams;

  if (copts.no_bamshrink)
  {
    shrinked_sams = sams;
  }
  else
  {
    print_log(log_severity::info, "Running BamShrink.");
    std::string bamshrink_ref_path;

    if (copts.force_use_input_ref_for_cram_reading)
      bamshrink_ref_path = ref_path;

    if (interval_fn.size() > 0)
      shrinked_sams = run_bamshrink_multi(sams, ref_path, interval_fn, avg_cov_by_readlen, tmp);
    else
      shrinked_sams = run_bamshrink(sams, sam_index_paths, bamshrink_ref_path, genomic_region, avg_cov_by_readlen, tmp);

    if (!copts.no_sample_name_reordering)
      std::sort(shrinked_sams.begin(), shrinked_sams.end()); // Sort by input filename

    run_samtools_merge(shrinked_sams, tmp);
  }

  GenomicRegion padded_region(genomic_region);
  GenomicRegion unpadded_region(padded_region);
  padded_region.pad(1000l);

  // Iteration 1 out of 1
  {
    print_log(log_severity::info, "Genotype calling step starting.");
    std::string const output_vcf = tmp + "/it1/final.vcf.gz";
    std::string const out_dir = tmp + "/it1";
    mkdir(out_dir.c_str(), 0755);
    print_log(log_severity::info, "Padded region is: ", padded_region.to_string());

    {
      bool constexpr is_sv_graph{false};
      bool constexpr use_index{true};

      print_log(log_severity::info, "Constructing graph.");
      gyper::construct_graph(ref_path, hla_vcf_fn, padded_region.to_string(), is_sv_graph, use_index);
      print_log(log_severity::info, "Calculating contig offsets.");
      absolute_pos.calculate_offsets(gyper::graph.contigs);
    }

    print_log(log_severity::info, "Reading input HLA VCF.");
    gyper::Vcf hla(READ_BGZF_MODE, hla_vcf_fn);
    hla.read();

    hla.set_filemode(WRITE_BGZF_MODE);
    hla.filename = tmp + "/graphtyper.hla.vcf.gz";
    hla.write(); // Just for debbuging, remove later

    std::unordered_map<long, std::pair<uint32_t, uint32_t>> event2hap_gt;

    {
      uint32_t v{0};
      uint32_t h{0};
      assert(graph.ref_nodes.size() > 0);

      for (uint32_t r{0}; r < (graph.ref_nodes.size() - 1); ++r)
      {
        assert(v < graph.var_nodes.size());
        auto const & ref_node = graph.ref_nodes[r];
        uint32_t const out_degree = ref_node.out_degree();

        for (uint32_t v_e{0}; v_e < out_degree; ++v_e)
        {
          assert(v + v_e < graph.var_nodes.size());
          auto const & var_node = graph.var_nodes[v + v_e];

          for (long const & event : var_node.events)
          {
            if (event > 0)
              event2hap_gt[event] = std::pair<uint32_t, uint32_t>(h, v_e);
          }
        }

        ++h;
        v += ref_node.out_degree();
      }
    }

    std::unordered_set<uint32_t> exon_haps;

    for (Variant const & var : hla.variants)
    {
      if (var.infos.count("FEATURE") == 0 || var.infos.count("GT_ID") == 0)
      {
        // Records without INFO/FEATURE or/and INFO/GT_ID should be ignored for calling (but are included in the graph)
        continue;
      }

      if (var.infos.at("FEATURE") != "exon")
      {
        continue;
      }

      long const gt_id = std::stol(var.infos.at("GT_ID"));
      auto find_it = event2hap_gt.find(gt_id);
      assert(find_it != event2hap_gt.end());
      exon_haps.insert(find_it->second.first);
    }

    print_log(log_severity::info, "Got ", exon_haps.size(), " exonic variant records.");
    std::vector<std::unordered_map<uint32_t, uint32_t>> allele_hap_gts(hla.sample_names.size());

    for (long s{0}; s < static_cast<long>(hla.sample_names.size()); ++s)
    {
      std::unordered_map<uint32_t, uint32_t> & allele_hap_gt = allele_hap_gts[s];

      for (Variant const & var : hla.variants)
      {
        auto find_feature_it = var.infos.find("FEATURE");
        auto find_gt_id_it = var.infos.find("GT_ID");

        if (find_feature_it == var.infos.end() || find_gt_id_it == var.infos.end())
        {
          // Records without INFO/FEATURE or/and INFO/GT_ID should be ignored for calling (but are included in the
          // graph)
          continue;
        }

        if (find_feature_it->second != "exon")
        {
          print_log(log_severity::debug, __HERE__, " Skipping ", var.infos.at("FEATURE"));
          continue;
        }

        long const gt_id = std::stol(find_gt_id_it->second);
        assert(s < static_cast<long>(var.calls.size()));
        SampleCall const & sample_call = var.calls[s];
        assert(sample_call.coverage.size() == 2);

        auto find_it = event2hap_gt.find(gt_id);
        assert(find_it != event2hap_gt.end());

        if (sample_call.coverage[0] == 0)
        {
          allele_hap_gt.insert(find_it->second);
        }
      }

      // Add reference genotypes
      for (auto const & exon_hap : exon_haps)
      {
        if (allele_hap_gt.count(exon_hap) == 0)
          allele_hap_gt.insert({exon_hap, 0});
      }
    }

    print_log(log_severity::info, "Read ", hla.sample_names.size(), " alleles.");

#ifndef NDEBUG
    // Save graph in debug mode
    save_graph(out_dir + "/graph");
#endif // NDEBUG

    PHIndex ph_index = index_graph(gyper::graph);

    print_log(log_severity::info, "Finished indexing graph.");

    std::string reference_fn{}; // empty by default
    std::map<std::pair<uint16_t, uint16_t>, std::map<std::pair<uint16_t, uint16_t>, int8_t>> ph;

    if (Options::const_instance()->force_use_input_ref_for_cram_reading)
      reference_fn = ref_path;

    std::vector<std::string> paths = gyper::call(shrinked_sams,
                                                 avg_cov_by_readlen,
                                                 "", // graph_path
                                                 ph_index,
                                                 out_dir,
                                                 "",  // reference
                                                 ".", // region
                                                 nullptr,
                                                 ph,
                                                 is_writing_calls_vcf,
                                                 is_writing_hap,
                                                 &allele_hap_gts);

    print_log(log_severity::info, "Merging output VCFs.");

    if (copts.force_ignore_segment)
    {
      vcf_merge_and_break(paths, tmp + "/graphtyper.vcf.gz", genomic_region.to_string(), false, false, true);
    }
    else
    {
      Vcf hla_vcf;
      vcf_merge_and_return(hla_vcf, paths, tmp + "/graphtyper.vcf.gz");
      assert(hla_vcf.variants.size() == 1);
      Variant & var = hla_vcf.variants[0];
      int constexpr MAX_ALLELES{80}; // Maximum number of allowed HLA alleles

      {
        // Fix sequences names
        auto is_pass_alt = var.generate_infos();
        assert(var.seqs.size() == hla.sample_names.size());
        assert((is_pass_alt.size() + 1u) == hla.sample_names.size());
        std::vector<std::vector<char>> new_seqs;

        for (long s{0}; s < static_cast<long>(hla.sample_names.size()); ++s)
        {
          if (s > 0 && is_pass_alt[s - 1] == 0)
            continue;

          std::string const & hla_allele = hla.sample_names[s];
          std::vector<char> seq;
          seq.push_back('<');
          std::copy(hla_allele.begin(), hla_allele.end(), std::back_inserter(seq));
          seq.push_back('>');
          new_seqs.push_back(seq);
        }

        if (new_seqs.size() == 1 && hla.sample_names.size() >= 2)
        {
          // special case where only reference allele is called, in this case also add some other sequence
          is_pass_alt[0] = true;
          std::string const & hla_allele = hla.sample_names[1];
          std::vector<char> seq;
          seq.push_back('<');
          std::copy(hla_allele.begin(), hla_allele.end(), std::back_inserter(seq));
          seq.push_back('>');
          new_seqs.push_back(seq);
        }

        assert(new_seqs.size() >= 2);

        for (auto const & new_seq : new_seqs)
        {
          print_log(log_severity::info, "Called contig sequence: ", std::string(new_seq.begin(), new_seq.end()));
        }

        long const old_cnum = var.seqs.size();
        long const cnum = new_seqs.size();

        for (auto & call : var.calls)
        {
          call.coverage.resize(cnum);

          std::vector<uint8_t> new_phred;
          new_phred.reserve((cnum * (cnum + 1)) / 2);

          for (long y{0}; y < old_cnum; ++y)
          {
            if (y > 0 && is_pass_alt[y - 1] == 0)
              continue;

            for (long x{0}; x <= y; ++x)
            {
              if (x > 0 && is_pass_alt[x - 1] == 0)
                continue;

              long const i = to_index(x, y);
              new_phred.push_back(call.phred[i]);
            }
          }

          assert(static_cast<long>(new_phred.size()) == ((cnum * (cnum + 1)) / 2));

          auto min_it = std::min_element(new_phred.begin(), new_phred.end());

          if (min_it != new_phred.end() && *min_it > 0)
          {
            for (auto & pl : new_phred)
              pl -= *min_it;
          }

          call.phred = new_phred;
        }

        var.seqs = std::move(new_seqs);
        var.stats = VarStats();
        var.infos.clear();
        var.generate_infos();
      }

      hla_vcf.open_for_writing(Options::const_instance()->threads);
      hla_vcf.write_header();
      assert(hla_vcf.variants.size() > 0);

      if (var.seqs.size() <= MAX_ALLELES)
      {
        hla_vcf.write_record(hla_vcf.variants[0], ".all", false, false);
      }
      else
      {
        print_log(log_severity::info,
                  "Skipping all HLA allele calling because there are more than ",
                  MAX_ALLELES,
                  " called HLA alleles (",
                  hla_vcf.variants[0].seqs.size(),
                  ")");
      }

      // Create a tree of HLA alleles
      std::unordered_set<std::string> common_4digit;
      int num_2digit_seqs{1};
      bool is_retry_4digit{false}; // Is set to try if we are retrying 4-digit calling

      // Add 2 and 4 digit variant
      for (int d{2}; d < 6; d += 2)
      {
        // Maps alleles to their index position in new_var.seqs vector
        std::unordered_map<std::string, long> seen_alleles;
        std::vector<long> old2new(var.seqs.size(), 0);
        Variant new_var;
        new_var.abs_pos = var.abs_pos;

        for (long a{0}; a < static_cast<long>(var.seqs.size()); ++a)
        {
          auto const & seq = var.seqs[a]; // i.e. <HLA-DRB1*15:01:01:01>
          std::optional<std::string> new_allele;

          if (d == 4 && is_retry_4digit)
          {
            // special case: We are retrying 4-digit HLA calling.
            // In this case we must first check if 4-digit allele is common
            assert(seq.begin() != seq.end());

            auto find_it = find_nth_occurence(seq.begin(), seq.end(), ':', 2);
            std::string four_digit_allele(seq.begin(), find_it);
            assert(!four_digit_allele.empty());

            if (four_digit_allele.empty() || four_digit_allele.back() != '>')
              four_digit_allele.push_back('>');

            if (common_4digit.count(four_digit_allele) > 0)
            {
              // 4digit HLA allele is common, just add 4-digit version
#ifndef NDEBUG
              print_log(log_severity::debug, __HERE__, " Adding 4digit ", four_digit_allele);
#endif // NDEBUG

              new_allele = std::move(four_digit_allele);
            }
            else
            {
              // 4digit HLA allele is not common, add XX version
              auto find_2digit_it = find_nth_occurence(seq.begin(), seq.end(), ':', 1);
              std::string two_digit_allele(seq.begin(), find_2digit_it);
              assert(!two_digit_allele.empty());
              assert(two_digit_allele.back() != '>');
              two_digit_allele += std::string(":XX>");
              // print_log(log_severity::info, __HERE__, " Adding 2digit ", two_digit_allele);
              new_allele = std::move(two_digit_allele);
            }
          }
          else
          {
            auto find_it = find_nth_occurence(seq.begin(), seq.end(), ':', d / 2);
            new_allele = std::string(seq.begin(), find_it);

            if (new_allele->empty() || new_allele->back() != '>')
              new_allele->push_back('>');
          }

          assert(new_allele);
          auto find_allele_it = seen_alleles.find(new_allele.value());

          if (find_allele_it == seen_alleles.end())
          {
            seen_alleles[new_allele.value()] = new_var.seqs.size();
            old2new[a] = new_var.seqs.size();
            new_var.seqs.push_back(std::vector<char>(new_allele->begin(), new_allele->end()));
          }
          else
          {
            old2new[a] = find_allele_it->second;
          }
        }

        assert(seen_alleles.size() == new_var.seqs.size());

        if (new_var.seqs.size() <= 1)
        {
          print_log(log_severity::info, "Skipping ", d, "-digit calling because there is only a single allele called.");
          continue;
        }

        new_var.calls.reserve(var.calls.size());

        // bin PHRED values in calls
        for (long s{0}; s < static_cast<long>(var.calls.size()); ++s)
        {
          SampleCall new_call = bin_phred(new_var, var, var.calls[s], old2new);
          new_var.calls.push_back(std::move(new_call));
        }

        new_var.generate_infos();
        bool const is_skipping = static_cast<long>(new_var.seqs.size()) > MAX_ALLELES;

        if (!is_skipping)
        {
          print_log(log_severity::info, d, "-digit calling with ", new_var.seqs.size(), " alleles.");
          hla_vcf.write_record(new_var, "." + std::to_string(d) + "digit", false, false);
        }
        else if (d == 2)
        {
          // If we are skipping 2-digit then what are we even doing
          print_log(log_severity::warning,
                    "In ",
                    d,
                    "-digit calling there are more than ",
                    MAX_ALLELES,
                    " alleles (",
                    new_var.seqs.size(),
                    ") but I will try to write the record anyway.");

          hla_vcf.write_record(new_var, "." + std::to_string(d) + "digit", false, false);
        }
        else
        {
          print_log(log_severity::info,
                    "Skipping ",
                    d,
                    "-digit calling because there are more than ",
                    MAX_ALLELES,
                    " alleles (",
                    new_var.seqs.size(),
                    ")");
        }

        if (d == 2)
        {
          // take note how many 2digit alleles are called altogether
          num_2digit_seqs = new_var.seqs.size();
        }
        else if (d == 4 && is_skipping && !is_retry_4digit && MAX_ALLELES > num_2digit_seqs)
        {
          // Find top MAX_ALLELES-num_2digit_seqs 4 digit alleles
          auto const & per_allele = new_var.stats.per_allele;
          int const n_alleles = static_cast<int>(per_allele.size());
          assert(n_alleles == static_cast<int>(new_var.seqs.size()));

          // Get allele count of each HLA allele
          std::vector<uint32_t> ac(per_allele.size());

          for (int i{0}; i < n_alleles; ++i)
            ac[i] = per_allele[i].pass_ac;

          assert(per_allele.size() == new_var.seqs.size());

          // Get a vector with indices for every allele
          std::vector<uint32_t> idx(new_var.seqs.size());
          std::iota(idx.begin(), idx.end(), 0);

          // Sort the indices by their allele count (highest first)
          std::sort(idx.begin(), idx.end(), [&ac](uint32_t idx1, uint32_t idx2) { return ac[idx1] > ac[idx2]; });

          for (int j{0}; j < (MAX_ALLELES - num_2digit_seqs); ++j)
          {
            if (j >= static_cast<int>(idx.size())) // edge case
              break;

            uint32_t const index = idx[j];
            assert(index < new_var.seqs.size());
            assert(index < ac.size());

            // If there aren't any reliable 4digit calls then don't add that allele
            if (ac[index] == 0)
              continue;

#ifndef NDEBUG
            print_log(
              log_severity::info,
              __HERE__,
              " Common allele: ",
              std::string_view(reinterpret_cast<char *>(new_var.seqs[index].data()), new_var.seqs[index].size()),
              " with ac = ",
              ac[index]);
#endif // NDEBUG

            // Add the most common HLA alleles to 'common_4digit' set
            common_4digit.insert(std::string(new_var.seqs[index].begin(), new_var.seqs[index].end()));
          }

          d -= 2;                 // Retry 4-digit calling
          is_retry_4digit = true; // Prevents endless loop of retries
        }
      }

      hla_vcf.close_vcf_file();
      hla_vcf.write_tbi_index();
    }
  }

  // Copy final VCFs
  auto copy_vcf_to_system = [&](std::string const & name, std::string const & extension) -> void
  {
    filesystem::path src = tmp + "/" + name + extension;
    filesystem::path dest = output_path + "/" + genomic_region.to_file_string() +
                            ((name == "graphtyper.hla.vcf.gz") ? ".hla" : "") + ".vcf.gz" + extension;

    filesystem::copy_file(src, dest, filesystem::copy_options::overwrite_existing);
  };

  std::string const index_ext = copts.is_csi ? ".csi" : ".tbi";

  copy_vcf_to_system("graphtyper.vcf.gz", "");        // Copy final VCF
  copy_vcf_to_system("graphtyper.vcf.gz", index_ext); // Copy tabix index for final VCF
  // copy_vcf_to_system("graphtyper.hla.vcf.gz", ""); // Copy final HLA VCF

  if (!copts.no_cleanup)
  {
    print_log(log_severity::info, "Cleaning up temporary files.");
    remove_file_tree(tmp);
  }
  else
  {
    print_log(log_severity::info, "Temporary files left: ", tmp);
  }

  {
    std::ostringstream ss;

    ss << output_path << "/" << genomic_region.chr << "/" << std::setw(9) << std::setfill('0')
       << (genomic_region.begin + 1) << '-' << std::setw(9) << std::setfill('0') << genomic_region.end << ".vcf.gz";

    print_log(log_severity::info, "Finished! Output written at: ", ss.str());
  }

  // free memory
  graph = Graph();
}

void genotype_hla_regions(std::string const & ref_path,
                          std::string const & hla_vcf,
                          std::string const & interval_fn,
                          std::vector<std::string> const & sams,
                          std::vector<std::string> const & sam_index_paths,
                          std::vector<double> const & avg_cov_by_readlen,
                          std::vector<gyper::GenomicRegion> const & regions,
                          std::string const & output_path,
                          bool const is_copy_reference)
{
  for (auto const & region : regions)
  {
    genotype_hla(ref_path,
                 hla_vcf,
                 interval_fn,
                 sams,
                 sam_index_paths,
                 avg_cov_by_readlen,
                 region,
                 output_path,
                 is_copy_reference);
  }
}

} // namespace gyper
