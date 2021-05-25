#include <fstream>
#include <string>
#include <sstream>
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
#include <graphtyper/utilities/options.hpp>
#include <graphtyper/utilities/genotype.hpp>
#include <graphtyper/utilities/graph_help_functions.hpp>
#include <graphtyper/utilities/hts_parallel_reader.hpp>
#include <graphtyper/utilities/system.hpp>

#include <boost/log/trivial.hpp>


namespace gyper
{

void
genotype_hla(std::string ref_path,
             std::string const & hla_vcf,
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
  BOOST_LOG_TRIVIAL(info) << "HLA genotyping region " << region;
  BOOST_LOG_TRIVIAL(info) << "Path to genome is '" << ref_path << "'";
  BOOST_LOG_TRIVIAL(info) << "Running with up to " << copts.threads << " threads.";
  BOOST_LOG_TRIVIAL(info) << "Copying data from " << NUM_SAMPLES << " input SAM/BAM/CRAMs to local disk.";
  std::string tmp = create_temp_dir(genomic_region);

  BOOST_LOG_TRIVIAL(info) << "Temporary folder is " << tmp;

  // Create directories
  mkdir(output_path.c_str(), 0755);
  mkdir((output_path + "/" + genomic_region.chr).c_str(), 0755);

  // Copy reference genome to temporary directory
  if (is_copy_reference)
  {
    BOOST_LOG_TRIVIAL(info) << "Copying reference genome FASTA and its index to temporary folder.";

    filesystem::copy_file(ref_path, tmp + "/genome.fa");
    filesystem::copy_file(ref_path + ".fai", tmp + "/genome.fa.fai");

    ref_path = tmp + "/genome.fa";
  }

  std::vector<std::string> shrinked_sams;

  if (copts.no_bamshrink)
  {
    shrinked_sams = std::move(sams);
  }
  else
  {
    BOOST_LOG_TRIVIAL(info) << __HERE__ << " Running BamShrink.";
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
    BOOST_LOG_TRIVIAL(info) << "Genotype calling step starting.";
    std::string const output_vcf = tmp + "/it1/final.vcf.gz";
    std::string const out_dir = tmp + "/it1";
    mkdir(out_dir.c_str(), 0755);
    BOOST_LOG_TRIVIAL(info) << "Padded region is: " << padded_region.to_string();

    {
      bool constexpr is_sv_graph{false};
      bool constexpr use_index{true};

      BOOST_LOG_TRIVIAL(info) << "Constructing graph.";

      gyper::construct_graph(ref_path,
                             hla_vcf,
                             padded_region.to_string(),
                             is_sv_graph,
                             use_index);

      BOOST_LOG_TRIVIAL(info) << "Calculating contig offsets.";

      absolute_pos.calculate_offsets(gyper::graph.contigs);
    }

    BOOST_LOG_TRIVIAL(info) << "Reading input HLA VCF.";
    gyper::Vcf hla(READ_BGZF_MODE, hla_vcf);
    hla.read();

    hla.set_filemode(WRITE_BGZF_MODE);
    hla.filename = tmp + "/graphtyper.hla.vcf.gz";
    hla.write(); // Just for debbuging, remove later

    std::unordered_map<long, std::pair<uint32_t, uint32_t> > event2hap_gt;

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
        //Haplotype hap;
        //hap.add_genotype(Genotype(var_node.get_label().order, ref_node.out_degree(), uint32_t(v)));
        //haplotypes.push_back(std::move(hap));

        ++h;
        v += ref_node.out_degree();
      }
    }

    std::unordered_set<uint32_t> exon_haps;

    for (Variant const & var : hla.variants)
    {
      if (var.infos.count("FEATURE") == 0 || var.infos.at("FEATURE") != "exon")
        continue;

      long const gt_id = std::stol(var.infos.at("GT_ID"));
      auto find_it = event2hap_gt.find(gt_id);
      assert(find_it != event2hap_gt.end());

      exon_haps.insert(find_it->second.first);
    }

    std::vector<std::unordered_map<uint32_t, uint32_t> > allele_hap_gts(hla.sample_names.size());
    //std::vector<std::unordered_set<long> > allele_events(hla.sample_names.size());

    for (long s{0}; s < static_cast<long>(hla.sample_names.size()); ++s)
    {
      std::unordered_map<uint32_t, uint32_t> & allele_hap_gt = allele_hap_gts[s];

      for (Variant const & var : hla.variants)
      {
        if (var.infos.count("FEATURE") == 0 || var.infos.at("FEATURE") != "exon")
        {
          //BOOST_LOG_TRIVIAL(info) << __HERE__ << " Skipping " << var.infos.at("FEATURE");
          continue;
        }

        long const gt_id = std::stol(var.infos.at("GT_ID"));
        assert(s < static_cast<long>(var.calls.size()));
        SampleCall const & sample_call = var.calls[s];
        assert(sample_call.coverage.size() == 2);

        auto find_it = event2hap_gt.find(gt_id);
        assert(find_it != event2hap_gt.end());

        if (sample_call.coverage[0] == 0)
        {
          allele_hap_gt.insert(find_it->second);
          //allele_events[s].insert(gt_id);
        }
        //else
        //{
        //  allele_events[s].insert(-gt_id);
        //}
      }

      // Add reference genotypes
      for (auto const & exon_hap : exon_haps)
      {
        if (allele_hap_gt.count(exon_hap) == 0)
          allele_hap_gt.insert({exon_hap, 0});
      }
    }


    BOOST_LOG_TRIVIAL(info) << "Read " << hla.sample_names.size() << " alleles.";

#ifndef NDEBUG
    // Save graph in debug mode
    save_graph(out_dir + "/graph");
#endif // NDEBUG

    PHIndex ph_index = index_graph(gyper::graph);

    BOOST_LOG_TRIVIAL(info) << "Finished indexing graph.";

    std::string reference_fn{}; // empty by default
    std::map<std::pair<uint16_t, uint16_t>, std::map<std::pair<uint16_t, uint16_t>, int8_t> > ph;

    if (Options::const_instance()->force_use_input_ref_for_cram_reading)
      reference_fn = ref_path;

    std::vector<std::string> paths;

    paths = gyper::call(shrinked_sams,
                        avg_cov_by_readlen,
                        "", // graph_path
                        ph_index,
                        out_dir,
                        "", // reference
                        ".", // region
                        nullptr,
                        ph,
                        is_writing_calls_vcf,
                        is_writing_hap,
                        &allele_hap_gts);

    BOOST_LOG_TRIVIAL(info) << "Merging output VCFs.";

    //for (auto & path : paths)
    //  path += "_calls.vcf.gz";

    // VCF merge
    {
      // Append _calls.vcf.gz
      //for (auto & path : paths)
      //  path += "_calls.vcf.gz";

      //> FILTER_ZERO_QUAL, force_no_variant_overlapping
      //vcf_merge_and_break(paths, tmp + "/graphtyper.vcf.gz", genomic_region.to_string(), false, false, true);
    }

    if (copts.force_ignore_segment)
    {
      vcf_merge_and_break(paths,
                          tmp + "/graphtyper.vcf.gz",
                          genomic_region.to_string(),
                          false,
                          false,
                          true);
    }
    else
    {
      Vcf hla_vcf;
      vcf_merge_and_return(hla_vcf, paths, tmp + "/graphtyper.vcf.gz");
      assert(hla_vcf.variants.size() == 1);
      Variant & var = hla_vcf.variants[0];

      {
        // Fix sequences names
        auto const is_pass_alt = var.generate_infos();
        assert(var.seqs.size() == hla.sample_names.size());
        assert((is_pass_alt.size() + 1u) == hla.sample_names.size());
        std::vector<std::vector<char> > new_seqs;

        for (long s{0}; s < static_cast<long>(hla.sample_names.size()); ++s)
        {
          if (s > 0 && is_pass_alt[s - 1] == 0)
          {
            //BOOST_LOG_TRIVIAL(info) << __HERE__ << " " << s << " " << static_cast<long>(is_pass_alt[s - 1]);
            continue;
          }

          std::string const & hla_allele = hla.sample_names[s];
          std::vector<char> seq; // = var.seqs[s];
          seq.push_back('<');
          std::copy(hla_allele.begin(), hla_allele.end(), std::back_inserter(seq));
          seq.push_back('>');
          new_seqs.push_back(seq);
        }

        // TODO Remove unused HLA alleles
        //if (new_seqs.size() < var.seqs.size())
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

              //BOOST_LOG_TRIVIAL(info) << __HERE__ << " " << x << "," << y;
              long const i = to_index(x, y);
              new_phred.push_back(call.phred[i]);
            }
          }

          //BOOST_LOG_TRIVIAL(warning) << __HERE__ << " " << new_phred.size() << " " << ((cnum * (cnum + 1)) / 2);
          assert(static_cast<long>(new_phred.size()) == ((cnum * (cnum + 1)) / 2));
          call.phred = new_phred;
        }

        var.seqs = new_seqs;
        var.stats = VarStats();
        var.generate_infos();
      }

      hla_vcf.open_for_writing(Options::const_instance()->threads);
      hla_vcf.write_header();

      if (hla_vcf.variants[0].seqs.size() < 100)
        hla_vcf.write_record(hla_vcf.variants[0], ".all", false, false);

      // Add 4 digit variant
      {
        std::unordered_map<std::string, long> seen_alleles;
        std::vector<long> old2new(var.seqs.size(), 0);
        Variant new_var;
        new_var.abs_pos = var.abs_pos;

        for (long a{0}; a < static_cast<long>(var.seqs.size()); ++a)
        {
          auto const & seq = var.seqs[a];
          auto find_it = std::find(seq.begin(), seq.end(), ':');
          assert(find_it != seq.end());

          if (find_it != seq.end())
            find_it = std::find(std::next(find_it), seq.end(), ':');

          std::string new_allele(seq.begin(), find_it);

          if (new_allele.back() != '>')
            new_allele.push_back('>');

          auto find_allele_it = seen_alleles.find(new_allele);

          if (find_allele_it == seen_alleles.end())
          {
            seen_alleles[new_allele] = new_var.seqs.size();
            old2new[a] = new_var.seqs.size();
            new_var.seqs.push_back(std::vector<char>(new_allele.begin(), new_allele.end()));
          }
          else
          {
            old2new[a] = find_allele_it->second;
          }
        }

        assert(seen_alleles.size() == new_var.seqs.size());

        if (new_var.seqs.size() > 1)
        {
          new_var.calls.reserve(var.calls.size());

          // fix calls
          for (long s{0}; s < static_cast<long>(var.calls.size()); ++s)
          {
            SampleCall new_call = bin_phred(new_var, var, var.calls[s], old2new);
            new_var.calls.push_back(std::move(new_call));
          }

          new_var.generate_infos();

          if (new_var.seqs.size() < 100)
            hla_vcf.write_record(new_var, ".4digit", false, false);
          //hla_vcf.variants.push_back(std::move(new_var));
        }
      }

      // Add 2 digit variant
      {
        std::unordered_map<std::string, long> seen_alleles;
        std::vector<long> old2new(var.seqs.size(), 0);
        Variant new_var;
        new_var.abs_pos = var.abs_pos;

        for (long a{0}; a < static_cast<long>(var.seqs.size()); ++a)
        {
          auto const & seq = var.seqs[a];
          auto find_it = std::find(seq.begin(), seq.end(), ':');

          std::string new_allele(seq.begin(), find_it);

          if (new_allele.back() != '>')
            new_allele.push_back('>');

          auto find_allele_it = seen_alleles.find(new_allele);

          if (find_allele_it == seen_alleles.end())
          {
            seen_alleles[new_allele] = new_var.seqs.size();
            old2new[a] = new_var.seqs.size();
            new_var.seqs.push_back(std::vector<char>(new_allele.begin(), new_allele.end()));
          }
          else
          {
            old2new[a] = find_allele_it->second;
          }
        }

        assert(seen_alleles.size() == new_var.seqs.size());

        if (new_var.seqs.size() > 1)
        {
          new_var.calls.reserve(var.calls.size());

          // fix calls
          for (long s{0}; s < static_cast<long>(var.calls.size()); ++s)
          {
            SampleCall new_call = bin_phred(new_var, var, var.calls[s], old2new);
            new_var.calls.push_back(std::move(new_call));
          }

          new_var.generate_infos();
          hla_vcf.write_record(new_var, ".2digit", false, false);
          //hla_vcf.variants.push_back(std::move(new_var));
        }
      }

      hla_vcf.close_vcf_file();
      hla_vcf.write_tbi_index();
    }
  }

  // gather segments from VCF

  // Copy final VCFs
  auto copy_vcf_to_system =
    [&](std::string const & name, std::string const & extension) -> void
    {
        filesystem::path src = tmp + "/" + name + extension;
        filesystem::path dest = output_path + "/" + genomic_region.to_file_string() +
          ((name == "graphtyper.hla.vcf.gz") ? ".hla" : "") + ".vcf.gz" + extension;

        filesystem::copy_file(src, dest, filesystem::copy_options::overwrite_existing);
    };

  std::string const index_ext = copts.is_csi ? ".csi" : ".tbi";

  copy_vcf_to_system("graphtyper.vcf.gz", ""); // Copy final VCF
  copy_vcf_to_system("graphtyper.vcf.gz", index_ext); // Copy tabix index for final VCF
  //copy_vcf_to_system("graphtyper.hla.vcf.gz", ""); // Copy final HLA VCF

  if (!copts.no_cleanup)
  {
    BOOST_LOG_TRIVIAL(info) << "Cleaning up temporary files.";
    remove_file_tree(tmp.c_str());
  }
  else
  {
    BOOST_LOG_TRIVIAL(info) << "Temporary files left: " << tmp;
  }

  {
    std::ostringstream ss;

    ss << output_path << "/" << genomic_region.chr << "/"
       << std::setw(9) << std::setfill('0') << (genomic_region.begin + 1)
       << '-'
       << std::setw(9) << std::setfill('0') << genomic_region.end
       << ".vcf.gz";

    BOOST_LOG_TRIVIAL(info) << "Finished! Output written at: " << ss.str();
  }

  // free memory
  graph = Graph();
}


void
genotype_hla_regions(std::string ref_path,
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
