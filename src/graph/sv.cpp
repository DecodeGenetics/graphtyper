#include <cstdint>
#include <numeric>
#include <string>
#include <sstream>
#include <unordered_set>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include <graphtyper/graph/graph.hpp>
#include <graphtyper/graph/sv.hpp>
#include <graphtyper/graph/absolute_position.hpp>
#include <graphtyper/typer/variant.hpp>
#include <graphtyper/utilities/options.hpp>
#include <graphtyper/utilities/graph_help_functions.hpp>

namespace gyper
{

std::string
SV::get_type() const
{
  switch (type)
  {
  case DEL:
    return "DEL";

  case DEL_ALU:
    return "DEL:ME:ALU";

  case DUP:
    return "DUP";

  case INS:
    return "INS";

  case INS_ALU:
    return "INS:ME:ALU";

  case INV:
    return "INV";

  case BND:
    return "BND";

  default:
    return "SV";

  } // switch
}


std::string
SV::get_allele() const
{
  std::ostringstream ss;

  ss << '<' << get_type() << ":SVSIZE=";

  if (size > 0)
    ss << size;
  else
    ss << (ins_seq_left.size() + ins_seq_right.size()) << "+";

  ss << '>';
  return ss.str();
}


std::string
SV::get_allele_with_model() const
{
  std::ostringstream ss;

  ss << '<' << get_type() << ":SVSIZE=";

  if (size > 0)
    ss << size;
  else
    ss << (ins_seq_left.size() + ins_seq_right.size()) << "+";

  if (model.size() > 0)
    ss << ':' << model;

  ss << '>';
  return ss.str();
}


template <typename Archive>
void
SV::serialize(Archive & ar, const unsigned int /*version*/)
{
  ar & chrom;
  ar & begin;
  ar & type;
  ar & length;
  ar & size;
  ar & end;
  ar & n_clusters;
  ar & num_merged_svs;
  ar & or_end;
  ar & or_start;
  ar & related_sv;
  ar & model;
  ar & old_variant_id;
  ar & inv_type;
  ar & seq;
  ar & hom_seq;
  ar & ins_seq;
  ar & ins_seq_left;
  ar & ins_seq_right;
  ar & original_alt;
}


/***************************
 * EXPLICIT INSTANTIATIONS *
 ***************************/

template void SV::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive &,
                                                             const unsigned int);

template void SV::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive &,
                                                             const unsigned int);


void
reformat_sv_vcf_records(std::vector<Variant> & variants, ReferenceDepth const & reference_depth)
{
  long const variants_original_size = variants.size();
  std::unordered_set<long> variant_ids_to_erase; // Index of variants to erase
  std::unordered_map<int32_t, int32_t> related_svs; // Map of all related SVs
  std::vector<Variant> new_vars; // Keep new variants here to keep var reference valid

  for (long v = 0; v < variants_original_size; ++v)
  {
    auto & var = variants[v];

    // sv_id == -1 indicates that the alternative allele is NOT an SV
    std::vector<long> sv_ids(var.seqs.size() - 1, -1l);
    assert(var.seqs.size() >= 2u);

    for (long a = 1; a < static_cast<long>(var.seqs.size()); ++a)
    {
      auto & seq = var.seqs[a];
      auto find_it = std::find(seq.cbegin(), seq.cend(), '<');

      if (std::distance(find_it, seq.cend()) > 11)
      {
        // It is an SV
        std::istringstream ss{std::string(find_it + 4, find_it + 11)};
        long sv_id;
        ss >> sv_id;

        // If we can't parse correctly the SV ID we ignore it
        assert(ss.eof());
        assert(sv_id < static_cast<long>(graph.SVs.size()));
        sv_ids[a - 1l] = sv_id;
      }
    }

    {
      bool const is_any_sv = std::find_if(sv_ids.cbegin(),
                                          sv_ids.cend(),
                                          [](long const id){
          return id != -1;
        }
                                          ) != sv_ids.cend();

      if (not is_any_sv)
        continue; // Nothing to do, there are no SVs here
    }

    /*
    auto merge_alt_info_lambda =
      [](Variant & new_var, std::string const & id, long const aa)
      {
#ifndef NDEBUG
        if (new_var.infos.count(id) == 0)
        {
          BOOST_LOG_TRIVIAL(error) << "[graphtyper::graph::sv] Unable to find ID " << id << ".";
          std::exit(1);
        }
#endif // NDEBUG

        std::vector<uint32_t> values = split_bias_to_numbers(new_var.infos[id]);
        assert(aa + 1l < static_cast<long>(values.size()));
        values[1] = values[aa + 1];
        values.resize(2);
        new_var.infos[id] = join_strand_bias(values);
      };
    */

    auto make_new_sv_var =
      [&](Variant const & old_var, long const aa) -> Variant
      {
        Variant new_var;
        new_var.abs_pos = old_var.abs_pos;
        new_var.seqs.reserve(2);
        new_var.seqs.push_back(old_var.seqs[0]);
        new_var.seqs.push_back(old_var.seqs[aa + 1]);
        new_var.infos = old_var.infos;
        new_var.stats = old_var.stats;

        // Set value at index 1 to the alternative allele stats and then shrink the vector with .resize()
        new_var.stats.per_allele[1] = old_var.stats.per_allele[aa + 1];
        new_var.stats.read_strand[1] = old_var.stats.read_strand[aa + 1];

        new_var.stats.per_allele.resize(2);
        new_var.stats.read_strand.resize(2);

        assert(new_var.stats.per_allele.size() == new_var.seqs.size());

        new_var.calls.reserve(old_var.calls.size());

        for (long i = 0; i < static_cast<long>(old_var.calls.size()); ++i)
        {
          auto & call = old_var.calls[i];
          new_var.calls.push_back(make_bi_allelic_call(call, aa));
        }

        //new_var.normalize();
        assert(sv_ids[aa] != -1);
        auto const & sv_of_new_var = graph.SVs[sv_ids[aa]];

        if (sv_of_new_var.n_clusters > 0)
          new_var.infos["NCLUSTERS"] = std::to_string(sv_of_new_var.n_clusters);

        if (sv_of_new_var.num_merged_svs > 0)
          new_var.infos["NUM_MERGED_SVS"] = std::to_string(sv_of_new_var.num_merged_svs);

        new_var.infos["SV_ID"] = std::to_string(sv_ids[aa]);

        if (sv_of_new_var.related_sv >= 0)
          new_var.infos["RELATED_SV_ID"] = std::to_string(sv_of_new_var.related_sv);

        // Move SV to its original position
        new_var.abs_pos = absolute_pos.get_absolute_position(sv_of_new_var.chrom, sv_of_new_var.begin);
        return new_var;
      };


    auto make_variant_with_combined_calls =
      [](Variant const & var1, Variant const & var2) -> Variant
      {
        Variant combined_var(var1);
        assert(var1.calls.size() == var2.calls.size());

        for (long i = 0; i < static_cast<long>(var1.calls.size()); ++i)
        {
          auto & combined_call = combined_var.calls[i];
          auto const & var2_call = var2.calls[i];

          std::pair<uint16_t, uint16_t> gt_call2 = var2_call.get_gt_call();
          std::pair<uint16_t, uint16_t> gt_call1 = combined_call.get_gt_call();

          long gq1 = var2_call.get_gq();
          long gq2 = combined_call.get_gq();
          long max_gq = gq1;
          long min_gq = gq2;

          // Get unique depth here since here it is call1
          uint32_t dp1 = combined_call.get_unique_depth();

          if (gq1 > gq2)
          {
            combined_call = var2_call;
            max_gq = gq1;
            min_gq = gq2;
          }

          if (var1.calls[i].filter > 0 && var2.calls[i].filter > 0)
          {
            combined_call.filter = 3;
          }
          else if (var1.calls[i].filter > 0)
          {
            combined_call.filter = var1.calls[i].filter;
          }
          else if (var2.calls[i].filter > 0)
          {
            combined_call.filter = var2.calls[i].filter;
          }
          else if (dp1 >= 10u && var2_call.get_unique_depth() >= 10u)
          {
            std::pair<uint16_t,
                      uint16_t> final_gt_call = combined_call.get_gt_call();
            long index = to_index(final_gt_call.first, final_gt_call.second);
            assert(index < static_cast<long>(var1.calls[i].phred.size()));
            assert(index < static_cast<long>(var2.calls[i].phred.size()));

            if ((final_gt_call == gt_call1 && final_gt_call == gt_call2) &&
                min_gq > 10)
            {
              combined_call.filter = 0; // PASS
            }
            else if (max_gq > 40 &&
                     (var1.calls[i].phred[index] + var2.calls[i].phred[index]) <= 20)
            {
              combined_call.filter = 0; // PASS
            }
            else if (max_gq > 30)
            {
              combined_call.filter = 1; // FAIL1
            }
            else
            {
              combined_call.filter = 2; // FAIL2
            }
          }
          else
          {
            combined_call.filter = 3; // FAIL3
          }
        }

        return combined_var;
      };

    auto add_sv_to_new_vars_vector =
      [&new_vars](Variant && var, SV const & sv, std::string const & model) -> void
      {
        /// Set genotyping model in alternative allele
        if (sv.type != BND && !model.empty())
        {
          std::vector<char> & an = var.seqs[1]; // an = allele name
          an[an.size() - 1] = ':';
          std::copy(model.begin(), model.end(), std::back_inserter(an));
          an.push_back('>');
        }
        else if (sv.type == BND)
        {
          var.seqs[1] = sv.original_alt;
        }

        var.infos["SVTYPE"] = sv.get_type();

        // Make sure we do not write an END which is smaller then the begin position - otherwise bcftools is not happy
        if (sv.end < sv.begin)
          var.infos["END"] = std::to_string(sv.begin);
        else
          var.infos["END"] = std::to_string(sv.end);

        if (sv.length != 0)
          var.infos["SVSIZE"] = std::to_string(sv.size);

        if (sv.length != 0)
          var.infos["SVLEN"] = std::to_string(sv.length);

        if (model.size() > 0)
          var.infos["SVMODEL"] = model;

        if (sv.or_start != -1)
          var.infos["ORSTART"] = std::to_string(sv.or_start);

        if (sv.or_start != -1)
          var.infos["OREND"] = std::to_string(sv.or_end);

        if (sv.seq.size() > 0)
          var.infos["SEQ"] = std::string(sv.seq.begin(), sv.seq.end());

        if (sv.n_clusters > 0)
          var.infos["NCLUSTERS"] = std::to_string(sv.n_clusters);

        if (sv.num_merged_svs >= 0)
          var.infos["NUM_MERGED_SVS"] = std::to_string(sv.num_merged_svs);

        if (sv.old_variant_id.size() > 0 && sv.old_variant_id != ".")
          var.infos["OLD_VARIANT_ID"] = sv.old_variant_id;

        if (sv.hom_seq.size() > 0)
          var.infos["HOMSEQ"] = std::string(sv.hom_seq.begin(), sv.hom_seq.end());

        if (sv.ins_seq.size() > 0)
          var.infos["SVINSSEQ"] = std::string(sv.ins_seq.begin(), sv.ins_seq.end());

        if (sv.ins_seq_left.size() > 0)
          var.infos["LEFT_SVINSSEQ"] = std::string(sv.ins_seq_left.begin(), sv.ins_seq_left.end());

        if (sv.ins_seq_right.size() > 0)
          var.infos["RIGHT_SVINSSEQ"] = std::string(sv.ins_seq_right.begin(), sv.ins_seq_right.end());

        if (sv.type == INV && sv.inv_type != NOT_INV)
        {
          switch (sv.inv_type)
          {
          case INV3:
            var.infos["INV3"] = "";
            break;

          case INV5:
            var.infos["INV5"] = "";
            break;

          case BOTH_BREAKPOINTS:
            var.infos["INV3"] = "";
            var.infos["INV5"] = "";
            break;

          default:
            break;
          }
        }


        new_vars.push_back(std::move(var));
      };


    /// Adds an SV variant to 'new_vars' vector
    auto add_sv_variant =
      [&](Variant const & var,
          long aa,
          std::unordered_map<int32_t, int32_t> & related_svs_map) -> void
      {
        // Put first SV allele in new VCF records
        Variant new_sv_var = make_new_sv_var(var, aa);
        assert(new_sv_var.seqs.size() == 2);
        assert(new_sv_var.stats.per_allele.size() == 2);
        assert(new_sv_var.stats.read_strand.size() == 2);
        SV const & sv = graph.SVs[sv_ids[aa]];

        if (sv.type != BND)
        {
          new_sv_var.seqs[0] = std::vector<char>(1, 'N');
          std::string sv_allele = sv.get_allele();
          new_sv_var.seqs[1] = std::vector<char>(begin(sv_allele), end(sv_allele));
        }

        /*// Fix position of first breakpoint for cases it is moved
        if ((sv.type == DUP && sv.related_sv != -1 && sv.model == "BREAKPOINT1") ||
            (sv.type == INV && sv.related_sv != -1 && sv.model == "BREAKPOINT2"))
        {
          new_sv_var.abs_pos -= sv.size;                     // Move back to original pos
          }*/

        // Fix genotype likelihood for duplication
        if (sv.type == DUP &&
            (sv.model == "BREAKPOINT1" || sv.model == "BREAKPOINT2"))
        {
          for (auto & call : new_sv_var.calls)
          {
            uint64_t constexpr ERROR = 25;
            double constexpr minus_10log10_one_third = 4.77121255;
            double constexpr minus_10log10_two_third = 1.76091259;

            uint64_t gt_00 = call.coverage[1] * ERROR;
            uint64_t gt_01 = static_cast<uint64_t>(0.499999999 +
                                                   minus_10log10_one_third *
                                                   static_cast<double>(call.coverage[1]) +
                                                   minus_10log10_two_third *
                                                   static_cast<double>(call.coverage[0]));

            uint64_t gt_11 = 3ul * (call.coverage[0] + static_cast<uint64_t>(call.coverage[1]));
            uint64_t min_gt = std::min(gt_00, std::min(gt_01, gt_11));

            gt_00 -= min_gt;
            gt_01 -= min_gt;
            gt_11 -= min_gt;

            call.phred[0] = static_cast<uint8_t>(std::min(static_cast<uint64_t>(0xFFu), gt_00));
            call.phred[1] = static_cast<uint8_t>(std::min(static_cast<uint64_t>(0xFFu), gt_01));
            call.phred[2] = static_cast<uint8_t>(std::min(static_cast<uint64_t>(0xFFu), gt_11));
          }
        }

        // If the new variant is an insertion and this is the second breakpoint, make a combined
        // variant. Also do the same if it is a duplication and we have coverage turned off
        if ((sv.type == INS && related_svs_map.count(sv_ids[aa]) == 1) ||
            (sv.type == INV && related_svs_map.count(sv_ids[aa]) == 1))
        {
          assert(new_vars.size() > 0);
          Variant const & var_bp1 = new_vars[related_svs_map.at(sv_ids[aa])];
          assert(var_bp1.seqs.size() == var_bp1.stats.per_allele.size());
          Variant combined_var = make_variant_with_combined_calls(new_sv_var, var_bp1);
          assert(combined_var.seqs.size() == combined_var.stats.per_allele.size());
          add_sv_to_new_vars_vector(std::move(combined_var), sv, "AGGREGATED");
        }

        // If the new variant is a deletion, check coverage
        if (graph.is_sv_graph)
        {
          if ((sv.type == DEL || sv.type == DEL_ALU))
          {
            Variant cov_var(new_sv_var);
            assert(cov_var.seqs.size() == 2);
            assert(cov_var.stats.per_allele.size() == 2);
            assert(cov_var.stats.read_strand.size() == 2);

            for (long pn_index = 0; pn_index < static_cast<long>(cov_var.calls.size()); ++pn_index)
              cov_var.calls[pn_index] = make_call_based_on_coverage(pn_index, sv, reference_depth);

            // Make a combined variant with both breakpoint and coverage
            Variant combined_var = make_variant_with_combined_calls(new_sv_var, cov_var);
            assert(combined_var.seqs.size() == 2);
            assert(combined_var.stats.per_allele.size() == 2);
            assert(combined_var.stats.read_strand.size() == 2);
            add_sv_to_new_vars_vector(std::move(combined_var), sv, "AGGREGATED");
            assert(cov_var.seqs.size() >= 2);
            add_sv_to_new_vars_vector(std::move(cov_var), sv, "COVERAGE");
          }
          else if (sv.type == DUP && related_svs_map.count(sv_ids[aa]) == 1)
          {
            // Second breakpoint of a duplication
            Variant cov_var(new_sv_var);
            assert(cov_var.seqs.size() == 2ul);
            assert(cov_var.seqs.size() == cov_var.stats.per_allele.size());

            for (long pn_index = 0; pn_index < static_cast<long>(cov_var.calls.size()); ++pn_index)
              cov_var.calls[pn_index] = make_call_based_on_coverage(pn_index, sv, reference_depth);

            // Make a combined variant with both breakpoint and coverage
            Variant combined_var = make_variant_with_combined_calls(new_sv_var, cov_var);
            Variant const & other_bp_variant = new_vars[related_svs.at(sv_ids[aa])];
            Variant combined_var2 =
              make_variant_with_combined_calls(combined_var, other_bp_variant);
            assert(combined_var.seqs.size() >= 2);
            add_sv_to_new_vars_vector(std::move(combined_var2), sv, "AGGREGATED");
            assert(cov_var.seqs.size() >= 2);
            add_sv_to_new_vars_vector(std::move(cov_var), sv, "COVERAGE");
          }
          else if (sv.type == BND)
          {
            if (new_sv_var.seqs[1][1] == '<')
            {
              new_sv_var.add_base_in_back();
              new_sv_var.remove_common_prefix();
            }
          }
        }

        // If the SV has a related SV, add it to the related_svs_map hashmap
        if (sv.related_sv != -1)
        {
          assert(related_svs_map.count(sv.related_sv) == 0);
          related_svs_map[static_cast<int>(sv.related_sv)] = static_cast<int>(new_vars.size());
          assert(related_svs_map.count(sv.related_sv) == 1);
        }

        assert(new_sv_var.seqs.size() >= 2);
        add_sv_to_new_vars_vector(std::move(new_sv_var), sv, sv.model);
      };

    assert(sv_ids.size() + 1u == var.seqs.size());
    //long aa = 0; // Index of alternative allele
    bool is_any_not_sv = false; //std::find(sv_ids.cbegin(), sv_ids.cend(), -1) != sv_ids.cend();

    for (long aa = 0; aa < static_cast<long>(sv_ids.size()); ++aa)
    {
      if (sv_ids[aa] == -1l)
      {
        is_any_not_sv = true;
        continue; // Not SV
      }

      assert(var.seqs.size() == var.stats.per_allele.size());
      add_sv_variant(var, aa, related_svs);
    }

    if (is_any_not_sv)
    {
      Variant non_sv_var;
      non_sv_var.abs_pos = var.abs_pos;
      non_sv_var.infos = var.infos;
//      non_sv_var.phase = var.phase;
      non_sv_var.suffix_id = var.suffix_id;

      // Copy reference everywhere except this SNP
      non_sv_var.seqs = std::vector<std::vector<char> >(var.seqs.size(), var.seqs[0]);

      for (long aa = 0; aa < static_cast<long>(sv_ids.size()); ++aa)
      {
        if (sv_ids[aa] == -1l)
          non_sv_var.seqs[aa + 1] = var.seqs[aa + 1];
      }

      find_variant_sequences(non_sv_var, var);
      non_sv_var.normalize();
      new_vars.push_back(std::move(non_sv_var));

      /*
      var.infos.clear();
      // Replace the old variant with a new one that contains no SV records
      std::vector<std::vector<char> > new_seqs;
      new_seqs.push_back(var.seqs[0]); // Add reference allele

      for (long aa = 0; aa < static_cast<long>(sv_ids.size()); ++aa)
      {
        if (sv_ids[aa] == -1l)
          new_seqs.push_back(var.seqs[aa + 1]);
      }

      assert(new_seqs.size() >= 2u);
      var.seqs = std::move(new_seqs);

      for (auto & call : var.calls)
      {
        call.ambiguous_depth = 0;
        call.ref_total_depth = 0;
        call.alt_total_depth = 0;
        call.alt_proper_pair_depth = 0;

        assert(call.coverage.size() >= 2);
        assert(sv_ids.size() + 1u == call.coverage.size());

        std::vector<uint16_t> new_coverage;
        new_coverage.push_back(call.coverage[0]);

        for (long aa = 0; aa < static_cast<long>(sv_ids.size()); ++aa)
        {
          if (sv_ids[aa] == -1l)
            new_coverage.push_back(call.coverage[aa + 1]);
        }

        long const cnum = var.seqs.size();
        assert(new_coverage.size() == cnum);
        long const total_coverage = std::accumulate(begin(new_coverage), end(new_coverage), 0l);
        std::vector<long> phreds(cnum * (cnum + 1) / 2, 255);
        long constexpr ERROR_PHRED = 25;

        for (long y = 0; y < cnum; ++y)
        {
          for (long x = 0; x <= y; ++x)
          {
            assert(to_index(x, y) < (long)phreds.size());
            assert(x < (long)new_coverage.size());
            assert(y < (long)new_coverage.size());

            if (x == y) // homozygous
            {
              assert(new_coverage[x] <= total_coverage);
              phreds[to_index(x, y)] = ERROR_PHRED * (total_coverage - new_coverage[x]);
            }
            else // heterozygous
            {
              uint64_t const allele_coverage = new_coverage[x] + new_coverage[y];
              assert(allele_coverage <= total_coverage);
              phreds[to_index(x, y)] = 3l * allele_coverage + ERROR_PHRED * (total_coverage - allele_coverage);
            }
          }
        }

        assert(phreds.size() > 0);
        long min_val = phreds.front();

        for (auto const val : phreds)
          min_val = std::min(min_val, val);

        call.phred.resize(static_cast<std::size_t>(cnum + (cnum + 1) / 2), 0);

        for (long i = 0; i < static_cast<long>(call.phred.size()); ++i)
        {
          assert(i < (long)phreds.size());
          long const new_phred = std::max(0l, std::min(255l, phreds[i] - min_val));
          call.phred[i] = static_cast<uint8_t>(new_phred);
        }

        call.coverage = std::move(new_coverage);
      }
      */
    }

    variant_ids_to_erase.insert(v); // Remove the original variant
  }

  // Erase variants that should be erased
  if (variant_ids_to_erase.size() > 0u)
  {
    std::size_t const old_size = new_vars.size() + variants.size();
    new_vars.reserve(old_size - variant_ids_to_erase.size());

    for (long v = 0l; v < variants_original_size; ++v)
    {
      if (variant_ids_to_erase.count(v) == 0)
        new_vars.push_back(std::move(variants[v]));
    }

    assert(new_vars.size() == old_size - variant_ids_to_erase.size());
    variants = std::move(new_vars);
  }
}


} // namespace gyper
