#include <climits>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>

#include <graphtyper/graph/constructor.hpp>
#include <graphtyper/graph/graph_serialization.hpp>
#include <graphtyper/index/indexer.hpp>
#include <graphtyper/index/ph_index.hpp>
#include <graphtyper/utilities/type_conversions.hpp>

#include "../help_functions.hpp" // create_test_graph
#include <catch2/catch.hpp>

TEST_CASE("Test index chr1")
{
  gyper::print_log(gyper::log_severity::debug, "[", __HERE__, "] Test index chr1");

  using namespace gyper;

  create_test_graph("/test/data/reference/index_test.fa", "/test/data/reference/index_test.vcf.gz", "chr1", true);

  std::stringstream my_graph;
  my_graph << gyper_SOURCE_DIRECTORY << "/test/data/graphs/index_test_chr1.grf";

  REQUIRE(graph.size() > 0);

  // All reference kmers are present in the graph
  {
    std::vector<char> ref = graph.get_all_ref();
    REQUIRE(ref == gyper::to_vec("AGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCCCCTTTGGA"));
  }

  gyper::PHIndex ph_index = gyper::index_graph(my_graph.str());
  REQUIRE(ph_index.check());

  // Test one var kmer are present in the graph
  {
    std::vector<char> ref = graph.get_first_var();
    REQUIRE(ref == gyper::to_vec("AGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTTTGCCCAGGTTTCCCCAGGTTTCCCCTTTGGA"));
  }

  // All kmers should be have the correct count
  {
    REQUIRE(ph_index.get(to_uint64("AGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAG")).size() == 3);
    REQUIRE(ph_index.get(to_uint64("AGGTTTCCCCAGGTTTCCCCAGGTTTCCCCTT")).size() == 1);
    REQUIRE(ph_index.get(to_uint64("TTCCCCAGGTTTCCCCAGGTTTCCCCTTTGGA")).size() == 1);
    REQUIRE(ph_index.get(to_uint64("GGTTTCCCCAGGTTTCCCCAGGTTTGCCCAGG")).size() == 1);
  }

  // All kmers should have the correct starting index
  {
    REQUIRE(ph_index.get(to_uint64("AGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAG"))[0].start_index == 1);
    REQUIRE(ph_index.get(to_uint64("AGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAG"))[1].start_index == 11);
    REQUIRE(ph_index.get(to_uint64("AGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAG"))[2].start_index == 21);
    REQUIRE(ph_index.get(to_uint64("AGGTTTCCCCAGGTTTCCCCAGGTTTCCCCTT"))[0].start_index == 31);
    REQUIRE(ph_index.get(to_uint64("GGTTTCCCCAGGTTTCCCCAGGTTTGCCCAGG"))[0].start_index == 12);
  }

  // All kmers should have the correct end index
  {
    REQUIRE(ph_index.get(to_uint64("AGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAG"))[0].end_index == 32);
    REQUIRE(ph_index.get(to_uint64("AGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAG"))[1].end_index == 42);
    REQUIRE(ph_index.get(to_uint64("AGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAG"))[2].end_index == 52);
    REQUIRE(ph_index.get(to_uint64("AGGTTTCCCCAGGTTTCCCCAGGTTTCCCCTT"))[0].end_index == 62);
    REQUIRE(ph_index.get(to_uint64("GGTTTCCCCAGGTTTCCCCAGGTTTGCCCAGG"))[0].end_index == 43);
  }

  // Correct variant id
  {
    // AGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCCCCTTTGGA
    // chr1 37  rs1 C G 0 . .
    REQUIRE(ph_index.get(to_uint64("AGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAG"))[0].variant_id == INVALID_ID);
    REQUIRE(ph_index.get(to_uint64("AGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAG"))[1].variant_id == 0);
    REQUIRE(ph_index.get(to_uint64("AGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAG"))[2].variant_id == 0);
    REQUIRE(ph_index.get(to_uint64("AGGTTTCCCCAGGTTTCCCCAGGTTTCCCCTT"))[0].variant_id == 0);
    REQUIRE(ph_index.get(to_uint64("GGTTTCCCCAGGTTTCCCCAGGTTTGCCCAGG"))[0].variant_id == 1);
  }
}

TEST_CASE("Test index chr2")
{
  gyper::print_log(gyper::log_severity::debug, "[", __HERE__, "] Test index chr2");

  using namespace gyper;

  create_test_graph("/test/data/reference/index_test.fa", "/test/data/reference/index_test.vcf.gz", "chr2", true);

  std::stringstream my_graph;
  my_graph << gyper_SOURCE_DIRECTORY << "/test/data/graphs/index_test_chr2.grf";

  REQUIRE(graph.size() > 0);

  // The reference is present in the graph
  {
    std::vector<char> ref = graph.get_all_ref();
    REQUIRE(ref == gyper::to_vec("CCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTTTGGACCC"));
  }

  gyper::PHIndex ph_index = gyper::index_graph(graph);
  REQUIRE(ph_index.check());

  // All kmers should be have the correct count
  {
    REQUIRE(ph_index.get(to_uint64("CCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCC")).size() == 4);
    REQUIRE(ph_index.get(to_uint64("CCCCAGGTTTCCCCAGGTTTCCCCAGGTTTGG")).size() == 1);
    REQUIRE(ph_index.get(to_uint64("CACCAGGTTTCCCCAGGTTTCCCCAGGTTTCC")).size() == 2);
    REQUIRE(ph_index.get(to_uint64("CCACAGGTTTCCCCAGGTTTCCCCAGGTTTCC")).size() == 2);
    REQUIRE(ph_index.get(to_uint64("CAACAGGTTTCCCCAGGTTTCCCCAGGTTTCC")).size() == 2);
  }

  // All kmers should have the correct starting index
  {
    REQUIRE(ph_index.get(to_uint64("CCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCC"))[0].start_index == 1);
    REQUIRE(ph_index.get(to_uint64("CCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCC"))[1].start_index == 1);
    REQUIRE(ph_index.get(to_uint64("CCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCC"))[2].start_index == 11);
    REQUIRE(ph_index.get(to_uint64("CCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCC"))[3].start_index == 21);

    REQUIRE(ph_index.get(to_uint64("CCCCAGGTTTCCCCAGGTTTCCCCAGGTTTGG"))[0].start_index == 31);
  }

  // All kmers should have the correct end index
  {
    REQUIRE(ph_index.get(to_uint64("CCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCC"))[0].end_index == 32);
    REQUIRE(ph_index.get(to_uint64("CCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCC"))[1].end_index == 32);
    REQUIRE(ph_index.get(to_uint64("CCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCC"))[2].end_index == 42);
    REQUIRE(ph_index.get(to_uint64("CCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCC"))[3].end_index == 52);

    REQUIRE(ph_index.get(to_uint64("CCCCAGGTTTCCCCAGGTTTCCCCAGGTTTGG"))[0].end_index == 62);
  }

  // Correct variant id and num
  {
    using gyper::INVALID_ID;
    using gyper::to_uint64;

    // CCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCCCCAGGTTTGGACCC
    // chr2 2 rs2 C A 0 . .
    // chr2  3 rs3 C A 0 . .
    REQUIRE(ph_index.get(to_uint64("CCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCC"))[0].variant_id == 0);
    REQUIRE(ph_index.get(to_uint64("CCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCC"))[1].variant_id == 2);
    REQUIRE(ph_index.get(to_uint64("CCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCC"))[2].variant_id == INVALID_ID);
    REQUIRE(ph_index.get(to_uint64("CCCCAGGTTTCCCCAGGTTTCCCCAGGTTTCC"))[3].variant_id == INVALID_ID);
    REQUIRE(ph_index.get(to_uint64("CCCCAGGTTTCCCCAGGTTTCCCCAGGTTTGG"))[0].variant_id == INVALID_ID);
  }
}

TEST_CASE("Test index chr3")
{
  gyper::print_log(gyper::log_severity::debug, __HERE__, " Test index chr3");

  // AAAACAAAATAAAACAAAATAAAAGAAAACAAAATAAAACAAAATAAAAGAAAACATTATAAAACA
  // chr3 31 rs4 A G,GA

  using namespace gyper;

  create_test_graph("/test/data/reference/index_test.fa", "/test/data/reference/index_test.vcf.gz", "chr3", true);

  std::stringstream my_graph;
  my_graph << gyper_SOURCE_DIRECTORY << "/test/data/graphs/index_test_chr3.grf";
  REQUIRE(graph.size() > 0);

  gyper::PHIndex ph_index = gyper::index_graph(my_graph.str());
  REQUIRE(ph_index.check());

  // All reference kmers are present in the graph
  {
    std::vector<char> ref = graph.get_all_ref();
    REQUIRE(ref == gyper::to_vec("AAAACAAAATAAAACAAAATAAAAGAAAACAAAATAAAACAAAATAAAAGAAAACATTATAAAACA"));
  }

  // Kmer that starts at zero is has matches positions with the next one (although this is not a real alignmment)
  {
    // kmer ref
    std::vector<gyper::KmerLabel> labels0 = ph_index.get(gyper::to_uint64("AAAACAAAATAAAACAAAATAAAAGAAAACAA"));
    REQUIRE(labels0.size() == 1);
    REQUIRE(labels0[0].start_index == 1);
    REQUIRE(labels0[0].end_index == 32);

    // kmer 1
    std::vector<gyper::KmerLabel> labels1 = ph_index.get(gyper::to_uint64("AAAACAAAATAAAACAAAATAAAAGAAAACGA"));
    REQUIRE(labels1.size() == 2);
    REQUIRE(labels1[0].start_index == 1);
    REQUIRE(labels1[0].end_index == gyper::SPECIAL_START);
    REQUIRE(labels1[0].variant_id == 2);
    REQUIRE(labels1[1].start_index == 1);
    REQUIRE(labels1[1].end_index == 32);
    REQUIRE(labels1[1].variant_id == 1);

    // kmer 2
    std::vector<gyper::KmerLabel> labels2 = ph_index.get(gyper::to_uint64("AAAATAAAACAAAATAAAAGAAAACATTATAA"));
    REQUIRE(labels2.size() == 2);
    REQUIRE(labels2[0].start_index == 31);
    REQUIRE(labels2[0].end_index == 62);
    REQUIRE(labels2[0].variant_id == 0);
    REQUIRE(labels2[1].start_index == gyper::SPECIAL_START);
    REQUIRE(labels2[1].end_index == 62);
    REQUIRE(labels2[1].variant_id == 2);

    // kmer 3
    std::vector<gyper::KmerLabel> labels3 = ph_index.get(gyper::to_uint64("AAATAAAACAAAATAAAAGAAAACATTATAAA"));
    REQUIRE(labels3.size() == 1);
    REQUIRE(labels3[0].start_index == 32);
    REQUIRE(labels3[0].end_index == 63);
    REQUIRE(labels3[0].variant_id == -1);
  }
}

TEST_CASE("Test index chr4")
{
  gyper::print_log(gyper::log_severity::debug, __HERE__, " Test index chr4");

  using namespace gyper;

  create_test_graph("/test/data/reference/index_test.fa", "/test/data/reference/index_test.vcf.gz", "chr4", true);

  REQUIRE(graph.size() > 0);

  // The reference is present in the graph
  {
    std::vector<char> ref = graph.get_all_ref();
    REQUIRE(ref == gyper::to_vec("AAAACAAAATAAAACAAAATAAAAGAAAACAAAATAAAACAAAATAANNNNNNNNNNNNNNNNNNN"));
  }

  gyper::PHIndex ph_index = gyper::index_graph(graph);
  REQUIRE(ph_index.check());

  // kmer ref
  {
    std::vector<gyper::KmerLabel> labels0 = ph_index.get(to_uint64("AAAACAAAATAAAACAAAATAAAAGAAAACAA", 0));
    REQUIRE(labels0.size() == 1);
    REQUIRE(labels0[0].start_index == 1);
    REQUIRE(labels0[0].end_index == 32);
    REQUIRE(labels0[0].variant_id == 0);

    std::vector<gyper::KmerLabel> labels1 = ph_index.get(to_uint64("ATAACAAAATAAAACAAAATAAAAGAAAACAA", 0));
    REQUIRE(labels1.size() == 1);
    REQUIRE(labels1[0].start_index == 1);
    REQUIRE(labels1[0].end_index == 32);
    REQUIRE(labels1[0].variant_id == 1);
  }
}

TEST_CASE("Test index chr5")
{
  gyper::print_log(gyper::log_severity::debug, __HERE__, " Test index chr5");

  using namespace gyper;

  // 70A 70C 70G 70T
  // chr5 70A70C SVTYPE=DEL,SVSIZE=70

  create_test_graph("/test/data/reference/index_test.fa", "/test/data/reference/index_test.vcf.gz", "chr5", true);

  std::stringstream my_graph;
  my_graph << gyper_SOURCE_DIRECTORY << "/test/data/graphs/index_test_chr5.grf";

  REQUIRE(graph.size() > 0);

  PHIndex ph_index = gyper::index_graph(my_graph.str());
  REQUIRE(ph_index.check());

  // All reference kmers are present in the graph
  {
    std::vector<char> ref = graph.get_all_ref();
    REQUIRE(ref == gyper::to_vec("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                                 "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
                                 "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
                                 "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
  }

  // k-mer that starts at zero is has matches positions with the next one (although this is not a real alignmment)
  {
    // kmer ref
    std::vector<KmerLabel> labels0 = ph_index.get(to_uint64("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 0));
    REQUIRE(labels0.size() == 40);

    std::vector<KmerLabel> labels1 = ph_index.get(to_uint64("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAG", 0));
    REQUIRE(labels1.size() == 1);
    REQUIRE(labels1[0].start_index == 40);
    REQUIRE(labels1[0].end_index == gyper::SPECIAL_START);

    std::vector<KmerLabel> labels2 = ph_index.get(to_uint64("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGG", 0));
    REQUIRE(labels2.size() == 1);
    REQUIRE(labels2[0].start_index == 41);
    REQUIRE(labels2[0].end_index == gyper::SPECIAL_START + 1);

    std::vector<KmerLabel> labels3 = ph_index.get(to_uint64("AGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG", 0));
    REQUIRE(labels3.size() == 1);
    REQUIRE(labels3[0].start_index == 70);
    REQUIRE(labels3[0].end_index == gyper::SPECIAL_START + 30);

    std::vector<KmerLabel> labels4 = ph_index.get(to_uint64("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG", 0));
    REQUIRE(labels4.size() == 2 * (71 - K)); // 78

    // There should be strictly one starting at gyper::SPECIAL_START
    {
      unsigned found_count = 0;

      for (auto const & label : labels4)
        found_count += (label.start_index == gyper::SPECIAL_START + 1);

      REQUIRE(found_count == 1);
    }

    std::vector<KmerLabel> labels6 = ph_index.get(to_uint64("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT", 0));
    REQUIRE(labels6.size() == 2 * (71 - K)); // 78
  }
}

/*
TEST_CASE("Test index chr6")
{
  using namespace gyper;

  std::stringstream my_graph;
  my_graph << gyper_SOURCE_DIRECTORY << "/test/data/graphs/index_test_chr6.grf";

  gyper::load_graph(my_graph.str().c_str());
  REQUIRE(graph.size() > 0);

  gyper::index_graph(my_graph.str());
  REQUIRE(gyper::index.check());

  // All reference kmers are present in the graph
  {
    std::vector<char> ref = graph.get_all_ref();
    REQUIRE(ref == gyper::to_vec(
              "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
              "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
              "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
              "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
              ));
  }

  {
    // kmer ref
    std::vector<gyper::KmerLabel> labels0 = gyper::index.get(gyper::to_uint64("CGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG", 0));
    REQUIRE(labels0.size() == 3);
    REQUIRE(labels0[0].start_index >= gyper::SPECIAL_START);
    REQUIRE(labels0[0].end_index >= gyper::SPECIAL_START);
    REQUIRE(labels0[1].start_index == 139);
    REQUIRE(labels0[1].end_index == 170);
    REQUIRE(labels0[2].start_index >= gyper::SPECIAL_START);
    REQUIRE(labels0[2].end_index >= gyper::SPECIAL_START);
  }
}
*/

TEST_CASE("Test index chr9 with anti event")
{
  using namespace gyper;

  std::stringstream my_graph;
  my_graph << gyper_SOURCE_DIRECTORY << "/test/data/graphs/index_test_chr9.grf";

  gyper::load_graph(my_graph.str().c_str());
  REQUIRE(graph.size() > 0);

  gyper::PHIndex ph_index = gyper::index_graph(my_graph.str());
  REQUIRE(ph_index.check());

  std::vector<gyper::KmerLabel> labels = ph_index.get(gyper::to_uint64("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG", 0));
  REQUIRE(labels.size() == 36);

  labels = ph_index.get(gyper::to_uint64("GGGGGAGTGGGGGGGGGGGGGGGGGGGGGGGG", 0));
  REQUIRE(labels.size() == 1);
  REQUIRE(labels[0].variant_id == 3); // insertion

  labels = ph_index.get(gyper::to_uint64("GGGGGGGTGGGGGGGGGGGGGGGGGGGGGGGG", 0));

  for (auto const & label : labels)
  {
    bool const is_either = label.variant_id == 0 || label.variant_id == 2;
    REQUIRE(is_either);
  }

  REQUIRE(labels.size() == 2);
  REQUIRE(labels[0].variant_id != labels[1].variant_id);

  labels = ph_index.get(gyper::to_uint64("AGGGGGGTGGGGGGGGGGGGGGGGGGGGGGGG", 0));

  for (auto const & label : labels)
  {
    bool const is_either = label.variant_id == 0 || label.variant_id == 2;
    std::cerr << __HERE__ << " " << label.to_string() << "\n";
  }

  REQUIRE(labels.size() == 2);
  REQUIRE(labels[0].variant_id != labels[1].variant_id);

  labels = ph_index.get(gyper::to_uint64("AGGGGAGTGGGGGGGGGGGGGGGGGGGGGGGG", 0));
  REQUIRE(labels.size() == 0);
}

TEST_CASE("Test index chr10 with parity event")
{
  using namespace gyper;

  std::stringstream my_graph;
  my_graph << gyper_SOURCE_DIRECTORY << "/test/data/graphs/index_test_chr10.grf";

  gyper::load_graph(my_graph.str().c_str());
  REQUIRE(graph.size() > 0);

  gyper::PHIndex ph_index = gyper::index_graph(my_graph.str());
  REQUIRE(ph_index.check());

  ph_index.print();

  std::vector<gyper::KmerLabel> labels = ph_index.get(gyper::to_uint64("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG", 0));
  REQUIRE(labels.size() == 36);

  labels = ph_index.get(gyper::to_uint64("GGGGGAGTGGGGGGGGGGGGGGGGGGGGGGGG", 0));
  REQUIRE(labels.size() == 1);
  REQUIRE(labels[0].variant_id == 3); // insertion

  labels = ph_index.get(gyper::to_uint64("GGGGGGGTGGGGGGGGGGGGGGGGGGGGGGGG", 0));

  for (auto const & label : labels)
  {
    bool const is_either = label.variant_id == 0 || label.variant_id == 2;
    REQUIRE(is_either);
  }

  REQUIRE(labels.size() == 2);
  REQUIRE(labels[0].variant_id != labels[1].variant_id);

  labels = ph_index.get(gyper::to_uint64("AGGGGGGTGGGGGGGGGGGGGGGGGGGGGGGG", 0));
  REQUIRE(labels.size() == 2);

  labels = ph_index.get(gyper::to_uint64("AGGGGGAGTGGGGGGGGGGGGGGGGGGGGGGG", 0));

  for (auto const & label : labels)
  {
    bool const is_either = label.variant_id == 1 || label.variant_id == 3;
    REQUIRE(is_either);
  }

  REQUIRE(labels.size() == 2);
  REQUIRE(labels[0].variant_id != labels[1].variant_id);
}
