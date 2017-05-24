#include <catch.hpp>

#include <stdio.h>
#include <limits.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>

#include <graphtyper/graph/graph_serialization.hpp>
#include <graphtyper/graph.hpp>
#include <graphtyper/options.hpp>
#include <graphtyper/constants.hpp>


using TConstruct = gyper::Constructor<gyper::SimpleGraph>;


TEST_CASE("SeqAn tabix support", "[seqanhts]")
{
  SECTION("BCF (binary vcf)")
  {
    SECTION("Same results as above with a compressed VCF file (BCF) (example_bgz.vcf.gz).")
    {
      std::stringstream base_path;
      base_path << gyper_SOURCE_DIRECTORY << "/test/reference/example_bgz.vcf.gz";
      std::string base_path_str = base_path.str();
      TConstruct * constructor = new TConstruct();
      constructor->open_tabix(base_path_str.c_str());

      SECTION("The getNextFunction")
      {
        // REQUIRE(constructor->read_tabix_record() == 0);
        seqan::String<char> line;
        REQUIRE(seqan::readRawRecord(line, constructor->tabix_file));
        REQUIRE(line == "chr1	2	rs6054257	G	A	0	.	.");
        REQUIRE(seqan::readRawRecord(line, constructor->tabix_file));
        REQUIRE(line == "chr1	4	ins504320	G	GT	0	.	.");
        REQUIRE(seqan::readRawRecord(line, constructor->tabix_file));
        REQUIRE(line == "chr1	6	del504320	TC	T	0	.	.");
        REQUIRE(seqan::readRawRecord(line, constructor->tabix_file));
        REQUIRE(line == "chr1	11	msat14123	GGT	G,GGTGTGT	0	.	.");
        REQUIRE(seqan::readRawRecord(line, constructor->tabix_file));
        REQUIRE(line == "chr2	3	rs0000001	A	T	0	.	.");
        REQUIRE(!seqan::readRawRecord(line, constructor->tabix_file));
      }

      SECTION("The readRecord method")
      {
        REQUIRE(constructor->read_tabix_record());
        REQUIRE(constructor->vcf_record.rID == 0);
        REQUIRE(constructor->vcf_record.beginPos == 1);
        REQUIRE(constructor->vcf_record.id == "rs6054257");
        REQUIRE(constructor->vcf_record.ref == "G");
        REQUIRE(constructor->vcf_record.alt == "A");
        REQUIRE(constructor->vcf_record.qual == 0.0);
        REQUIRE(constructor->vcf_record.filter == ".");
        REQUIRE(constructor->vcf_record.info == ".");
        REQUIRE(constructor->vcf_record.format == "");
        REQUIRE(seqan::length(constructor->vcf_record.genotypeInfos) == 0);
        seqan::String<seqan::Dna> ref_dna = "G";
        REQUIRE(static_cast<seqan::String<seqan::Dna> >(constructor->vcf_record.ref) == ref_dna);

        REQUIRE(constructor->read_tabix_record());
        REQUIRE(constructor->vcf_record.rID == 0);
        REQUIRE(constructor->vcf_record.id == "ins504320");
        REQUIRE(constructor->vcf_record.ref == "G");
        REQUIRE(constructor->vcf_record.alt == "GT");
        seqan::String<seqan::Dna> alt_dna = "GT";
        REQUIRE(static_cast<seqan::String<seqan::Dna> >(constructor->vcf_record.alt) == alt_dna);

        REQUIRE(constructor->read_tabix_record());
        REQUIRE(constructor->vcf_record.rID == 0);
        REQUIRE(constructor->vcf_record.id == "del504320");
        REQUIRE(constructor->vcf_record.ref == "TC");
        REQUIRE(constructor->vcf_record.alt == "T");

        REQUIRE(constructor->read_tabix_record());
        REQUIRE(constructor->vcf_record.rID == 0);
        REQUIRE(constructor->vcf_record.id == "msat14123");
        REQUIRE(constructor->vcf_record.ref == "GGT");
        REQUIRE(constructor->vcf_record.alt == "G,GGTGTGT");

        REQUIRE(constructor->read_tabix_record());
        REQUIRE(constructor->vcf_record.rID == 1);
        REQUIRE(constructor->vcf_record.id == "rs0000001");
        REQUIRE(constructor->vcf_record.ref == "A");
        REQUIRE(constructor->vcf_record.alt == "T");

        REQUIRE(!constructor->read_tabix_record());
      }

      SECTION("The readRegion method")
      {
        seqan::String<seqan::VcfRecord> records;

        SECTION("Using extractions of all variants in whole chromosomes")
        {
          // Read chr1
          clear(records);
          REQUIRE(seqan::length(records) == 0);
          readRegion(records, constructor->tabix_file, "chr1");
          REQUIRE(seqan::length(records) == 4);

          // Read chr2
          clear(records);
          REQUIRE(seqan::length(records) == 0);
          readRegion(records, constructor->tabix_file, "chr2");
          REQUIRE(seqan::length(records) == 1);

          // Read both chromosomes
          clear(records);
          REQUIRE(seqan::length(records) == 0);
          readRegion(records, constructor->tabix_file, "chr1");
          readRegion(records, constructor->tabix_file, "chr2");
          REQUIRE(seqan::length(records) == 5);
        }

        SECTION("Using extraction of partial chromosomes")
        {
          // Read chr1:2-2. This should only read the variant on position 2 of chr 1.
          clear(records);
          REQUIRE(seqan::length(records) == 0);
          readRegion(records, constructor->tabix_file, "chr1:2-2");
          REQUIRE(seqan::length(records) == 1);

          // Read chr1:2-3. This does not read the variant on position 4.
          clear(records);
          REQUIRE(seqan::length(records) == 0);
          readRegion(records, constructor->tabix_file, "chr1:2-3");
          REQUIRE(seqan::length(records) == 1);

          // Read chr1:2-3
          clear(records);
          REQUIRE(seqan::length(records) == 0);
          readRegion(records, constructor->tabix_file, "chr1:2-4");
          REQUIRE(seqan::length(records) == 2);

          // Read chr1:4-4
          clear(records);
          REQUIRE(seqan::length(records) == 0);
          readRegion(records, constructor->tabix_file, "chr1:4-4");
          REQUIRE(seqan::length(records) == 1);

          // Read chr1:4-6
          clear(records);
          REQUIRE(seqan::length(records) == 0);
          readRegion(records, constructor->tabix_file, "chr1:4-6");
          REQUIRE(seqan::length(records) == 2);

          // Read chr1:0-10000000000. This reads the entire chromosome 1
          clear(records);
          REQUIRE(seqan::length(records) == 0);
          readRegion(records, constructor->tabix_file, "chr1:1-10000000000");
          REQUIRE(seqan::length(records) == 4);

          // Read chr2:0-0.
          clear(records);
          REQUIRE(seqan::length(records) == 0);
          readRegion(records, constructor->tabix_file, "chr2:0-0");
          REQUIRE(seqan::length(records) == 0);

          // Read chr3:3-3.
          clear(records);
          REQUIRE(seqan::length(records) == 0);
          readRegion(records, constructor->tabix_file, "chr2:3-3");
          REQUIRE(seqan::length(records) == 1);
        }

        seqan::VcfRecord record;
        seqan::Tabix index = constructor->tabix_file;

        SECTION("Reading one record")
        {
          seqan::setRegion(index, "chr2:3-3");

          while (seqan::readRegion(record, index))
          {
            // Do stuff with record
          }
        }
      }

      delete constructor;
    }

    SECTION("Another BCF example with 2 variants on chr1 and 3 on chr2. (example2_bgz.vcf.gz)")
    {
      std::stringstream base_path;
      base_path << gyper_SOURCE_DIRECTORY << "/test/reference/example2_bgz.vcf.gz";
      std::string base_path_str = base_path.str();
      TConstruct * constructor = new TConstruct();
      constructor->open_tabix(base_path_str.c_str());

      SECTION("The getNextFunction")
      {
        seqan::String<char> line;
        REQUIRE(seqan::readRawRecord(line, constructor->tabix_file));
        REQUIRE(line == "chr1	2	rs6054257	G	A	0	.	.");
        REQUIRE(seqan::readRawRecord(line, constructor->tabix_file));
        REQUIRE(line == "chr1	4	ins504320	T	GT	0	.	.");
        REQUIRE(seqan::readRawRecord(line, constructor->tabix_file));
        REQUIRE(line == "chr2	6	del504320	T	TA	0	.	.");
        REQUIRE(seqan::readRawRecord(line, constructor->tabix_file));
        REQUIRE(line == "chr2	13	del0000001	AG	A	0	.	.");
        REQUIRE(seqan::readRawRecord(line, constructor->tabix_file));
        REQUIRE(line == "chr2	14	mix14123	G	GTT,T	0	.	.");
        REQUIRE(!seqan::readRawRecord(line, constructor->tabix_file));
      }

      SECTION("The readRecord method")
      {
        REQUIRE(constructor->read_tabix_record());
        REQUIRE(constructor->vcf_record.rID == 0);
        REQUIRE(constructor->vcf_record.beginPos == 1);
        REQUIRE(constructor->vcf_record.id == "rs6054257");
        REQUIRE(constructor->vcf_record.ref == "G");
        REQUIRE(constructor->vcf_record.alt == "A");
        REQUIRE(constructor->vcf_record.qual == 0.0);
        REQUIRE(constructor->vcf_record.filter == ".");
        REQUIRE(constructor->vcf_record.info == ".");
        REQUIRE(constructor->vcf_record.format == "");
        REQUIRE(length(constructor->vcf_record.genotypeInfos) == 0);
        seqan::String<seqan::Dna> ref_dna = "G";
        REQUIRE(static_cast<seqan::String<seqan::Dna> >(constructor->vcf_record.ref) == ref_dna);

        REQUIRE(constructor->read_tabix_record());
        REQUIRE(constructor->vcf_record.rID == 0);
        REQUIRE(constructor->vcf_record.id == "ins504320");
        REQUIRE(constructor->vcf_record.ref == "T");
        REQUIRE(constructor->vcf_record.alt == "GT");

        REQUIRE(constructor->read_tabix_record());
        REQUIRE(constructor->vcf_record.rID == 1);
        REQUIRE(constructor->vcf_record.id == "del504320");
        REQUIRE(constructor->vcf_record.ref == "T");
        REQUIRE(constructor->vcf_record.alt == "TA");

        REQUIRE(constructor->read_tabix_record());
        REQUIRE(constructor->vcf_record.rID == 1);
        REQUIRE(constructor->vcf_record.id == "del0000001");
        REQUIRE(constructor->vcf_record.ref == "AG");
        REQUIRE(constructor->vcf_record.alt == "A");

        REQUIRE(constructor->read_tabix_record());
        REQUIRE(constructor->vcf_record.rID == 1);
        REQUIRE(constructor->vcf_record.id == "mix14123");
        REQUIRE(constructor->vcf_record.ref == "G");
        REQUIRE(constructor->vcf_record.alt == "GTT,T");

        REQUIRE(!constructor->read_tabix_record());
      }

      delete constructor;
    }

    SECTION("Another BCF example with 2 variants on chr1. (example3.vcf.gz)")
    {
      std::stringstream base_path;
      base_path << gyper_SOURCE_DIRECTORY << "/test/reference/example3.vcf.gz";
      std::string base_path_str = base_path.str();
      TConstruct * constructor = new TConstruct();
      constructor->open_tabix(base_path_str.c_str());

      SECTION("The 'readRawRecord' method")
      {
        seqan::String<char> line;
        REQUIRE(seqan::readRawRecord(line, constructor->tabix_file));
        REQUIRE(line == "chr1	2	rs6054257	G	A	0	.	.");
        REQUIRE(seqan::readRawRecord(line, constructor->tabix_file));
        REQUIRE(line == "chr1	4	ins504320	T	GT	0	.	.");
        REQUIRE(!seqan::readRawRecord(line, constructor->tabix_file));
      }

      SECTION("The readRecord method")
      {
        REQUIRE(constructor->read_tabix_record());
        REQUIRE(constructor->vcf_record.rID == 0);
        REQUIRE(constructor->vcf_record.beginPos == 1);
        REQUIRE(constructor->vcf_record.id == "rs6054257");
        REQUIRE(constructor->vcf_record.ref == "G");
        REQUIRE(constructor->vcf_record.alt == "A");
        REQUIRE(constructor->vcf_record.qual == 0.0);
        REQUIRE(constructor->vcf_record.filter == ".");
        REQUIRE(constructor->vcf_record.info == ".");
        REQUIRE(constructor->vcf_record.format == "");
        REQUIRE(length(constructor->vcf_record.genotypeInfos) == 0);
        seqan::String<seqan::Dna> ref_dna = "G";
        REQUIRE(static_cast<seqan::String<seqan::Dna> >(constructor->vcf_record.ref) == ref_dna);

        REQUIRE(constructor->read_tabix_record());
        REQUIRE(constructor->vcf_record.rID == 0);
        REQUIRE(constructor->vcf_record.id == "ins504320");
        REQUIRE(constructor->vcf_record.ref == "T");
        REQUIRE(constructor->vcf_record.alt == "GT");

        REQUIRE(!constructor->read_tabix_record());
      }

      delete constructor;
    }


    SECTION("Another BCF example with and 3 on chr2. (example4.vcf.gz)")
    {
      std::stringstream base_path;
      base_path << gyper_SOURCE_DIRECTORY << "/test/reference/example4.vcf.gz";
      std::string base_path_str = base_path.str();
      TConstruct * constructor = new TConstruct();
      constructor->open_tabix(base_path_str.c_str());

      SECTION("The getNextFunction")
      {
        seqan::String<char> line;
        REQUIRE(seqan::readRawRecord(line, constructor->tabix_file));
        REQUIRE(line == "chr2	6	del504320	T	TA	0	.	.");
        REQUIRE(seqan::readRawRecord(line, constructor->tabix_file));
        REQUIRE(line == "chr2	13	del0000001	AG	A	0	.	.");
        REQUIRE(seqan::readRawRecord(line, constructor->tabix_file));
        REQUIRE(line == "chr2	14	mix14123	G	GTT,T	0	.	.");
        REQUIRE(!seqan::readRawRecord(line, constructor->tabix_file));
      }

      SECTION("The readRecord method")
      {
        REQUIRE(constructor->read_tabix_record());
        REQUIRE(constructor->vcf_record.rID == 0);
        REQUIRE(constructor->vcf_record.id == "del504320");
        REQUIRE(constructor->vcf_record.ref == "T");
        REQUIRE(constructor->vcf_record.alt == "TA");

        REQUIRE(constructor->read_tabix_record());
        REQUIRE(constructor->vcf_record.rID == 0);
        REQUIRE(constructor->vcf_record.id == "del0000001");
        REQUIRE(constructor->vcf_record.ref == "AG");
        REQUIRE(constructor->vcf_record.alt == "A");

        REQUIRE(constructor->read_tabix_record());
        REQUIRE(constructor->vcf_record.rID == 0);
        REQUIRE(constructor->vcf_record.id == "mix14123");
        REQUIRE(constructor->vcf_record.ref == "G");
        REQUIRE(constructor->vcf_record.alt == "GTT,T");

        REQUIRE(!constructor->read_tabix_record());
      }

      delete constructor;
    }

    SECTION("Another BCF example with variants on chr1 and chr3. (example5.vcf.gz)")
    {
      std::stringstream base_path;
      base_path << gyper_SOURCE_DIRECTORY << "/test/reference/example5.vcf.gz";
      std::string base_path_str = base_path.str();
      TConstruct * constructor = new TConstruct();
      constructor->open_tabix(base_path_str.c_str());

      SECTION("The getNextFunction")
      {
        seqan::String<char> line;
        REQUIRE(seqan::readRawRecord(line, constructor->tabix_file));
        REQUIRE(line == "chr1	6	del504320	T	TA	0	.	.");
        REQUIRE(seqan::readRawRecord(line, constructor->tabix_file));
        REQUIRE(line == "chr3	13	del0000001	AG	A	0	.	.");
        REQUIRE(seqan::readRawRecord(line, constructor->tabix_file));
        REQUIRE(line == "chr3	14	mix14123	G	GTT,T	0	.	.");
        REQUIRE(!seqan::readRawRecord(line, constructor->tabix_file));
      }

      SECTION("The readRecord method")
      {
        REQUIRE(constructor->read_tabix_record());
        REQUIRE(constructor->vcf_record.rID == 0);
        REQUIRE(constructor->vcf_record.beginPos == 5);
        REQUIRE(constructor->vcf_record.id == "del504320");
        REQUIRE(constructor->vcf_record.ref == "T");
        REQUIRE(constructor->vcf_record.alt == "TA");
        REQUIRE(constructor->vcf_record.qual == 0.0);
        REQUIRE(constructor->vcf_record.filter == ".");
        REQUIRE(constructor->vcf_record.info == ".");
        REQUIRE(constructor->vcf_record.format == "");
        REQUIRE(length(constructor->vcf_record.genotypeInfos) == 0);

        REQUIRE(constructor->read_tabix_record());
        REQUIRE(constructor->vcf_record.rID == 1);
        REQUIRE(constructor->vcf_record.beginPos == 12);
        REQUIRE(constructor->vcf_record.id == "del0000001");
        REQUIRE(constructor->vcf_record.ref == "AG");
        REQUIRE(constructor->vcf_record.alt == "A");

        REQUIRE(constructor->read_tabix_record());
        REQUIRE(constructor->vcf_record.rID == 1);
        REQUIRE(constructor->vcf_record.beginPos == 13);
        REQUIRE(constructor->vcf_record.id == "mix14123");
        REQUIRE(constructor->vcf_record.ref == "G");
        REQUIRE(constructor->vcf_record.alt == "GTT,T");

        REQUIRE(!constructor->read_tabix_record());
      }

      delete constructor;
    }
  }
}

TEST_CASE("CRAM read support", "[seqanhts]")
{
  SECTION("Sample file test")
  {
    SECTION("HtsSequenceRecord")
    {
      std::stringstream base_path;

      SECTION("SAM")
      {
        base_path << gyper_SOURCE_DIRECTORY << "/test/reference/test.sam";
      }

      SECTION("BAM")
      {
        base_path << gyper_SOURCE_DIRECTORY << "/test/reference/test.bam";
      }

      SECTION("CRAM")
      {
        base_path << gyper_SOURCE_DIRECTORY << "/test/reference/test.cram";
      }

      std::string base_path_str = base_path.str();
      seqan::HtsSequenceRecord hts_record;
      seqan::HtsFile hts_file(base_path_str.c_str(), "r");
      // seqan::open(hts_file);
      REQUIRE(seqan::readRecord(hts_record, hts_file));
      REQUIRE(hts_record.qName == "B7_591:4:96:693:509");
      seqan::String<seqan::Iupac> expected_sequence = "CACTAGTGGCTCATTGTAAATGTGTGGTTTAACTCG";
      REQUIRE(hts_record.seq == expected_sequence);

      for (unsigned i = 0; i < 3305; ++i)
      {
        seqan::readRecord(hts_record, hts_file);
      }

      REQUIRE(seqan::readRecord(hts_record, hts_file));
      REQUIRE(hts_record.qName == "EAS114_26:7:37:79:581");
      expected_sequence = "TTTTTTTTTTTTTTTTTTTTTTTCATGCCAGAAAA";
      REQUIRE(hts_record.seq == expected_sequence);
      REQUIRE(!seqan::readRecord(hts_record, hts_file));
    }
  }

  SECTION("HTS indexing")
  {
    std::stringstream base_path;

    SECTION("BAM")
    {
      base_path << gyper_SOURCE_DIRECTORY << "/test/reference/test.bam";
    }

    SECTION("CRAM")
    {
      base_path << gyper_SOURCE_DIRECTORY << "/test/reference/test.cram";
    }

    std::string base_path_str = base_path.str();
    seqan::HtsSequenceRecord hts_record;
    seqan::HtsFile hts_file(base_path_str.c_str(), "r");
    seqan::loadIndex(hts_file);
    seqan::setRegion(hts_file, "chr1:238-239");
    REQUIRE(seqan::readRegion(hts_record, hts_file));

    for (unsigned i = 0; i < 39; ++i)
    {
      seqan::readRegion(hts_record, hts_file);
    }

    REQUIRE(seqan::readRegion(hts_record, hts_file));
    REQUIRE(!seqan::readRegion(hts_record, hts_file));
  }
}

TEST_CASE("CRAM write support", "[seqanhts]")
{
  SECTION("Copy test cram file")
  {
    // No requirements but .../test/reference/test_write.cram should now have the same data as ../test/reference/test.cram
    std::stringstream read_path;
    read_path << gyper_SOURCE_DIRECTORY << "/test/reference/test.cram";
    std::string read_path_str = read_path.str();
    seqan::HtsFileIn hts_read_file(read_path_str.c_str());

    std::stringstream write_path;
    write_path << gyper_SOURCE_DIRECTORY << "/test/reference/test_write.cram";
    std::string write_path_str = write_path.str();
    seqan::HtsFileOut hts_write_file(write_path_str.c_str());

    seqan::copyHeader(hts_write_file, hts_read_file);
    seqan::writeHeader(hts_write_file);

    while (seqan::readRecord(hts_read_file))
    {
      seqan::copyRecord(hts_write_file, hts_read_file);
      // Here you could change the hts_write_file.hts_record if you want.
      seqan::writeRecord(hts_write_file);
    }
  }
}
