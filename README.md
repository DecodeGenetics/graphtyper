# Graphtyper
Graphtyper is a highly scalable genotyping software. It represents a reference genome and known variants of a genomic region using an acyclic mathematical graph structure (a "pangenome reference"), which high-throughput sequence reads are re-aligned to for the purpose of discovering and genotyping SNPs and small indels.

First author and software maintainer: Hannes Pétur Eggertsson (Hannes.Eggertsson@decode.is)

## Citation
Hannes P. Eggertsson, Hakon Jonsson, Snaedis Kristmundsdottir, Eirikur Hjartarson, Birte Kehr, Gisli Masson, Florian Zink, Kristjan E. Hjorleifsson, Aslaug Jonasdottir, Adalbjorg Jonasdottir, Ingileif Jonsdottir, Daniel F. Gudbjartsson, Pall Melsted, Kari Stefansson, Bjarni V. Halldorsson. Graphtyper enables population-scale genotyping using pangenome graphs. *Nature Genetics*. [doi:10.1038/ng.3964](http://dx.doi.org/10.1038/ng.3964) (2017)

## Getting started
### Dependencies
Before installing Graphtyper you'll need the following:
* C++ compiler with C++11 supported (we mostly use gcc 4.8.5, gcc 5.4.0, and the latest gcc 6 and 7 releases)
* Boost>=1.57.0
* zlib>=1.2.8
* libbz2
* liblzma
* Autotools, Automake, libtool, Make, and CMake>=2.8.8 (if you want to use our build system)

All other dependencies are submodules of this repository.

### Installation
Make sure have the `CXX` environment variable set as the same compiler as `which g++` returns (because some of the submodules use the compiler directed by the `CXX` variable while other ignore it). Also set the `BOOST_ROOT` variable to the root of BOOST which should already be compiled with the same compiler. Graphtyper is linked with BOOST dynamically, but other libraries statically.

For the purpose of demonstration, we assume you want to clone Graphtyper to `~/git/graphtyper` and build it in `~/git/graphtyper/release-build`.

```sh
mkdir -p ~/git && cd ~/git
git clone --recursive https://github.com/DecodeGenetics/graphtyper.git graphtyper && cd graphtyper
mkdir -p release-build && cd release-build
cmake ..
make -j4 graphtyper # The 'j' argument specifies how many compilation threads to use, you can change this if you have more threads available. Also, the compilation will take awhile... consider getting coffee at this point.
bin/graphtyper # Will run Graphtyper for the very first time!
```
And that's all. If you are lucky enough to have administrative access, you can run `sudo make install` to install Graphtyper system-wide.

### Testing Graphtyper
To make sure you have a good release and everything compiled as expected, you can run Graphtyper's unit tests with:
```sh
make -j4 # Will compile everything, including tests
make test # Will run the unit tests
```

## Usage
Please use the [recommended pipeline](https://github.com/DecodeGenetics/graphtyper-pipelines).

## Variant filtering
Graphtyper is designed to have high specificity. That also means it often has many false positive sequence variant calls, in particular when genotyping large datasets. If you want to reduce false positives, you can post-process the VCF files using a tool called `vcffilter` which is a part of [vcflib](https://github.com/vcflib/vcflib).

In our manuscript we used the following commands to filter:
* Single sample filtering: `vcffilter -f "ABHet < 0.0 | ABHet > 0.30" -f "MQ > 30" -f "QD > 6.0" <VCF>`
* Multi sample filtering: `vcffilter -f "ABHet < 0.0 | ABHet > 0.33" -f "ABHom < 0.0 | ABHom > 0.97" -f "MaxAASR > 0.4" -f "MQ > 30" <VCF>`

## Graphtyper commands (advanced)
All Graphtyper commands are in the following format:
```sh
graphtyper <COMMAND> [OPTIONS]
```
Run `graphtyper` to get a list of commands and running `graphtyper <COMMAND> --help` to get all available options for that command.

### Graph construction
Assuming you have variants in file 'variants.vcf' (not required) and your reference genome in 'genome.fa'. then you can create a graph of the region chr22:35,000,000-35,005,000 using
```sh
samtools faidx index genome.fa # Only required if the FASTA file is not already indexed
bgzip -c variants.vcf > variants.vcf.gz && tabix variants.vcf.gz # Only required if the VCF file is not already BGZF compressed and tabix indexed
graphtyper construct my_new_graph genome.fa --vcf variants.vcf.gz chr22:35,000,000-35,005,000
```

### Index construction
To index the graph use
```sh
graphtyper index my_new_graph
```
By default the index directory will be 'my_new_graph_gti'. Here, 'gti' stands for GraphTyper Index.

### Genotype call variants
Assuming you have SAM/BAM/CRAM files listed in a file called 'my_samples', you can call the samples by using only reads which are aligned to the graph's region using

```sh
graphtyper call my_new_graph --sams my_samples chr22:35,000,000-35,005,000
```

This will output the following files in the current working directory (if SAMPLE1 is your first sample):
 * SAMPLE1_calls.vcf.gz (only if you provided a VCF) - Genotype calls of all the samples with no overlapping variants. This VCF can be broken down into smaller variant calls.
 * SAMPLE1_variants.vcf.gz - A list of new variant candidates.
 * SAMPLE1.haps (only if you provided a VCF) - A BOOST serialization format of haplotype calls made.

All the files are prefixed with the first sample they have, to distinguish calls of pooled sample sets.

### Haplotype extraction
Extracting haplotypes is useful to limit the complexity of the variation in the graph

```sh
ls *.haps > "all.haps"
graphtyper haplotypes my_new_graph --haplotypes all.haps --output all_called_haplotypes.vcf.gz
```

### Merging Graphtyper VCF files
Graphtyper has implemented some commonly used VCF operations (note that they are only designed to work for Graphtyper VCFs!). These operations can be faster than similar operations using other tools because it makes assmptions which are true for Graphtyper VCF files, but may not be true for any VCF file.

#### VCF concatenate
Concatenate is useful if you want to concatenate haplotypes and new variants to generate a single VCF with all unique variants.

```sh
graphtyper vcf_concatenate all_called_haplotypes.vcf.gz *_variants.vcf.gz | bgzip -c > new_variants.vcf.gz && tabix new_variants.vcf.gz
```

#### VCF merge
Merging is useful if you call samples in pools and want to merge the sample genotype calls into a single VCF file.

```sh
graphtyper vcf_merge *_calls.vcf.gz | bgzip -c > calls_merged.vcf.gz && tabix calls_merged.vcf.gz
```

#### VCF break down
Breaking down VCF files is useful if you want to split complex variants into smaller sets of variants

```sh
graphtyper vcf_break_down my_new_graph calls_merged.vcf.gz | bgzip -c > calls_merged_broken_down.vcf.gz && tabix calls_merged_broken_down.vcf.gz
```

## License
GNU GPLv3
