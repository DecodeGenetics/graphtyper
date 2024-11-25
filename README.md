[![Build Status](https://github.com/DecodeGenetics/graphtyper/actions/workflows/ci_linux.yaml/badge.svg?branch=master)](https://github.com/DecodeGenetics/graphtyper/actions/workflows/ci_linux.yaml?query=branch%3Amaster) [![Format Status](https://github.com/DecodeGenetics/graphtyper/actions/workflows/clang-format.yaml/badge.svg?branch=master)](https://github.com/DecodeGenetics/graphtyper/actions/workflows/clang-format.yaml?query=branch%3Amaster) [![Conda](https://img.shields.io/conda/pn/bioconda/graphtyper?color=green)](http://bioconda.github.io/recipes/graphtyper/README.html) [![Conda](https://img.shields.io/conda/v/bioconda/graphtyper?color=green)](http://bioconda.github.io/recipes/graphtyper/README.html)

# graphtyper
graphtyper is a graph-based variant caller capable of genotyping population-scale short read data sets. It represents a reference genome and known variants of a genomic region using an acyclic graph structure (a "pangenome reference"), which high-throughput sequence reads are re-aligned to for the purpose of discovering and genotyping SNPs, small indels, and structural variants.

Maintainer: Hannes Pétur Eggertsson (Hannes.Eggertsson@decode.is)

## Installation
### Static binary release
The easiest way to install GraphTyper is go to "Releases" and download the latest binary, here: https://github.com/DecodeGenetics/graphtyper/releases

The binary is linked statically and therefore does not require any runtime libraries. If you prefer, you can also install graphtyper via bioconda: http://bioconda.github.io/recipes/graphtyper/README.html

### Building from source
You may also want to build graphtyper from source, for example if you want to make changes to the code. In this case, you'll first need the following:
* C++ compiler with full AVX512 support (GCC 8+)
* Boost>=1.57.0
* zlib>=1.2.8
* libbz2
* liblzma
* Autotools, Automake, libtool, Make, and CMake>=3.2 (if you want to use our build system)

All other dependencies are submodules of this repository. Make sure have the `CXX` environment variable set as the same compiler as `which g++` returns (because some of the submodules use the compiler directed by the `CXX` variable while other ignore it). Also set the `BOOST_ROOT` variable to the root of BOOST which should already be compiled with the same compiler. Graphtyper is linked with BOOST dynamically, but other libraries statically.

For the purpose of demonstration, we assume you want to clone graphtyper to `~/git/graphtyper` and build it in `~/git/graphtyper/release-build`.

```sh
mkdir -p ~/git && cd ~/git
git clone --recursive https://github.com/DecodeGenetics/graphtyper.git graphtyper && cd graphtyper
mkdir -p release-build && cd release-build
cmake ..
make -j4 graphtyper # The 'j' argument specifies how many compilation threads to use, you can change this if you have more threads available. Also, the compilation will take awhile... consider getting coffee at this point.
bin/graphtyper # Will run graphtyper for the very first time!
```
And that's all. If you are lucky enough to have administrative access, you can run `sudo make install` to install graphtyper system-wide.


## Usage

The recommended way of genotyping small variants (SNP+indels) is using the `genotype` subcommand

```sh
./graphtyper genotype <REFERENCE.fa> --sams=<BAMLIST_OR_CRAMLIST> --region=<chrA:begin-end> --threads=<T>
```

and use the `genotype_sv` subcommand for genotyping structural variants

```sh
./graphtyper genotype_sv <REFERENCE.fa> <input.vcf.gz> --sams=<BAMLIST_OR_CRAMLIST> --region=<chrA:begin-end> --threads=<T>
```

See the [graphtyper user guide](https://github.com/DecodeGenetics/graphtyper/wiki/User-guide) for more details.


## Citation
### Small variant genotyping
Hannes P. Eggertsson, Hakon Jonsson, Snaedis Kristmundsdottir, Eirikur Hjartarson, Birte Kehr, Gisli Masson, Florian Zink, Kristjan E. Hjorleifsson, Aslaug Jonasdottir, Adalbjorg Jonasdottir, Ingileif Jonsdottir, Daniel F. Gudbjartsson, Pall Melsted, Kari Stefansson, Bjarni V. Halldorsson. Graphtyper enables population-scale genotyping using pangenome graphs. *Nature Genetics* **49**, 1654–1660 (2017). [doi:10.1038/ng.3964](http://dx.doi.org/10.1038/ng.3964)

### Structural variant genotyping
Eggertsson, H.P., Kristmundsdottir, S., Beyter, D. et al. GraphTyper2 enables population-scale genotyping of structural variation using pangenome graphs. *Nature Communications* **10**, 5402 (2019) [doi:10.1038/s41467-019-13341-9](https://www.nature.com/articles/s41467-019-13341-9)


## License
MIT License
