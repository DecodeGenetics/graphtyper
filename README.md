[![Build Status](https://travis-ci.org/DecodeGenetics/graphtyper.svg?branch=master)](https://travis-ci.org/DecodeGenetics/graphtyper)

# GraphTyper
GraphTyper is a highly scalable genotyping software. It represents a reference genome and known variants of a genomic region using an acyclic mathematical graph structure (a "pangenome reference"), which high-throughput sequence reads are re-aligned to for the purpose of discovering and genotyping SNPs, small indels, and structural variants.

Maintainer: Hannes Pétur Eggertsson (Hannes.Eggertsson@decode.is)

## Citation
### Small variant calling
Hannes P. Eggertsson, Hakon Jonsson, Snaedis Kristmundsdottir, Eirikur Hjartarson, Birte Kehr, Gisli Masson, Florian Zink, Kristjan E. Hjorleifsson, Aslaug Jonasdottir, Adalbjorg Jonasdottir, Ingileif Jonsdottir, Daniel F. Gudbjartsson, Pall Melsted, Kari Stefansson, Bjarni V. Halldorsson. Graphtyper enables population-scale genotyping using pangenome graphs. *Nature Genetics* **49**, 1654–1660 (2017). [doi:10.1038/ng.3964](http://dx.doi.org/10.1038/ng.3964)

### Strucural variant genotyping
Eggertsson, H.P., Kristmundsdottir, S., Beyter, D. et al. GraphTyper2 enables population-scale genotyping of structural variation using pangenome graphs. *Nature Communications* **10**, 5402 (2019) [doi:10.1038/s41467-019-13341-9](https://www.nature.com/articles/s41467-019-13341-9)

## Getting started
### Download binary
The easiest way to get GraphTyper is to download the latest 64bit linux binary release
```sh
wget https://github.com/DecodeGenetics/graphtyper/releases/download/v2.0/graphtyper
chmod a+x graphtyper
```

The binary is linked statically and therefore does not require any runtime libraries. Refer to our [recommended pipeline](https://github.com/DecodeGenetics/graphtyper-pipelines) for running the tool.

### Build from source
Alternatively you may want to build GraphTyper from source. First, you'll need the following:
* C++ compiler with C++11 supported (we mostly use gcc 4.8.5, gcc 5.4.0, and the latest gcc 6 and 7 releases)
* Boost>=1.57.0
* zlib>=1.2.8
* libbz2
* liblzma
* Autotools, Automake, libtool, Make, and CMake>=2.8.8 (if you want to use our build system)

All other dependencies are submodules of this repository. Make sure have the `CXX` environment variable set as the same compiler as `which g++` returns (because some of the submodules use the compiler directed by the `CXX` variable while other ignore it). Also set the `BOOST_ROOT` variable to the root of BOOST which should already be compiled with the same compiler. Graphtyper is linked with BOOST dynamically, but other libraries statically.

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

## License
GNU GPLv3
