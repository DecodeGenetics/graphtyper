cmake_minimum_required(VERSION 3.2)

set(graphtyper_SOURCE_FILES
  graph/absolute_position.cpp
  graph/constructor.cpp
  graph/genomic_region.cpp
  graph/genotype.cpp
  graph/graph.cpp
  graph/graph_serialization.cpp
  graph/haplotype.cpp
  graph/haplotype_calls.cpp
  graph/haplotype_extractor.cpp
  graph/label.cpp
  graph/read_strand.cpp
  graph/reference_depth.cpp
  graph/ref_node.cpp
  graph/sequence_extractor.cpp
  graph/sv.cpp
  graph/var_node.cpp
  graph/var_record.cpp
  index/indexer.cpp
  index/kmer_label.cpp
  index/ph_index.cpp
  typer/alignment.cpp
  typer/caller.cpp
  typer/genotype_paths.cpp
  typer/path.cpp
  typer/primers.cpp
  typer/sample_call.cpp
  typer/segment.cpp
  typer/var_stats.cpp
  typer/variant.cpp
  typer/variant_candidate.cpp
  typer/variant_map.cpp
  typer/variant_support.cpp
  typer/vcf.cpp
  typer/vcf_operations.cpp
  typer/vcf_writer.cpp
  utilities/bamshrink.cpp
  utilities/genotype.cpp
  utilities/genotype_camou.cpp
  utilities/genotype_sv.cpp
  utilities/hts_parallel_reader.cpp
  utilities/hts_reader.cpp
  utilities/hts_writer.cpp
  utilities/io.cpp
  utilities/kmer_help_functions.cpp
  utilities/options.cpp
  utilities/type_conversions.cpp
  utilities/sam_reader.cpp
  utilities/system.cpp
)

# Object libarary
add_library(graphtyper_objects OBJECT ${graphtyper_SOURCE_FILES})
add_dependencies(graphtyper_objects project_htslib) # Depends on building htslib

# Binary
add_executable(graphtyper main.cpp $<TARGET_OBJECTS:graphtyper_objects>)
target_link_libraries(graphtyper ${graphtyper_all_libraries})

# Installs graphtyper to /usr/local/bin using [sudo] make install
install (TARGETS graphtyper DESTINATION bin)
