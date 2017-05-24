cmake_minimum_required(VERSION 2.8.8)

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
  graph/reference_depth.cpp
  graph/ref_node.cpp
  graph/sequence_extractor.cpp
  graph/var_node.cpp
  graph/var_record.cpp
  index/indexer.cpp
  index/mem_index.cpp
  index/rocksdb.cpp
  typer/aligner.cpp
  typer/alignment.cpp
  typer/caller.cpp
  typer/discovery.cpp
  typer/genotype_paths.cpp
  typer/graph_swapper.cpp
  typer/path.cpp
  typer/sample_call.cpp
  typer/segment.cpp
  typer/segment_calling.cpp
  typer/var_stats.cpp
  typer/variant.cpp
  typer/variant_candidate.cpp
  typer/variant_map.cpp
  typer/variant_support.cpp
  typer/vcf.cpp
  typer/vcf_operations.cpp
  typer/vcf_writer.cpp
  utilities/adapter_removal.cpp
  utilities/io.cpp
  utilities/kmer_help_functions.cpp
  utilities/type_conversions.cpp
  utilities/sam_reader.cpp
  utilities/options.cpp
)

# Object libarary
add_library(graphtyper_objects OBJECT ${graphtyper_SOURCE_FILES})
add_dependencies(graphtyper_objects project_htslib) # Depends on building htslib
add_dependencies(graphtyper_objects project_sparsehash) # Depends on building sparsehash
add_dependencies(graphtyper_objects project_statgen) # Depends on building statgen
add_dependencies(graphtyper_objects project_snappy) # Depends on building snappy
# add_dependencies(graphtyper_objects project_rocksdb) # Depends on building rocksdb

# Static library
add_library(static STATIC $<TARGET_OBJECTS:graphtyper_objects>)
SET_TARGET_PROPERTIES(static PROPERTIES OUTPUT_NAME graphtyper)

# Binary
add_executable(graphtyper main.cpp $<TARGET_OBJECTS:graphtyper_objects>)
target_link_libraries(graphtyper ${graphtyper_all_libraries})

# Installs graphtyper to /usr/local/bin using [sudo] make install
install (TARGETS graphtyper DESTINATION bin)
