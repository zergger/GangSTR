set(GANGSTR_SOURCES
    bam_info_extract.cpp
    bam_io.cpp
    common.cpp
    enclosing_class.cpp
    flanking_class.cpp
    frr_class.cpp
    gc_region_reader.cpp
    genotyper.cpp
    likelihood_maximizer.cpp
    locus.cpp
    main_gangstr.cpp
    mathops.cpp
    options.cpp
    read_class.cpp
    read_extractor.cpp
    read_pair.cpp
    realignment.cpp
    ref_genome.cpp
    region_reader.cpp
    sample_info.cpp
    spanning_class.cpp
    ssw.c
    ssw_cpp.cpp
    str_info.cpp
    stringops.cpp
    vcf_writer.cpp
    hipstr_models/hip_em_learner.cpp
    hipstr_models/hip_stutter_model.cpp
)

add_library(gstr ${GANGSTR_SOURCES})
target_link_libraries(gstr ${htslib_static} ${zlib_static} ${nlopt_static})
target_include_directories(gstr PUBLIC "${PROJECT_BINARY_DIR}")
