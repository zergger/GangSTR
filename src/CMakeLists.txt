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

if(BUILD_TESTING)
  find_package(PkgConfig QUIET)
  if(PkgConfig_FOUND)
    pkg_check_modules(CPPUNIT QUIET cppunit)
  endif()
  find_path(CPPUNIT_INCLUDE_DIR cppunit/TestCase.h
      /usr/include
      /usr/local/include
  )
  find_library(CPPUNIT_LIBRARY cppunit
      /usr/lib
      /usr/lib/x86_64-linux-gnu
      /usr/local/lib
  )

  add_executable(gangstr_integration_smoke
      tests/gangstr_integration_smoke.cpp
  )
  target_link_directories(gangstr_integration_smoke PRIVATE
      ${CMAKE_BINARY_DIR}/src
      ${CMAKE_BINARY_DIR}/thirdparty/zlib/lib
      ${CMAKE_BINARY_DIR}/thirdparty/gsl/lib
      ${CMAKE_BINARY_DIR}/thirdparty/nlopt/lib64
      ${CMAKE_BINARY_DIR}/thirdparty/nlopt/lib
      ${CMAKE_BINARY_DIR}/thirdparty/htslib/lib
  )
  target_link_libraries(gangstr_integration_smoke PRIVATE gstr Threads::Threads hts z m nlopt gsl)
  target_include_directories(gangstr_integration_smoke PRIVATE "${PROJECT_BINARY_DIR}")
  add_test(NAME gangstr_integration_smoke COMMAND gangstr_integration_smoke)

  add_executable(hip_em_learner_unit
      tests/hip_em_learner_unit.cpp
  )
  target_link_directories(hip_em_learner_unit PRIVATE
      ${CMAKE_BINARY_DIR}/src
      ${CMAKE_BINARY_DIR}/thirdparty/zlib/lib
      ${CMAKE_BINARY_DIR}/thirdparty/gsl/lib
      ${CMAKE_BINARY_DIR}/thirdparty/nlopt/lib64
      ${CMAKE_BINARY_DIR}/thirdparty/nlopt/lib
      ${CMAKE_BINARY_DIR}/thirdparty/htslib/lib
  )
  target_link_libraries(hip_em_learner_unit PRIVATE gstr Threads::Threads hts z m nlopt gsl)
  target_include_directories(hip_em_learner_unit PRIVATE "${PROJECT_BINARY_DIR}")
  add_test(NAME hip_em_learner_unit COMMAND hip_em_learner_unit)

  if(CPPUNIT_FOUND OR (CPPUNIT_INCLUDE_DIR AND CPPUNIT_LIBRARY))
    if(CMAKE_CXX_COMPILER MATCHES "conda")
      message(WARNING "Skipping legacy CPPUNIT suite under the Conda toolchain; build with CC=/usr/bin/gcc CXX=/usr/bin/g++ to run gangstr_cppunit.")
    else()
      add_executable(gangstr_cppunit
          tests/main_test.cpp
          tests/ReadClass_test.cpp
          tests/LikelihoodMaximizer_test.cpp
          tests/ReadExtractor_test.cpp
      )
      target_compile_definitions(gangstr_cppunit PRIVATE
          GANGSTR_TEST_DATA_DIR="${PROJECT_SOURCE_DIR}/unused/tests"
      )
      target_link_directories(gangstr_cppunit PRIVATE
          ${CMAKE_BINARY_DIR}/src
          ${CMAKE_BINARY_DIR}/thirdparty/zlib/lib
          ${CMAKE_BINARY_DIR}/thirdparty/gsl/lib
          ${CMAKE_BINARY_DIR}/thirdparty/nlopt/lib64
          ${CMAKE_BINARY_DIR}/thirdparty/nlopt/lib
          ${CMAKE_BINARY_DIR}/thirdparty/htslib/lib
      )
      target_include_directories(gangstr_cppunit PRIVATE
          "${PROJECT_BINARY_DIR}"
          ${CPPUNIT_INCLUDE_DIRS}
          ${CPPUNIT_INCLUDE_DIR}
      )
      target_link_libraries(gangstr_cppunit PRIVATE
          gstr
          Threads::Threads
          ${CPPUNIT_LIBRARIES}
          ${CPPUNIT_LIBRARY}
          hts
          z
          m
          nlopt
          gsl
      )
      target_link_directories(gangstr_cppunit PRIVATE
          ${CPPUNIT_LIBRARY_DIRS}
          /usr/lib
          /usr/lib/x86_64-linux-gnu
          /usr/local/lib
      )
      add_test(NAME gangstr_cppunit COMMAND gangstr_cppunit)
    endif()
  endif()
endif()
