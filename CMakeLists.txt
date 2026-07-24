cmake_minimum_required(VERSION 3.21)
project(OpenDIAlyzer VERSION 0.0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# OpenMS is a prerequisite, not a host: we link against an installation and
# live outside its source tree. Point at it with
#   -DOpenMS_DIR=/path/to/openms/lib/cmake/OpenMS
# An OpenMS built against its own (static) contrib exports link-interface
# targets like CURL::libcurl_static. A consumer resolving dependencies from
# elsewhere (system, conda) only has the shared CURL::libcurl, and CMake aborts
# on the missing target even though the code is already inside libOpenMS.
# ponytail: alias the name so the target exists — the symbols are statically
# linked into libOpenMS already, so nothing extra actually needs linking.
find_package(CURL REQUIRED)
if(TARGET CURL::libcurl AND NOT TARGET CURL::libcurl_static)
  add_library(CURL::libcurl_static ALIAS CURL::libcurl)
endif()

find_package(OpenMS REQUIRED)
find_package(OpenMP)

add_executable(odia-info src/odia_info.cpp)
target_link_libraries(odia-info PRIVATE OpenMS)
if(OpenMP_CXX_FOUND)
  target_link_libraries(odia-info PRIVATE OpenMP::OpenMP_CXX)
endif()

# Streaming indexer. Uses std::thread directly rather than OpenMP: the work is
# partitioned by isolation window so no two threads share a cell, and OpenMP's
# default spin-waiting measurably hurt on this hardware (224 threads were 2x
# slower than 1 for mzML ingest).
find_package(Threads REQUIRED)
add_executable(odia-index src/odia_index.cpp)
target_link_libraries(odia-index PRIVATE OpenMS Threads::Threads)

# Scoring kernel (header-only) + its oracle test. No OpenMS link: the kernel is
# self-contained C++ so the same source can later compile for CUDA. Correctness
# is defined by OpenSWATH's reference values, asserted in the test.
add_executable(odia-score-test src/odia_score_test.cpp)

# Stage-1 library generator: peptide list -> PeptDeep-predicted RT/MS2 ->
# OpenSWATH transition TSV. Requires the ONNX-enabled OpenMS build
# (-DOpenMS_DIR=.../openms-onnx/lib/cmake/OpenMS). The PeptDeep ONNX models
# are not installed; OpenMS's CMake downloads them into its build tree under
# share/OpenMS/models, so the selftest needs that directory. The default below
# assumes OpenMS_DIR points into that build tree; override with
# -DOPENDIALIBGEN_MODEL_DIR=... if the layout differs.
add_executable(OpenDIALibGen src/opendialibgen.cpp)
target_link_libraries(OpenDIALibGen PRIVATE OpenMS)

# P2 mod_x featurization validation. Links libOpenMS for AASequence/EmpiricalFormula;
# the featurization itself is a header-only vendored patch (PeptDeepModX.h), so
# this needs no OpenMS rebuild. Validated against AlphaPeptDeep's reference fixture.
add_executable(odia-modx-test src/odia_modx_test.cpp)
target_link_libraries(odia-modx-test PRIVATE OpenMS)

set(OPENDIALIBGEN_MODEL_DIR "${OpenMS_DIR}/../../../share/OpenMS/models"
    CACHE PATH "Directory containing peptdeep_*_dynamic.onnx (used by the OpenDIALibGen selftest)")

enable_testing()
add_test(NAME odia-info-selftest COMMAND odia-info --selftest)
add_test(NAME odia-index-selftest COMMAND odia-index --selftest)
add_test(NAME odia-score-oracle COMMAND odia-score-test)
add_test(NAME OpenDIALibGen-selftest
         COMMAND OpenDIALibGen --selftest -model_dir "${OPENDIALIBGEN_MODEL_DIR}")
add_test(NAME odia-modx-p2 COMMAND odia-modx-test)
