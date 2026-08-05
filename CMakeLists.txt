cmake_minimum_required(VERSION 3.21)
project(OpenDIAlyzer VERSION 0.3.0 LANGUAGES CXX)

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

# In-process semi-supervised LDA (pyprophet-like) + its oracle test. Header-only,
# pure C++ (no OpenMS/Eigen), so the same header links into OpenDIAlyzer and the
# standalone test. Implementation delegated to codex; verified by the oracle.
add_executable(odia-lda-test src/odia_lda_test.cpp)
add_executable(odia-lda-adversarial-test src/odia_lda_adversarial_test.cpp)

# Error control above the precursor: peptide/protein context FDR, picked protein-group
# competition, entrapment FDP. Header-only and dependency-free, same as the LDA.
add_executable(odia-fdr-test src/odia_fdr_test.cpp)

# Compact XIC storage replacing vector<MSChromatogram>: shared RT grid + float32 intensities.
# Measured 4.09x smaller (225 GB -> 55 GB on the Astral benchmark). Header-only, standalone.
add_executable(odia-chromatogram-test src/odia_chromatogram_test.cpp)

# Histogram gradient-boosted trees as an alternative learner inside the same semi-supervised
# loop as the LDA. In-house because no GBT library is available to link and a vendored one would
# bring its own RNG and threading into a path that must stay deterministic (it feeds an FDR).
add_executable(odia-gbt-test src/odia_gbt_test.cpp)

# Interned/substring sequence storage, and the library object graph built on it. Header-only and
# dependency-free like the other kernels. These exist because a phase-resolved memory probe showed
# 83% of peak RSS was glibc arena fragmentation from ~471M string allocations, not any data
# structure -- a peptide is a substring of its protein and needs no characters of its own.
add_executable(odia-seqstore-test src/odia_seqstore_test.cpp)
add_executable(odia-library-test src/odia_library_test.cpp)
add_executable(odia-features-test src/odia_features_test.cpp)
add_executable(odia-prefilter-model-test src/odia_prefilter_model_test.cpp)
add_executable(odia-nn-test src/odia_nn_test.cpp)
add_executable(odia-anchor-training-test src/odia_anchor_training_test.cpp)
add_executable(odia-ablate src/odia_ablate.cpp)
add_executable(odia-nn-diag src/odia_nn_diag.cpp)
add_executable(odia-chromstore-test src/odia_chromstore_test.cpp)
add_executable(odia-rtaxis-test src/odia_rtaxis_test.cpp)

# The semi-supervised loop on a REAL score matrix (44,568 peak groups x 29 sub-scores extracted
# from the Astral benchmark), not synthetic data. Guarded on the fixture: the test is skipped
# rather than failed when testdata/lda_fixture.txt is absent, so a clone without it still passes.
add_executable(odia-lda-realdata-test src/odia_lda_realdata_test.cpp)
# Two test files existed in src/ and were never added here, so they had never been compiled, let
# alone run. Both cover code that is LIVE in the tool: odia_split.h drives the extraction split and
# odia_mzpeak_access.h is the mzPeak input layer. "19/19 tests pass" was true and misleading.
add_executable(odia-split-test src/odia_split_test.cpp)
add_executable(odia-spectrumstore-test src/odia_spectrumstore_test.cpp)
add_executable(odia-scored-test src/odia_scored_test.cpp)

# The kernels in odia_lda.h / odia_gbt.h contain `#pragma omp` regions. Without linking OpenMP,
# _OPENMP is undefined, every pragma is ignored, and the tests that set n_threads = 1/8/64 run
# THREE IDENTICAL SERIAL FITS -- reporting "bit-identical across thread counts" while never having
# started a second thread. Link it, so the determinism claim is actually exercised.
if(OpenMP_CXX_FOUND)
  target_link_libraries(odia-lda-test            PRIVATE OpenMP::OpenMP_CXX)
  target_link_libraries(odia-lda-adversarial-test PRIVATE OpenMP::OpenMP_CXX)
  target_link_libraries(odia-lda-realdata-test   PRIVATE OpenMP::OpenMP_CXX)
  target_link_libraries(odia-gbt-test            PRIVATE OpenMP::OpenMP_CXX)
  # The NN test asserts BIT-IDENTICAL results across 1/2/7/16/64 threads. Unlinked, _OPENMP is
  # undefined, every pragma vanishes, and it runs five identical serial fits -- the precise failure
  # this block was written to prevent, reproduced one target lower. odia-ablate needs it to finish
  # this century.
  target_link_libraries(odia-nn-test             PRIVATE OpenMP::OpenMP_CXX)
  target_link_libraries(odia-anchor-training-test PRIVATE OpenMP::OpenMP_CXX)
  target_link_libraries(odia-ablate              PRIVATE OpenMP::OpenMP_CXX)
  target_link_libraries(odia-nn-diag             PRIVATE OpenMP::OpenMP_CXX)
endif()

# The GBT's determinism guarantee is "the threaded result equals the serial one". That is only a
# real assertion if a genuinely serial build exists to compare against, so build one explicitly.
add_executable(odia-gbt-test-serial src/odia_gbt_test.cpp)

# Stage-1 library generator: peptide list -> PeptDeep-predicted RT/MS2 ->
# OpenSWATH transition TSV. Requires the ONNX-enabled OpenMS build
# (-DOpenMS_DIR=.../openms-onnx/lib/cmake/OpenMS). The PeptDeep ONNX models
# are not installed; OpenMS's CMake downloads them into its build tree under
# share/OpenMS/models, so the selftest needs that directory. The default below
# assumes OpenMS_DIR points into that build tree; override with
# -DOPENDIALIBGEN_MODEL_DIR=... if the layout differs.
add_executable(OpenDIALibGen src/opendialibgen.cpp)
target_link_libraries(OpenDIALibGen PRIVATE OpenMS)

# OpenDIAlyzer — the streaming DIA search tool (skeleton). TOPPBase-derived,
# links libOpenMS. It will additionally link libmzpeak for streaming input once
# the reader stage lands (target_link_libraries(... mzpeak) + include dir).
# Companion to OpenDIALibGen; the two form the OpenDIAlyzer tool suite that USES
# OpenMS but does not live in the OpenMS repository.
add_executable(OpenDIAlyzer src/opendialyzer.cpp)
target_link_libraries(OpenDIAlyzer PRIVATE OpenMS)

# mzPeak streaming input. Optional: point -DMZPEAK_ROOT at an mzpeak-cpp build (headers in
# <root>/include, libmzpeak in <root>/<builddir>). Without it OpenDIAlyzer still builds and
# simply does not accept .mzpeak.
if(MZPEAK_ROOT)
  find_library(MZPEAK_LIB NAMES mzpeak PATHS ${MZPEAK_ROOT} PATH_SUFFIXES buildfix build build-fix REQUIRED)
  find_package(Arrow REQUIRED)
  find_package(Parquet REQUIRED)
  target_include_directories(OpenDIAlyzer PRIVATE ${MZPEAK_ROOT}/include ${CMAKE_CURRENT_SOURCE_DIR}/src)
  target_link_libraries(OpenDIAlyzer PRIVATE ${MZPEAK_LIB} Arrow::arrow_shared Parquet::parquet_shared)
  target_compile_definitions(OpenDIAlyzer PRIVATE WITH_MZPEAK=1)
  add_executable(odia-mzpeak-access-test src/odia_mzpeak_access_test.cpp)
  target_include_directories(odia-mzpeak-access-test PRIVATE ${MZPEAK_ROOT}/include ${CMAKE_CURRENT_SOURCE_DIR}/src)
  target_link_libraries(odia-mzpeak-access-test PRIVATE ${MZPEAK_LIB} Arrow::arrow_shared Parquet::parquet_shared OpenMS)
  target_compile_definitions(odia-mzpeak-access-test PRIVATE WITH_MZPEAK=1)
  # Needs a real .mzpeak; exits 2 with a usage message when given none. Point it at a fixture if the
  # environment has one and let ctest SKIP otherwise, rather than fail on a missing file -- the same
  # arrangement odia-lda-realdata uses. ODIA_MZPEAK_FIXTURE overrides for a different environment.
  if(NOT DEFINED ODIA_MZPEAK_FIXTURE)
    set(ODIA_MZPEAK_FIXTURE "$ENV{ODIA_MZPEAK_FIXTURE}")
  endif()
  if(ODIA_MZPEAK_FIXTURE AND EXISTS "${ODIA_MZPEAK_FIXTURE}")
    add_test(NAME odia-mzpeak-access COMMAND odia-mzpeak-access-test "${ODIA_MZPEAK_FIXTURE}")
    set_tests_properties(odia-mzpeak-access PROPERTIES SKIP_RETURN_CODE 77)
    message(STATUS "OpenDIAlyzer: mzPeak access test fixture ${ODIA_MZPEAK_FIXTURE}")
  else()
    message(STATUS "OpenDIAlyzer: mzPeak access test NOT registered -- set ODIA_MZPEAK_FIXTURE to a .mzpeak")
  endif()
  message(STATUS "OpenDIAlyzer: mzPeak streaming input ENABLED (${MZPEAK_LIB})")
endif()
if(OpenMP_CXX_FOUND)
  target_link_libraries(OpenDIAlyzer PRIVATE OpenMP::OpenMP_CXX)
endif()
# OpenDIAlyzer reads the pass-1 .osw with sqlite3 to harvest RT anchors. libOpenMS
# already exports the sqlite3_* symbols (it vendors sqlite3), so we need only the
# header on the include path -- no extra link. The header ships with OpenMS (its
# installed include/, or the vendored extern/ tree). ponytail: header-only hookup.
find_path(SQLITE3_INCLUDE_DIR sqlite3.h
  HINTS
    "${OpenMS_DIR}/_deps/opentims-src/src/opentims++/sqlite"
    "${OpenMS_DIR}/../../include"
    "${OpenMS_DIR}/../../ext/OpenMS/src/openms/extern/SQLiteCpp/sqlite3"
  DOC "Directory containing sqlite3.h (ships with OpenMS)")
if(NOT SQLITE3_INCLUDE_DIR)
  message(FATAL_ERROR "sqlite3.h not found; pass -DSQLITE3_INCLUDE_DIR=/path/to/dir")
endif()
target_include_directories(OpenDIAlyzer PRIVATE ${SQLITE3_INCLUDE_DIR})

# P2 mod_x featurization validation. Links libOpenMS for AASequence/EmpiricalFormula;
# the featurization itself is a header-only vendored patch (PeptDeepModX.h), so
# this needs no OpenMS rebuild. Validated against AlphaPeptDeep's reference fixture.
add_executable(odia-modx-test src/odia_modx_test.cpp)
target_link_libraries(odia-modx-test PRIVATE OpenMS)

add_executable(odia-modpredict-test src/odia_modpredict_test.cpp)
target_link_libraries(odia-modpredict-test PRIVATE OpenMS)

set(OPENDIALIBGEN_MODEL_DIR "${OpenMS_DIR}/../../../share/OpenMS/models"
    CACHE PATH "Directory containing peptdeep_*_dynamic.onnx (used by the OpenDIALibGen selftest)")

enable_testing()
add_test(NAME odia-info-selftest COMMAND odia-info --selftest)
add_test(NAME odia-index-selftest COMMAND odia-index --selftest)
add_test(NAME odia-score-oracle COMMAND odia-score-test)
add_test(NAME odia-lda COMMAND odia-lda-test)
add_test(NAME odia-lda-adversarial COMMAND odia-lda-adversarial-test)
add_test(NAME odia-fdr COMMAND odia-fdr-test)
add_test(NAME odia-chromatogram COMMAND odia-chromatogram-test)
add_test(NAME odia-gbt COMMAND odia-gbt-test)
add_test(NAME odia-seqstore COMMAND odia-seqstore-test)
add_test(NAME odia-library COMMAND odia-library-test)
add_test(NAME odia-gbt-serial COMMAND odia-gbt-test-serial)
# Built but never RUN by ctest until now, so a broken minibatch implementation left the suite green.
add_test(NAME odia-nn COMMAND odia-nn-test)
add_test(NAME odia-anchor-training COMMAND odia-anchor-training-test)
add_test(NAME odia-chromstore COMMAND odia-chromstore-test)
add_test(NAME odia-rtaxis COMMAND odia-rtaxis-test)
# BUILT BUT NEVER RUN until now -- the same trap the odia-nn comment above records. These two are
# the designed fixes for the two largest measured problems (feature memory, prefilter recall), so a
# silent regression in either would have gone unnoticed precisely while they were being wired in.
add_test(NAME odia-features COMMAND odia-features-test)
add_test(NAME odia-prefilter-model COMMAND odia-prefilter-model-test)
add_test(NAME odia-split COMMAND odia-split-test)
add_test(NAME odia-spectrumstore COMMAND odia-spectrumstore-test)
add_test(NAME odia-scored COMMAND odia-scored-test)
add_test(NAME odia-lda-realdata
         COMMAND odia-lda-realdata-test "${CMAKE_CURRENT_SOURCE_DIR}/testdata/lda_fixture.txt")
set_tests_properties(odia-lda-realdata PROPERTIES SKIP_RETURN_CODE 77)
add_test(NAME OpenDIALibGen-selftest
         COMMAND OpenDIALibGen --selftest -model_dir "${OPENDIALIBGEN_MODEL_DIR}")
add_test(NAME odia-modx-p2 COMMAND odia-modx-test)
add_test(NAME odia-modpredict-p2 COMMAND odia-modpredict-test "${OPENDIALIBGEN_MODEL_DIR}")
