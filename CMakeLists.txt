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

enable_testing()
add_test(NAME odia-info-runs COMMAND odia-info --selftest)
