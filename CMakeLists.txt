cmake_minimum_required(VERSION 3.21)
project(OpenDIAlyzer VERSION 0.0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# OpenMS is a prerequisite, not a host: we link against an installation and
# live outside its source tree. Point at it with
#   -DOpenMS_DIR=/path/to/openms/lib/cmake/OpenMS
find_package(OpenMS REQUIRED)
find_package(OpenMP)

add_executable(odia-info src/odia_info.cpp)
target_link_libraries(odia-info PRIVATE OpenMS)
if(OpenMP_CXX_FOUND)
  target_link_libraries(odia-info PRIVATE OpenMP::OpenMP_CXX)
endif()

enable_testing()
add_test(NAME odia-info-runs COMMAND odia-info --selftest)
