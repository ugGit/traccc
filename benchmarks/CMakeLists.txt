cmake_minimum_required(VERSION 3.14)

include(FetchContent)

# Project include(s).
include( traccc-compiler-options-cpp )

# Set up a common library, shared by all of the tests.
add_library( traccc_benchmarks_common INTERFACE )
target_include_directories( traccc_benchmarks_common
    INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/common )

# Add cpu benchmarks
add_subdirectory(cpu)

# Add CUDA benchmarks
if( TRACCC_BUILD_CUDA )
  add_subdirectory(cuda)
endif()

# Add stdpar benchmarks
add_subdirectory(stdpar)
