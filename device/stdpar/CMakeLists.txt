# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2021-2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Project include(s).
include( traccc-compiler-options-cpp )

# Set up the build of the traccc::stdpar library.
traccc_add_library( traccc_stdpar stdpar TYPE SHARED
  # header files
  "include/traccc/stdpar/clusterization/test.hpp"
  # implementation files
  "src/clusterization/test.cpp" )
target_link_libraries( traccc_stdpar
  PUBLIC traccc::core)

set (STDPAR_TARGET "gpu")
target_compile_options(traccc_stdpar PRIVATE -stdpar=${STDPAR_TARGET})
target_link_options(traccc_stdpar PRIVATE -stdpar=${STDPAR_TARGET})
