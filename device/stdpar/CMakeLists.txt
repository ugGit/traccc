# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2021-2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

cmake_minimum_required(VERSION 3.18)

# Project include(s).
include( traccc-compiler-options-cpp )

# Set library name
set( LIB "traccc_stdpar" )

# Set up the build of the traccc::stdpar library.
traccc_add_library( ${LIB} stdpar TYPE SHARED
  # utils
  "include/traccc/stdpar/utils/CountingIterator.hpp"
  # header files
  "include/traccc/stdpar/clusterization/test.hpp"
  "include/traccc/stdpar/clusterization/spacepoint_formation.hpp"
  "include/traccc/stdpar/clusterization/measurement_creation.hpp"
  "include/traccc/stdpar/clusterization/component_connection_sparseccl.hpp"
  "include/traccc/stdpar/clusterization/component_connection_fastsv.hpp"
  "include/traccc/stdpar/clusterization/clusterization_algorithm.hpp"
  "include/traccc/stdpar/clusterization/detail/sparse_ccl.hpp"
  "include/traccc/stdpar/clusterization/detail/cluster_element.hpp"
  # implementation files
  "src/clusterization/test.cpp"
  "src/clusterization/component_connection_fastsv.cpp"
)
target_link_libraries( ${LIB}
  PUBLIC traccc::core)

target_compile_options(${LIB} PRIVATE ${TRACCC_STDPAR_BACKEND_FLAG})
target_link_options(${LIB} PRIVATE ${TRACCC_STDPAR_BACKEND_FLAG})
