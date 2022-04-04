# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2021-2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Project include(s).
include( traccc-compiler-options-cpp )

traccc_add_executable( stdpar_par_example "stdpar_par_example.cpp"
   LINK_LIBRARIES vecmem::core traccc::core traccc::io 
   traccc::options traccc::stdpar Boost::program_options)

set (EXE_NAME "traccc_stdpar_par_example")
set (STDPAR_TARGET "gpu")
target_compile_options(${EXE_NAME} PRIVATE -stdpar=${STDPAR_TARGET})
target_link_options(${EXE_NAME} PRIVATE -stdpar=${STDPAR_TARGET})
