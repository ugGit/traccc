# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2021-2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Project include(s).
include( traccc-compiler-options-cpp )

traccc_add_executable( stdpar_sparseccl_example "stdpar_sparseccl.cpp"
   LINK_LIBRARIES vecmem::core traccc::core traccc::io 
   traccc::options traccc::stdpar Boost::program_options)

set (EXE_NAME "traccc_stdpar_sparseccl_example")
set (STDPAR_TARGET "gpu")
target_compile_options(${EXE_NAME} PRIVATE ${TRACCC_STDPAR_BACKEND_FLAG})
target_link_options(${EXE_NAME} PRIVATE ${TRACCC_STDPAR_BACKEND_FLAG})
target_compile_options(${EXE_NAME} PRIVATE -gpu=managed:intercept)
target_link_options(${EXE_NAME} PRIVATE -gpu=managed:intercept)



traccc_add_executable( stdpar_fastsv_example "stdpar_fastsv.cpp"
   LINK_LIBRARIES vecmem::core traccc::core traccc::io 
   traccc::options traccc::stdpar Boost::program_options)

set (EXE_NAME2 "traccc_stdpar_fastsv_example")
set (STDPAR_TARGET "gpu")
target_compile_options(${EXE_NAME2} PRIVATE ${TRACCC_STDPAR_BACKEND_FLAG})
target_link_options(${EXE_NAME2} PRIVATE ${TRACCC_STDPAR_BACKEND_FLAG})
target_compile_options(${EXE_NAME2} PRIVATE -gpu=managed:intercept)
target_link_options(${EXE_NAME2} PRIVATE -gpu=managed:intercept)
