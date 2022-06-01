# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2021-2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Declare the cpu algorithm test(s).
traccc_add_test( stdpar_fastsv "test_cca_fastsv.cpp"
LINK_LIBRARIES GTest::gtest_main vecmem::core traccc_tests_common
               traccc::core traccc::io traccc::stdpar)

traccc_add_test( stdpar_sparseccl "test_cca_sparseccl.cpp"
   LINK_LIBRARIES GTest::gtest_main vecmem::core traccc_tests_common
                  traccc::core traccc::io traccc::stdpar)

set (LIB "traccc_test_stdpar_fastsv")
target_compile_options(${LIB} PRIVATE -stdpar=${STDPAR_TARGET})
target_link_options(${LIB} PRIVATE -stdpar=${STDPAR_TARGET})
target_compile_options(${LIB} PRIVATE -gpu=managed:intercept)
target_link_options(${LIB} PRIVATE -gpu=managed:intercept)

set (LIB2 "traccc_test_stdpar_sparseccl")
target_compile_options(${LIB2} PRIVATE -stdpar=${STDPAR_TARGET})
target_link_options(${LIB2} PRIVATE -stdpar=${STDPAR_TARGET})
target_compile_options(${LIB2} PRIVATE -gpu=managed:intercept)
target_link_options(${LIB2} PRIVATE -gpu=managed:intercept)
                  