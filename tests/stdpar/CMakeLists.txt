# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2021-2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Declare the cpu algorithm test(s).
traccc_add_test( stdpar "test_cca.cpp"
   LINK_LIBRARIES GTest::gtest_main vecmem::core traccc_tests_common
                  traccc::core traccc::io traccc::stdpar)

set (LIB "traccc_test_stdpar")
target_compile_options(${LIB} PRIVATE -stdpar=${STDPAR_TARGET})
target_link_options(${LIB} PRIVATE -stdpar=${STDPAR_TARGET})
target_compile_options(${LIB} PRIVATE -gpu=managed:intercept)
target_link_options(${LIB} PRIVATE -gpu=managed:intercept)
                  