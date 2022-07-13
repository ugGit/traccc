cmake_minimum_required(VERSION 3.14)

enable_testing()

find_package(CUDA REQUIRED)
enable_language(CUDA)

traccc_add_executable( benchmark_stdpar_cca "stdpar_cca.cpp"
    LINK_LIBRARIES traccc::core traccc::io traccc_benchmarks_common
                   traccc::stdpar benchmark::benchmark)

set (LIB "traccc_benchmark_stdpar_cca")
target_compile_options(${LIB} PRIVATE ${TRACCC_STDPAR_BACKEND_FLAG})
target_link_options(${LIB} PRIVATE ${TRACCC_STDPAR_BACKEND_FLAG})
target_compile_options(${LIB} PRIVATE -gpu=managed:intercept)
target_link_options(${LIB} PRIVATE -gpu=managed:intercept)
