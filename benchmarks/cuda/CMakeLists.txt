cmake_minimum_required(VERSION 3.14)

enable_testing()

find_package(CUDA REQUIRED)
enable_language(CUDA)

traccc_add_executable( benchmark_cuda_cca "cuda_cca.cpp"
    LINK_LIBRARIES traccc::core traccc::io traccc_benchmarks_common
                   traccc::cuda benchmark::benchmark)
