cmake_minimum_required(VERSION 3.14)

enable_testing()

traccc_add_library( traccc_benchmarks benchmarks TYPE SHARED
    "cca_benchmark.cpp")

target_link_libraries( traccc_benchmarks
PUBLIC traccc::io traccc::core traccc::options benchmark::benchmark )
