cmake_minimum_required(VERSION 3.14)

enable_testing()

traccc_add_executable( benchmark_seq_cca "seq_cca.cpp"
    LINK_LIBRARIES traccc::core traccc::io traccc_benchmarks_common benchmark::benchmark)
