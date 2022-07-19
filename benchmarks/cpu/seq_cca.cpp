  /** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Algorithms
#include "traccc/clusterization/clusterization_algorithm.hpp"

// Benchmark 
#include <benchmark/benchmark.h>
#include "benchmarks/cca_benchmark.cpp"
#include <chrono>

// Wrapper for coherent usage of benchmark
auto seq_cca_wrapper = [](traccc::cell_container_types::host &cells_per_event, double* kernel_execution_duration){
    // Memory resource used by the EDM.
    vecmem::host_memory_resource host_mr;  

    // Init algorithm
    traccc::clusterization_algorithm ca(host_mr);

    // Start measuring
    const auto start = std::chrono::high_resolution_clock::now();
    
    // Execute CCA
    ca(cells_per_event);

    // End measuring
    const auto end = std::chrono::high_resolution_clock::now();
    auto elpased_time =
      std::chrono::duration_cast<std::chrono::duration<double>>(
        end - start);
    *kernel_execution_duration = elpased_time.count();
};

BENCHMARK_CAPTURE(BM_CCA, seq_cca, seq_cca_wrapper)
  ->Apply(cca_benchmark_config);

BENCHMARK_MAIN();
