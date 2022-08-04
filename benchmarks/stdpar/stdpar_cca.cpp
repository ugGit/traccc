  /** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// algorithms
#include "traccc/stdpar/clusterization/component_connection_fastsv.hpp"
#include "traccc/stdpar/clusterization/clusterization_algorithm.hpp"

// Benchmark 
#include <benchmark/benchmark.h>
#include "benchmarks/cca_benchmark.cpp"

// generator that returns a wrapper to call the clusterization algorithm with the SparseCCL algorithm
auto stdpar_sparseccl_cca_wrapper_generator(){
  return [=](traccc::cell_container_types::host &cells_per_event, double* kernel_execution_duration){
    vecmem::host_memory_resource host_mr;  
    traccc::stdpar::clusterization_algorithm ca(host_mr);
    ca(cells_per_event, kernel_execution_duration);
  };
}

/*
 * Benchmark stdpar SparseCCL
 */
BENCHMARK_CAPTURE(BM_CCA, cca_stdpar_sparse_ccl, stdpar_sparseccl_cca_wrapper_generator())
  ->Apply(cca_benchmark_config);


static constexpr std::size_t DEFAULT_MIN_CELLS_PER_PARTITION = 1024;

// generator that returns a wrapper to call the clusterization algorithm with the desired cuda cc_algorithm
auto stdpar_fastsv_cca_wrapper_generator(traccc::stdpar::cc_algorithm algo, std::size_t min_cells_per_partition){
  return [=](traccc::cell_container_types::host &cells_per_event, double* kernel_execution_duration){
    traccc::stdpar::component_connection_fastsv ca;
    return ca(cells_per_event, kernel_execution_duration, algo, min_cells_per_partition);
  };
}

/*
 * Benchmark default versions
 */
BENCHMARK_CAPTURE(BM_CCA, cca_stdpar_simplified_sv, 
                  stdpar_fastsv_cca_wrapper_generator(traccc::stdpar::cc_algorithm::simplified_sv, DEFAULT_MIN_CELLS_PER_PARTITION))
  ->Apply(cca_benchmark_config);

BENCHMARK_CAPTURE(BM_CCA, cca_stdpar_fast_sv_1, 
                  stdpar_fastsv_cca_wrapper_generator(traccc::stdpar::cc_algorithm::fast_sv_1, DEFAULT_MIN_CELLS_PER_PARTITION))
  ->Apply(cca_benchmark_config);

BENCHMARK_CAPTURE(BM_CCA, cca_stdpar_fast_sv_2, 
                  stdpar_fastsv_cca_wrapper_generator(traccc::stdpar::cc_algorithm::fast_sv_2, DEFAULT_MIN_CELLS_PER_PARTITION))
  ->Apply(cca_benchmark_config);


/*
 * Benchmark hyper-parameters
 *
 * Note: in theory, we could write a new benchmark that takes an second argument range as input for the partition size. 
 * This would simplify the instantiation, but force us to write a second benchmark that is basically the same as the BM_CCA.
 * Weighting pro's and con's, I decided to take the duplication of code in the instantiation of the code. Especially,
 * because the CUDA version could not use the parametrized benchmark, as the partition size needs to be known at compile time!
 */
BENCHMARK_CAPTURE(BM_CCA, cca_stdpar_fast_sv_2_partition_128, 
                  stdpar_fastsv_cca_wrapper_generator(traccc::stdpar::cc_algorithm::fast_sv_2, 128))
  ->Apply(cca_benchmark_config);

BENCHMARK_CAPTURE(BM_CCA, cca_stdpar_fast_sv_2_partition_256, 
                  stdpar_fastsv_cca_wrapper_generator(traccc::stdpar::cc_algorithm::fast_sv_2, 256))
  ->Apply(cca_benchmark_config);

BENCHMARK_CAPTURE(BM_CCA, cca_stdpar_fast_sv_2_partition_512, 
                  stdpar_fastsv_cca_wrapper_generator(traccc::stdpar::cc_algorithm::fast_sv_2, 512))
  ->Apply(cca_benchmark_config);

BENCHMARK_MAIN();
