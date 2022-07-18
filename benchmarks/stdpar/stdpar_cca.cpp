  /** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// algorithms
#include "traccc/stdpar/clusterization/component_connection_fastsv.hpp"

// Benchmark 
#include <benchmark/benchmark.h>
#include "benchmarks/cca_benchmark.cpp"

static constexpr std::size_t DEFAULT_MIN_CELLS_PER_PARTITION = 1024;

// generator that returns a wrapper to call the clusterization algorithm with the desired cuda cc_algorithm
auto stdpar_cca_wrapper_generator(traccc::stdpar::cc_algorithm algo, std::size_t min_cells_per_partition){
  return [=](traccc::cell_container_types::host &cells_per_event, double* kernel_execution_duration){
    traccc::stdpar::component_connection_fastsv ca;
    return ca(cells_per_event, kernel_execution_duration, algo, min_cells_per_partition);
  };
}

/*
 * Benchmark default versions
 */
BENCHMARK_CAPTURE(BM_CCA, cca_stdpar_simplified_sv, 
                  stdpar_cca_wrapper_generator(traccc::stdpar::cc_algorithm::simplified_sv, DEFAULT_MIN_CELLS_PER_PARTITION))
  ->Unit(benchmark::kMillisecond)
  ->Apply(parameter_space)
  ->ArgNames({"DatasetFileIndex"})
  ->UseManualTime();

BENCHMARK_CAPTURE(BM_CCA, cca_stdpar_fast_sv_1, 
                  stdpar_cca_wrapper_generator(traccc::stdpar::cc_algorithm::fast_sv_1, DEFAULT_MIN_CELLS_PER_PARTITION))
  ->Unit(benchmark::kMillisecond)
  ->Apply(parameter_space)
  ->ArgNames({"DatasetFileIndex"})
  ->UseManualTime();

BENCHMARK_CAPTURE(BM_CCA, cca_stdpar_fast_sv_2, 
                  stdpar_cca_wrapper_generator(traccc::stdpar::cc_algorithm::fast_sv_2, DEFAULT_MIN_CELLS_PER_PARTITION))
  ->Unit(benchmark::kMillisecond)
  ->Apply(parameter_space)
  ->ArgNames({"DatasetFileIndex"})
  ->UseManualTime();


/*
 * Benchmark hyper-parameters
 *
 * Note: in theory, we could write a new benchmark that takes an second argument range as input for the partition size. 
 * This would simplify the instantiation, but force us to write a second benchmark that is basically the same as the BM_CCA.
 * Weighting pro's and con's, I decided to take the duplication of code in the instantiation of the code. Especially,
 * because the CUDA version could not use the parametrized benchmark, as the partition size needs to be known at compile time!
 */
BENCHMARK_CAPTURE(BM_CCA, cca_stdpar_fast_sv_2_partition_128, 
                  stdpar_cca_wrapper_generator(traccc::stdpar::cc_algorithm::fast_sv_2, 128))
  ->Unit(benchmark::kMillisecond)
  ->Apply(parameter_space)
  ->ArgNames({"DatasetFileIndex"})
  ->UseManualTime();

BENCHMARK_CAPTURE(BM_CCA, cca_stdpar_fast_sv_2_partition_256, 
                  stdpar_cca_wrapper_generator(traccc::stdpar::cc_algorithm::fast_sv_2, 256))
  ->Unit(benchmark::kMillisecond)
  ->Apply(parameter_space)
  ->ArgNames({"DatasetFileIndex"})
  ->UseManualTime();

BENCHMARK_CAPTURE(BM_CCA, cca_stdpar_fast_sv_2_partition_512, 
                  stdpar_cca_wrapper_generator(traccc::stdpar::cc_algorithm::fast_sv_2, 512))
  ->Unit(benchmark::kMillisecond)
  ->Apply(parameter_space)
  ->ArgNames({"DatasetFileIndex"})
  ->UseManualTime();

BENCHMARK_MAIN();
