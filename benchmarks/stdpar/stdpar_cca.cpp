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

// generator that returns a wrapper to call the clusterization algorithm with the desired cuda cc_algorithm
auto stdpar_cca_wrapper_generator(traccc::stdpar::cc_algorithm algo){
  return [=](traccc::cell_container_types::host &cells_per_event, double* kernel_execution_duration){
    traccc::stdpar::component_connection_fastsv ca;
    return ca(cells_per_event, kernel_execution_duration, algo);
  };
}

BENCHMARK_CAPTURE(BM_CCA, cca_stdpar_simplified_sv, stdpar_cca_wrapper_generator(traccc::stdpar::cc_algorithm::simplified_sv))
  ->Unit(benchmark::kMillisecond)
  ->Apply(parameter_space)
  ->ArgNames({"DatasetFileIndex"})
  ->UseManualTime();

BENCHMARK_CAPTURE(BM_CCA, cca_stdpar_fast_sv_1, stdpar_cca_wrapper_generator(traccc::stdpar::cc_algorithm::fast_sv_1))
  ->Unit(benchmark::kMillisecond)
  ->Apply(parameter_space)
  ->ArgNames({"DatasetFileIndex"})
  ->UseManualTime();

BENCHMARK_CAPTURE(BM_CCA, cca_stdpar_fast_sv_2, stdpar_cca_wrapper_generator(traccc::stdpar::cc_algorithm::fast_sv_2))
  ->Unit(benchmark::kMillisecond)
  ->Apply(parameter_space)
  ->ArgNames({"DatasetFileIndex"})
  ->UseManualTime();

BENCHMARK_MAIN();
