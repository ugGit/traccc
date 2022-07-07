  /** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// io
#include "traccc/io/csv.hpp"
#include "traccc/io/reader.hpp"
#include "traccc/io/utils.hpp"
#include "traccc/io/writer.hpp"
#include "traccc/io/data_format.hpp"

// algorithms
#include "traccc/cuda/cca/component_connection.hpp"

// System include(s).
#include <exception>
#include <iostream>

// Benchmark 
#include <benchmark/benchmark.h>
#include "../common/cca_benchmark.cpp"


// wrapper to add the algorithm variant to the function
auto cuda_cca_fast_sv_1 = [](traccc::cell_container_types::host &cells_per_event, double* kernel_execution_duration){
  traccc::cuda::component_connection ca;
  return ca(cells_per_event, kernel_execution_duration, traccc::cuda::cc_algorithm::fast_sv_1);
};

// generator that returns a wrapper to call the clusterization algorithm with the desired cuda cc_algorithm
auto cuda_cca_wrapper_generator(traccc::cuda::cc_algorithm algo){
  return [=](traccc::cell_container_types::host &cells_per_event, double* kernel_execution_duration){
    traccc::cuda::component_connection ca;
    return ca(cells_per_event, kernel_execution_duration, algo);
  };
}

BENCHMARK_CAPTURE(BM_CA, my_name, cuda_cca_wrapper_generator(traccc::cuda::cc_algorithm::fast_sv_1))
  ->Unit(benchmark::kMillisecond)
  ->Apply(parameter_space)
  ->UseManualTime();

BENCHMARK_MAIN();
