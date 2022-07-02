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


// Time measurements (temporary) // TODO: remove again (20220401)
#include <chrono>
#include <iostream>
#include <stdio.h>

// Benchmark 
#include <benchmark/benchmark.h>
#include <vector>

/*
execute_and_print "TRACCC_TEST_DATA_DIR=/home/nwachuch/bld6/traccc/data
/home/nwachuch/bld6/traccc/build/bin/traccc_stdpar_par_example 
--detector_file=tml_detector/trackml-detector.csv 
--cell_directory=tml_pixels/ 
--events=10
--digitization_config_file=tml_detector/default-geometric-config-generic.json"
*/

// define the cell directories in an indexable manner since google benchmark only accepts integer arguments...
const std::vector<std::string> cell_directories = {
  "tml_full/ttbar_mu20/",
  "tml_full/ttbar_mu40/",
  "tml_full/ttbar_mu60/",
  // "tml_full/ttbar_mu100/",
  // "tml_full/ttbar_mu200/",
  // "tml_full/ttbar_mu300/",
  "tml_pixels/",
};
// define the algorithms that shall be measured
const std::vector<traccc::cuda::cc_algorithm> cc_algorithms = {
  traccc::cuda::cc_algorithm::simplified_sv, 
  traccc::cuda::cc_algorithm::fast_sv_1, 
  traccc::cuda::cc_algorithm::fast_sv_2
};

// generate a sequence from 0 to last index of test files
static void parameter_space(benchmark::internal::Benchmark* b) {
  for (int i = 0; i < cell_directories.size(); ++i){
    for (int j = 0; j < cc_algorithms.size(); ++j){
      b->Args({i, j});
    }
  }
}

static void BM_CudaCCA(benchmark::State& state){
  // configuration for the text instead of CLI args
  const char* detector_file = "tml_detector/trackml-detector.csv";
  const short number_of_events = 10;
  const char* digitization_config_file = "tml_detector/default-geometric-config-generic.json";

  // Read the surface transforms
  auto surface_transforms = traccc::read_geometry(detector_file);

  // Read the digitization configuration file
  auto digi_cfg =
      traccc::read_digitization_config(digitization_config_file);

  // Memory resource used by the EDM.
  vecmem::host_memory_resource host_mr;  

  traccc::cuda::component_connection ca;

  for (auto _ : state){
    // Init manual timer
    double total_elpased_time = 0;
    double kernel_execution_duration;
    // Read params for this iteration of the benchmark
    auto cell_directory = cell_directories.at(state.range(0)); // returns only the file name
    auto cc_algorithm = cc_algorithms.at(state.range(1));
    // Loop over events
    for (unsigned int event = 0;
         event < number_of_events; ++event) {

        // Read the cells from the relevant event file
        traccc::cell_container_types::host cells_per_event =
            traccc::read_cells_from_event(
                event, cell_directory, traccc::data_format::csv,
                surface_transforms, digi_cfg, host_mr);

        /*-------------------
            Clusterization
          -------------------*/
        // Kernel execution time is going to be stored in the kernel_execution_duration variable
        auto measurements_per_event = ca(cells_per_event, &kernel_execution_duration, cc_algorithm);

        // Update manual timer
        total_elpased_time += kernel_execution_duration;
    }

    // Set execution time for this dataset
    state.SetIterationTime(total_elpased_time);
  }
}

BENCHMARK(BM_CudaCCA)
  ->Unit(benchmark::kMillisecond)
  ->Apply(parameter_space)
  ->UseManualTime();

BENCHMARK_MAIN();
