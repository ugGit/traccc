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

/*
execute_and_print "TRACCC_TEST_DATA_DIR=/home/nwachuch/bld6/traccc/data
/home/nwachuch/bld6/traccc/build/bin/traccc_stdpar_par_example 
--detector_file=tml_detector/trackml-detector.csv 
--cell_directory=tml_pixels/ 
--events=10
--digitization_config_file=tml_detector/default-geometric-config-generic.json"
*/

// configuration for the text instead of CLI args
const char* detector_file = "tml_detector/trackml-detector.csv";
const char* cell_directory = "tml_pixels/";
const short number_of_events = 10;
const char* digitization_config_file = "tml_detector/default-geometric-config-generic.json";

static void BM_CudaCCA(benchmark::State& state){
  // Read the surface transforms
  auto surface_transforms = traccc::read_geometry(detector_file);

  // Read the digitization configuration file
  auto digi_cfg =
      traccc::read_digitization_config(digitization_config_file);

  // Output stats
  uint64_t n_cells = 0;
  uint64_t n_modules = 0;
  uint64_t n_measurements = 0;
  uint64_t n_spacepoints = 0;

  // Memory resource used by the EDM.
  vecmem::host_memory_resource host_mr;  

  traccc::cuda::component_connection ca;

  for (auto _ : state){
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
        // start crono
        const auto start = std::chrono::high_resolution_clock::now();
        // run algorithm
        auto measurements_per_event = ca(cells_per_event);
        // stop crono
        const auto end = std::chrono::high_resolution_clock::now();
        // register elpased time as iteration duration
        auto elapsed_seconds =
          std::chrono::duration_cast<std::chrono::duration<double>>(
            end - start);
        state.SetIterationTime(elapsed_seconds.count());

        /*----------------------------
          Statistics
          ----------------------------*/
        n_modules += cells_per_event.size();
        n_cells += cells_per_event.total_size();
        n_measurements += measurements_per_event.total_size();
    }
  }

  std::cout << "==> Statistics ... " << std::endl;
  std::cout << "----------    -\n";
  std::cout << "- read    " << n_cells << " cells from " << n_modules
            << " modules" << std::endl;
  std::cout << "- created " << n_measurements << " measurements. "
            << std::endl;
  std::cout << "- created " << n_spacepoints << " space points. "
            << std::endl;     
}

BENCHMARK(BM_CudaCCA)->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();
