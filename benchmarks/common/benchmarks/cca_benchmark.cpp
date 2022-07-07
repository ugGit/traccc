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
#include "traccc/io/data_format.hpp"

// System include(s).
#include <exception>
#include <iostream>

// Benchmark 
#include <benchmark/benchmark.h>
#include <vector>


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

// generate a sequence from 0 to last index of test files
static void parameter_space(benchmark::internal::Benchmark* b) {
  for (unsigned int i = 0; i < cell_directories.size(); ++i){
    b->Arg(i);
  }
}

// caching structure to readuce amount of times where the data file has to be read (I/O)
// the map uses the cell_directory and event id as key, and holds the read in cell container
static std::map<std::tuple<std::string, int>, traccc::cell_container_types::host> cell_data_cache;

// generic function to measure a clusterization algorithm
template<typename T>
static void BM_CA(benchmark::State& state, T&& wrapper_ca_function){
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

  for (auto _ : state){
    // Init manual timer
    double total_elpased_time = 0;
    double kernel_execution_duration;
    // Read params for this iteration of the benchmark
    auto cell_directory = cell_directories.at(state.range(0)); // returns only the file name
    // Loop over events
    for (unsigned int event = 0;
         event < number_of_events; ++event) {

        // Read the cells from the relevant event file or the cache
        traccc::cell_container_types::host cells_per_event;
        auto key = std::make_tuple(cell_directory, event);
        if (cell_data_cache.count(key)){
            cells_per_event = cell_data_cache[key];
        } else {
            cells_per_event =
              traccc::read_cells_from_event(
                  event, cell_directory, traccc::data_format::csv,
                  surface_transforms, digi_cfg, host_mr);
            cell_data_cache[key] = cells_per_event;
        }
        
        /*-------------------
            Clusterization
          -------------------*/
        // Kernel execution time is going to be stored in the kernel_execution_duration variable
        wrapper_ca_function(cells_per_event, &kernel_execution_duration);

        // Update manual timer
        total_elpased_time += kernel_execution_duration;
    }

    // Set execution time for this dataset
    state.SetIterationTime(total_elpased_time);
    state.SetLabel(cell_directory);
  }
}
