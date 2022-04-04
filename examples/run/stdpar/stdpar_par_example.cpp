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

// algorithms
#include "traccc/stdpar/clusterization/test.hpp"
#include "traccc/stdpar/clusterization/clusterization_algorithm.hpp"

// options
#include "traccc/options/full_tracking_input_options.hpp"
#include "traccc/options/handle_argument_errors.hpp"

// System include(s).
#include <exception>
#include <iostream>

// Time measurements (temporary) // TODO: remove again (20220401)
#include <chrono>
#include <iostream>
#include <stdio.h>

namespace po = boost::program_options;


int par_run(const traccc::full_tracking_input_config& i_cfg) {
    // start crono
    const auto t1 = std::chrono::high_resolution_clock::now();

    // Read the surface transforms
    auto surface_transforms = traccc::read_geometry2(i_cfg.detector_file);

    // Output stats
    uint64_t n_cells = 0;
    uint64_t n_modules = 0;
    uint64_t n_measurements = 0;
    uint64_t n_spacepoints = 0;

    // Memory resource used by the EDM.
    vecmem::host_memory_resource host_mr;  

    // check if clusterization algo must be called created on heap or not

    traccc::stdpar::clusterization_algorithm *ca = new traccc::stdpar::clusterization_algorithm(host_mr);

    // Loop over events
    for (unsigned int event = i_cfg.skip; event < i_cfg.events + i_cfg.skip;
         ++event) {

        // Read the cells from the relevant event file
        traccc::host_cell_container cells_per_event =
            traccc::read_cells_from_event(event, i_cfg.cell_directory,
                                          surface_transforms, host_mr);

        /*-------------------
            Clusterization
          -------------------*/

        auto ca_result = (*ca)(cells_per_event);
        auto& measurements_per_event = ca_result.first;
        auto& spacepoints_per_event = ca_result.second;

        /*----------------------------
          Statistics
          ----------------------------*/

        n_modules += cells_per_event.size();
        n_cells += cells_per_event.total_size();
        n_measurements += measurements_per_event.total_size();
        n_spacepoints += spacepoints_per_event.total_size();

        break;
    }

    // stop crono
    const auto t2 = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double, std::milli> ms = t2 - t1;

    std::cout << "==> Statistics ... " << std::endl;
    std::cout << "Execution time [ms]: " << ms.count() << "\n";
    std::cout << "----------    -\n";
    std::cout << "- read    " << n_cells << " cells from " << n_modules
              << " modules" << std::endl;
    std::cout << "- created " << n_measurements << " measurements. "
              << std::endl;
    std::cout << "- created " << n_spacepoints << " space points. "
              << std::endl;      

    return 0;   
}

int main(int argc, char* argv[]){    
    // Set up the program options
    po::options_description desc("Allowed options");

    // Add options
    desc.add_options()("help,h", "Give some help with the program's options");
    traccc::full_tracking_input_config full_tracking_input_cfg(desc);

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    // Read options
    full_tracking_input_cfg.read(vm);

    // Check errors
    traccc::handle_argument_errors(vm, desc);

    std::cout << "Running " << argv[0] << " "
                << full_tracking_input_cfg.detector_file << " "
                << full_tracking_input_cfg.cell_directory << " "
                << full_tracking_input_cfg.events << std::endl;

    return par_run(full_tracking_input_cfg);
}
