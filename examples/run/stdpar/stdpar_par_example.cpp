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
#include "traccc/clusterization/spacepoint_formation.hpp"

// options
#include "traccc/options/common_options.hpp"
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

int seq_run(const traccc::full_tracking_input_config& i_cfg,
            const traccc::common_options& common_opts) {
    const auto t1 = std::chrono::high_resolution_clock::now();

    // Read the surface transforms
    auto surface_transforms = traccc::read_geometry2(i_cfg.detector_file); // to run successfully using std par, the modified version of read_geometry is needed (related to string deconstruction problem)

    // Read the digitization configuration file
    auto digi_cfg =
        traccc::read_digitization_config(i_cfg.digitization_config_file);

    // Output stats
    uint64_t n_cells = 0;
    uint64_t n_modules = 0;
    uint64_t n_measurements = 0;
    uint64_t n_spacepoints = 0;

    // Memory resource used by the EDM.
    vecmem::host_memory_resource host_mr;  

    traccc::stdpar::clusterization_algorithm ca(host_mr);
    traccc::spacepoint_formation sf(host_mr);

    // Loop over events
    for (unsigned int event = common_opts.skip;
         event < common_opts.events + common_opts.skip; ++event) {

        // Read the cells from the relevant event file
        traccc::cell_container_types::host cells_per_event =
            traccc::read_cells_from_event(
                event, i_cfg.cell_directory, common_opts.input_data_format,
                surface_transforms, digi_cfg, host_mr);

        /*-------------------
            Clusterization
          -------------------*/

        auto measurements_per_event = ca(cells_per_event);

        /*------------------------
            Spacepoint formation
          ------------------------*/

        auto spacepoints_per_event = sf(measurements_per_event);

        std::cout << "----------\n";
        std::cout << "Data of spacepoint for validation:\n";
        std::cout << "x: " << spacepoints_per_event[0].items[0].global[0] << std::endl;
        std::cout << "y: " << spacepoints_per_event[0].items[0].global[1] << std::endl;
        std::cout << "z: " << spacepoints_per_event[0].items[0].global[2] << std::endl;
        std::cout << "----------\n";

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
    traccc::common_options common_opts(desc);
    traccc::full_tracking_input_config full_tracking_input_cfg(desc);

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    // Check errors
    traccc::handle_argument_errors(vm, desc);

    // Read options
    common_opts.read(vm);
    full_tracking_input_cfg.read(vm);

    std::cout << "Running " << argv[0] << " "
              << full_tracking_input_cfg.detector_file << " "
              << full_tracking_input_cfg.cell_directory << " "
              << common_opts.events << std::endl;

    return seq_run(full_tracking_input_cfg, common_opts);
}
