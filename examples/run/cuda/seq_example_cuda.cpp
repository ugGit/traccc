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
#include "traccc/clusterization/clusterization_algorithm.hpp"
#include "traccc/clusterization/spacepoint_formation.hpp"
#include "traccc/cuda/seeding/seeding_algorithm.hpp"
#include "traccc/cuda/seeding/track_params_estimation.hpp"
#include "traccc/seeding/seeding_algorithm.hpp"
#include "traccc/seeding/track_params_estimation.hpp"

// performance
#include "traccc/efficiency/seeding_performance_writer.hpp"

// options
#include "traccc/options/common_options.hpp"
#include "traccc/options/full_tracking_input_options.hpp"
#include "traccc/options/handle_argument_errors.hpp"

// vecmem
#include <vecmem/memory/cuda/managed_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s).
#include <chrono>
#include <exception>
#include <iomanip>
#include <iostream>

namespace po = boost::program_options;

int seq_run(const traccc::full_tracking_input_config& i_cfg,
            const traccc::common_options& common_opts, bool run_cpu) {

    // Read the surface transforms
    auto surface_transforms = traccc::read_geometry(i_cfg.detector_file);

    // Read the digitization configuration file
    auto digi_cfg =
        traccc::read_digitization_config(i_cfg.digitization_config_file);

    // Output stats
    uint64_t n_cells = 0;
    uint64_t n_modules = 0;
    // uint64_t n_clusters = 0;
    uint64_t n_measurements = 0;
    uint64_t n_spacepoints = 0;
    uint64_t n_seeds = 0;
    uint64_t n_seeds_cuda = 0;

    // Elapsed time
    float wall_time(0);
    float file_reading_cpu(0);
    float clusterization_cpu(0);
    float sp_formation_cpu(0);
    float seeding_cpu(0);
    // float clusterization_cuda(0);
    float seeding_cuda(0);
    float tp_estimating_cpu(0);
    float tp_estimating_cuda(0);

    // Memory resource used by the EDM.
    vecmem::host_memory_resource host_mr;
    vecmem::cuda::managed_memory_resource mng_mr;

    traccc::clusterization_algorithm ca(mng_mr);
    traccc::spacepoint_formation sf(mng_mr);
    traccc::seeding_algorithm sa(host_mr);
    traccc::track_params_estimation tp(host_mr);

    traccc::cuda::seeding_algorithm sa_cuda(mng_mr);
    traccc::cuda::track_params_estimation tp_cuda(mng_mr);

    // performance writer
    traccc::seeding_performance_writer sd_performance_writer(
        traccc::seeding_performance_writer::config{});
    sd_performance_writer.add_cache("CPU");
    sd_performance_writer.add_cache("CUDA");

    /*time*/ auto start_wall_time = std::chrono::system_clock::now();

    // Loop over events
    for (unsigned int event = common_opts.skip;
         event < common_opts.events + common_opts.skip; ++event) {

        /*time*/ auto start_file_reading_cpu = std::chrono::system_clock::now();

        // Read the cells from the relevant event file
        traccc::cell_container_types::host cells_per_event =
            traccc::read_cells_from_event(
                event, i_cfg.cell_directory, common_opts.input_data_format,
                surface_transforms, digi_cfg, host_mr);

        /*time*/ auto end_file_reading_cpu = std::chrono::system_clock::now();
        /*time*/ std::chrono::duration<double> time_file_reading_cpu =
            end_file_reading_cpu - start_file_reading_cpu;
        /*time*/ file_reading_cpu += time_file_reading_cpu.count();

        /*-----------------------------
              Clusterization (cpu)
          -----------------------------*/

        /*time*/ auto start_clusterization_cpu =
            std::chrono::system_clock::now();

        auto measurements_per_event = ca(cells_per_event);

        /*time*/ auto end_clusterization_cpu = std::chrono::system_clock::now();
        /*time*/ std::chrono::duration<double> time_clusterization_cpu =
            end_clusterization_cpu - start_clusterization_cpu;
        /*time*/ clusterization_cpu += time_clusterization_cpu.count();

        /*---------------------------------
               Spacepoint formation (cpu)
          ---------------------------------*/

        /*time*/ auto start_sp_formation_cpu = std::chrono::system_clock::now();

        auto spacepoints_per_event = sf(measurements_per_event);

        /*time*/ auto end_sp_formation_cpu = std::chrono::system_clock::now();
        /*time*/ std::chrono::duration<double> time_sp_formation_cpu =
            end_sp_formation_cpu - start_sp_formation_cpu;
        /*time*/ sp_formation_cpu += time_sp_formation_cpu.count();

        /*----------------------------
             Seeding algorithm
          ----------------------------*/

        // CUDA

        /*time*/ auto start_seeding_cuda = std::chrono::system_clock::now();

        auto seeds_cuda = sa_cuda(std::move(spacepoints_per_event));

        /*time*/ auto end_seeding_cuda = std::chrono::system_clock::now();
        /*time*/ std::chrono::duration<double> time_seeding_cuda =
            end_seeding_cuda - start_seeding_cuda;
        /*time*/ seeding_cuda += time_seeding_cuda.count();

        // CPU

        /*time*/ auto start_seeding_cpu = std::chrono::system_clock::now();

        traccc::seeding_algorithm::output_type seeds;

        if (run_cpu) {
            seeds = sa(spacepoints_per_event);
        }

        /*time*/ auto end_seeding_cpu = std::chrono::system_clock::now();
        /*time*/ std::chrono::duration<double> time_seeding_cpu =
            end_seeding_cpu - start_seeding_cpu;
        /*time*/ seeding_cpu += time_seeding_cpu.count();

        /*----------------------------
          Track params estimation
          ----------------------------*/

        // CUDA

        /*time*/ auto start_tp_estimating_cuda =
            std::chrono::system_clock::now();

        auto params_cuda =
            tp_cuda(std::move(spacepoints_per_event), std::move(seeds_cuda));

        /*time*/ auto end_tp_estimating_cuda = std::chrono::system_clock::now();
        /*time*/ std::chrono::duration<double> time_tp_estimating_cuda =
            end_tp_estimating_cuda - start_tp_estimating_cuda;
        /*time*/ tp_estimating_cuda += time_tp_estimating_cuda.count();

        // CPU

        /*time*/ auto start_tp_estimating_cpu =
            std::chrono::system_clock::now();

        traccc::track_params_estimation::output_type params;
        if (run_cpu) {
            params = tp(spacepoints_per_event, seeds);
        }

        /*time*/ auto end_tp_estimating_cpu = std::chrono::system_clock::now();
        /*time*/ std::chrono::duration<double> time_tp_estimating_cpu =
            end_tp_estimating_cpu - start_tp_estimating_cpu;
        /*time*/ tp_estimating_cpu += time_tp_estimating_cpu.count();

        /*----------------------------------
          compare cpu and cuda result
          ----------------------------------*/

        if (run_cpu) {
            // seeding
            int n_match = 0;

            std::vector<std::array<traccc::spacepoint, 3>> sp3_vector =
                traccc::get_spacepoint_vector(seeds, spacepoints_per_event);

            std::vector<std::array<traccc::spacepoint, 3>> sp3_vector_cuda =
                traccc::get_spacepoint_vector(seeds_cuda,
                                              spacepoints_per_event);

            for (const auto& sp3 : sp3_vector) {
                if (std::find(sp3_vector_cuda.cbegin(), sp3_vector_cuda.cend(),
                              sp3) != sp3_vector_cuda.cend()) {
                    n_match++;
                }
            }

            float matching_rate = float(n_match) / seeds.size();
            std::cout << "event " << std::to_string(event) << std::endl;
            std::cout << " number of seeds (cpu): " << seeds.size()
                      << std::endl;
            std::cout << " number of seeds (cuda): " << seeds_cuda.size()
                      << std::endl;
            std::cout << " seed matching rate: " << matching_rate << std::endl;

            // track parameter estimation
            n_match = 0;
            for (auto& param : params) {
                if (std::find(params_cuda.begin(), params_cuda.end(), param) !=
                    params_cuda.end()) {
                    n_match++;
                }
            }
            matching_rate = float(n_match) / params.size();
            std::cout << " track parameters matching rate: " << matching_rate
                      << std::endl;
        }

        /*----------------
             Statistics
          ---------------*/

        n_modules += cells_per_event.size();
        n_cells += cells_per_event.total_size();
        n_measurements += measurements_per_event.total_size();
        n_spacepoints += spacepoints_per_event.total_size();
        n_seeds_cuda += seeds_cuda.size();
        n_seeds += seeds.size();

        /*------------
             Writer
          ------------*/

        if (i_cfg.check_seeding_performance) {
            traccc::event_map evt_map(event, i_cfg.detector_file,
                                      i_cfg.digitization_config_file,
                                      i_cfg.cell_directory, i_cfg.hit_directory,
                                      i_cfg.particle_directory, host_mr);
            sd_performance_writer.write("CUDA", seeds_cuda,
                                        spacepoints_per_event, evt_map);

            if (run_cpu) {
                sd_performance_writer.write("CPU", seeds, spacepoints_per_event,
                                            evt_map);
            }
        }
    }

    /*time*/ auto end_wall_time = std::chrono::system_clock::now();
    /*time*/ std::chrono::duration<double> time_wall_time =
        end_wall_time - start_wall_time;

    /*time*/ wall_time += time_wall_time.count();

    sd_performance_writer.finalize();

    std::cout << "==> Statistics ... " << std::endl;
    std::cout << "- read    " << n_spacepoints << " spacepoints from "
              << n_modules << " modules" << std::endl;
    std::cout << "- created        " << n_cells << " cells           "
              << std::endl;
    std::cout << "- created        " << n_measurements << " meaurements     "
              << std::endl;
    std::cout << "- created        " << n_spacepoints << " spacepoints     "
              << std::endl;

    std::cout << "- created (cpu)  " << n_seeds << " seeds" << std::endl;
    std::cout << "- created (cuda) " << n_seeds_cuda << " seeds" << std::endl;
    std::cout << "==> Elpased time ... " << std::endl;
    std::cout << "wall time           " << std::setw(10) << std::left
              << wall_time << std::endl;
    std::cout << "file reading (cpu)        " << std::setw(10) << std::left
              << file_reading_cpu << std::endl;
    std::cout << "clusterization_time (cpu) " << std::setw(10) << std::left
              << clusterization_cpu << std::endl;
    std::cout << "spacepoint_formation_time (cpu) " << std::setw(10)
              << std::left << sp_formation_cpu << std::endl;
    std::cout << "seeding_time (cpu)        " << std::setw(10) << std::left
              << seeding_cpu << std::endl;
    std::cout << "seeding_time (cuda)       " << std::setw(10) << std::left
              << seeding_cuda << std::endl;
    std::cout << "tr_par_esti_time (cpu)    " << std::setw(10) << std::left
              << tp_estimating_cpu << std::endl;
    std::cout << "tr_par_esti_time (cuda)   " << std::setw(10) << std::left
              << tp_estimating_cuda << std::endl;

    return 0;
}

// The main routine
//
int main(int argc, char* argv[]) {
    // Set up the program options
    po::options_description desc("Allowed options");

    // Add options
    desc.add_options()("help,h", "Give some help with the program's options");
    traccc::common_options common_opts(desc);
    traccc::full_tracking_input_config full_tracking_input_cfg(desc);
    desc.add_options()("run_cpu", po::value<bool>()->default_value(false),
                       "run cpu tracking as well");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    // Check errors
    traccc::handle_argument_errors(vm, desc);

    // Read options
    common_opts.read(vm);
    full_tracking_input_cfg.read(vm);
    auto run_cpu = vm["run_cpu"].as<bool>();

    std::cout << "Running " << argv[0] << " "
              << full_tracking_input_cfg.detector_file << " "
              << full_tracking_input_cfg.cell_directory << " "
              << common_opts.events << std::endl;

    return seq_run(full_tracking_input_cfg, common_opts, run_cpu);
}
