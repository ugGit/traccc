/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "traccc/edm/measurement.hpp"
#include "traccc/edm/spacepoint.hpp"
#include "traccc/utils/algorithm.hpp"

// Time measurements (temporary) // TODO: remove again (20220401)
#include <chrono>
#include <iostream>

namespace traccc {

/// Connected component labeling.
struct spacepoint_formation
    : public algorithm<host_spacepoint_collection(
          const cell_module&, const host_measurement_collection&)> {

    public:
    /// Constructor for spacepoint_formation
    ///
    /// @param mr is the memory resource
    spacepoint_formation(vecmem::memory_resource& mr) : m_mr(mr) {}

    /// Callable operator for the space point formation, based on one single
    /// module
    ///
    /// @param measurements are the input measurements, in this pixel
    /// demonstrator it one space
    ///    point per measurement
    ///
    /// C++20 piping interface
    ///
    /// @return a measurement collection - size of input/output container is
    /// identical
    output_type operator()(
        const cell_module& c,
        const host_measurement_collection& m) const override {
        output_type spacepoints;
        this->operator()(c, m, spacepoints);
        return spacepoints;
    }

    /// Callable operator for the space point formation, based on one single
    /// module
    ///
    /// @param measurements are the input measurements, in this pixel
    /// demonstrator it one space
    ///    point per measurement
    ///
    /// void interface
    ///
    /// @return a measurement collection - size of input/output container is
    /// identical
    void operator()(const cell_module& module,
                    const host_measurement_collection& measurements,
                    output_type& spacepoints) const {
        // Run the algorithm
        spacepoints.reserve(measurements.size());
        
        // start crono
        const auto t1 = std::chrono::high_resolution_clock::now();

        for (const auto& m : measurements) {
            point3 local_3d = {m.local[0], m.local[1], 0.};
            point3 global = module.placement.point_to_global(local_3d);
            variance3 variance = {0, 0, 0};
            spacepoint s({global, variance, m});

            // @todo add variance estimation
            spacepoints.push_back(std::move(s));
        }

        // stop crono
        const auto t2 = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double, std::milli> ms = t2 - t1;
        std::cout << "Execution time [ms]: " << ms.count() << "\n";
        std::cout << "-----------\n";
    }

    private:
    std::reference_wrapper<vecmem::memory_resource> m_mr;
};

}  // namespace traccc
