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

namespace traccc {

/// Connected component labeling.
struct spacepoint_formation : public algorithm<host_spacepoint_container(
                                  const host_measurement_container&)> {

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
    output_type operator()(const host_measurement_container&
                               measurements_per_event) const override {
        output_type spacepoints_per_event(measurements_per_event.size(),
                                          &m_mr.get());

        // Run the algorithm
        spacepoints.reserve(measurements.size());

        for (const auto& m : measurements) {
            point3 local_3d = {m.local[0], m.local[1], 0.};
            point3 global = module.placement.point_to_global(local_3d);
            variance3 variance = {0, 0, 0};
            spacepoint s({global, variance, m});
            spacepoints_per_module.push_back(std::move(s));
            }
        }
<<<<<<< HEAD

        return spacepoints_per_event;
=======
>>>>>>> validate results with seq algo
    }

    private:
    std::reference_wrapper<vecmem::memory_resource> m_mr;
};

}  // namespace traccc
