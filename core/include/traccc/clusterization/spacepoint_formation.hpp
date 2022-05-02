/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Library include(s).
#include "traccc/edm/measurement.hpp"
#include "traccc/edm/spacepoint.hpp"
#include "traccc/utils/algorithm.hpp"

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>

// System include(s).
#include <functional>

namespace traccc {

/// Algorithm forming space points out of measurements
///
/// This algorithm performs the local-to-global transformation of the 2D
/// measurements made on every detector module, into 3D spacepoint coordinates.
///
class spacepoint_formation : public algorithm<spacepoint_container_types::host(
                                 const measurement_container_types::host&)> {

    public:
    /// Constructor for spacepoint_formation
    ///
    /// @param mr is the memory resource
    ///
    spacepoint_formation(vecmem::memory_resource& mr);

    /// Callable operator for the space point formation, based on one single
    /// module
    ///
    /// @param measurements are the input measurements
    /// @return A spacepoint container, with one spacepoint for every
    ///         measurement
    ///
<<<<<<< HEAD
    output_type operator()(
        const measurement_container_types::host& measurements) const override;
=======
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
>>>>>>> validate results with seq algo

    private:
    std::reference_wrapper<vecmem::memory_resource> m_mr;
};

}  // namespace traccc
