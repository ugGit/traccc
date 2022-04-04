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

// STDPAR depencencies
#include "traccc/stdpar/utils/CountingIterator.h"
#include "traccc/stdpar/clusterization/test.hpp"
#include <algorithm>
#include <execution>

#include <iostream>
#include <vector>

namespace traccc::stdpar {

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
        measurement *measurements_array = new measurement[measurements.size()];
        spacepoint *spacepoints_array = new spacepoint[measurements.size()];

        int number_of_measurements = measurements.size();

        // populate the array by copying the values
        for (int i = 0; i < number_of_measurements; i++){
            measurements_array[i] = measurements.at(i);
        }

        std::for_each_n(std::execution::par_unseq, counting_iterator(0), number_of_measurements, 
            [=](unsigned int i){
                const auto m = measurements_array[i];
                point3 local_3d = {m.local[0], m.local[1], 0.};
                
                point3 global = module.placement.point_to_global(local_3d);
                variance3 variance = {0, 0, 0};
                spacepoint s({global, variance, m});
                
                // @todo add variance estimation
                spacepoints_array[i] = s;
            }
        ); 

        // store the values in the output
        for (int i = 0; i < number_of_measurements; i++){
            spacepoints.push_back(spacepoints_array[i]);
        }

        // TODO: test when the next steps are included before keeping this
        // delete[] measurements_array;
        // delete[] spacepoints_array;
    }

    private:
    std::reference_wrapper<vecmem::memory_resource> m_mr;
};

}  // namespace traccc::stdpar
