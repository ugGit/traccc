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
#include "traccc/stdpar/utils/CountingIterator.hpp"
#include "traccc/stdpar/clusterization/test.hpp"
#include <algorithm>
#include <execution>

namespace traccc::stdpar {

/// Connected component labeling.
struct spacepoint_formation
    : public algorithm<spacepoint_container_types::host(const measurement_container_types::host&)> {

    public:
    /// Constructor for spacepoint_formation
    ///
    /// @param mr is the memory resource
    spacepoint_formation(vecmem::memory_resource& mr) : m_mr(mr) {}

    // Placeholder to comply with the algorithm interface
    output_type operator()(
        const measurement_container_types::host& measurements) const override {
        output_type empty_result;
        return empty_result;
    }

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
    std::vector<spacepoint> operator()( // TODO the output type is substituted currently to prevent compilation problems.
        const cell_module& c,
        const measurement_container_types::host::item_vector::const_reference& m) const {        
        std::vector<spacepoint> spacepoints;
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
                    const measurement_container_types::host::item_vector::const_reference& measurements,
                    std::vector<spacepoint>& spacepoints) const {
        // Run the algorithm
        int number_of_measurements = measurements.size();
        measurement *measurements_array = new measurement[number_of_measurements];
        spacepoint *spacepoints_array = new spacepoint[number_of_measurements];
        
        // populate the array by copying the values
        for (int i = 0; i < number_of_measurements; i++){
            measurements_array[i] = measurements.at(i);
        }

        std::for_each_n(std::execution::par_unseq, counting_iterator(0), number_of_measurements, 
            [=](unsigned int i){
                const auto m = measurements_array[i];

                point3 local_3d = {m.local[0], m.local[1], 0.};
                point3 global = module.placement.point_to_global(local_3d);
                spacepoint s({global, m});

                spacepoints_array[i] = std::move(s);
            }
        ); 

        // store the values in the output
        spacepoints.reserve(number_of_measurements);
        spacepoints.insert(spacepoints.end(), &spacepoints_array[0], &spacepoints_array[number_of_measurements]);

        // TODO: test when the next steps are included before keeping this
        // delete[] measurements_array;
        // delete[] spacepoints_array;

        // store the values in the output
    }

    private:
    std::reference_wrapper<vecmem::memory_resource> m_mr;
};

}  // namespace traccc::stdpar
