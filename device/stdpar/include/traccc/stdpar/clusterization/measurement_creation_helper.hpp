/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/common.hpp"
#include "traccc/definitions/primitives.hpp"
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/cluster.hpp"

namespace traccc::stdpar {

/// Function used for retrieving the cell signal based on the module id
inline scalar signal_cell_modelling(scalar signal_in,
                                    const cell_module& /*cl_id*/) {
    return signal_in;
}

/// Function for pixel segmentation
inline vector2 position_from_cell(cell c, cell_module cl_id) {
    // Retrieve the specific values based on module idx
    return {c.channel0,
            c.channel1};
}

inline vector2 get_pitch(cell_module cl_id) {
    // return the values based on the module idx
    return {cl_id.pixel.pitch_x, cl_id.pixel.pitch_y};
}

/// Function used for calculating the properties of the cluster inside
/// measurement creation
template <template <typename> class vector_type, typename cell_t>
inline void calc_cluster_properties(
    cell* cluster, unsigned int cluster_size, const cell_module& cl_id, point2& mean,
    point2& var, scalar& totalWeight) {
    for(int j=0; j < cluster_size; j++){
        const auto& cell = cluster[j];
        scalar weight = signal_cell_modelling(cell.activation, cl_id);
        if (weight > cl_id.threshold) {
            totalWeight += cell.activation;
            const point2 cell_position = position_from_cell(cell, cl_id);
            const point2 prev = mean;
            const point2 diff = cell_position - prev;

            mean = prev + (weight / totalWeight) * diff;

            var[0] = var[0] + weight * (diff[0]) * (cell_position[0] - mean[0]);
            var[1] = var[1] + weight * (diff[1]) * (cell_position[1] - mean[1]);
        }
    }
}

}  // namespace traccc
