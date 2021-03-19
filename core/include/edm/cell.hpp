/** TRACCC library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "definitions/algebra.hpp"
#include "definitions/primitives.hpp"
#include <vector>
#include <limits>

namespace traccc {

    using channel_id = unsigned int;

    /// A cell definition: 
    ///
    /// maximum two channel identifiers
    /// and one activiation value, such as a time stamp
    struct cell {
        channel_id channel0 = std::numeric_limits<channel_id>::max();
        channel_id channel1 = std::numeric_limits<channel_id>::max();
        scalar activation = std::numeric_limits<scalar>::quiet_NaN();
        scalar time = std::numeric_limits<scalar>::quiet_NaN();
    };

    /// A cell collection: 
    ///
    /// it remembers the moduleentifier and also 
    /// keeps track of the cell ranges for chosing optimal
    /// algorithm.
    struct cell_collection { 

        event_id event = 0;
        geometry_id module = 0;
        transform3 placement = transform3{};

        std::vector<cell> items = {};
        std::array<channel_id,2> range0 = {std::numeric_limits<channel_id>::max(), 0};
        std::array<channel_id,2> range1 = {std::numeric_limits<channel_id>::max(), 0};         
    };

    using cell_container = std::vector<cell_collection>;

}

