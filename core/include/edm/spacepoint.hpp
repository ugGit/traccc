/** TRACCC library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "definitions/primitives.hpp"
#include "definitions/algebra.hpp"

#include <vector>
#include <limits>

namespace traccc {
    
    /// A cell definition: maximum two channel identifiers
    /// and one activiation value;
    struct spacepoint {
        point3 global = { 0., 0., 0.};
        variance3 variance = { 0., 0., 0.};

        scalar phi = std::numeric_limits<scalar>::quiet_NaN();
        scalar rho = std::numeric_limits<scalar>::quiet_NaN();

        /// Constructor that calculated phi / rho 
        ///
        /// @param global_ global position
        /// @param variance_ variance of the 
        spacepoint(point3&& global_, variance3&& variance_) :
         global(std::move(global_)), variance(std::move(variance_))
        {
            phi = getter::phi(global);
            rho = getter::perp(global);
        }
    };

    struct spacepoint_collection {     

        event_id event = 0;
        geometry_id module = 0;

        std::vector<spacepoint> items = {};
    };

    using spacepoint_container = std::vector<spacepoint_collection>;
}
