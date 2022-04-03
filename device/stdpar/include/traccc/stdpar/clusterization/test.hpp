/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 LBNL for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "traccc/edm/measurement.hpp"
#include "traccc/edm/spacepoint.hpp"

namespace traccc {
namespace stdpar {

void local_to_global(const cell_module& module, 
                     measurement *measurements_array, 
                     spacepoint *spacepoints_array, 
                     const int number_of_measurements);
void execute();

}  // namespace stdpar
}  // namespace traccc
