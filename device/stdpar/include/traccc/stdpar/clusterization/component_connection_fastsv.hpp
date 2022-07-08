/*
 * TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "traccc/edm/cell.hpp"
#include "traccc/edm/measurement.hpp"
#include "traccc/utils/algorithm.hpp"

namespace traccc::stdpar {

enum class cc_algorithm{ simplified_sv, fast_sv_1, fast_sv_2};

struct component_connection_fastsv : algorithm<measurement_container_types::host(
                                  const cell_container_types::host& cells)> {
    // pass the benchmark state as nullptr to detect when none is passed, marks a non-breaking code extension
    output_type operator()(const cell_container_types::host& cells) const{
      return this->operator()(cells, nullptr, cc_algorithm::simplified_sv);
    };
    
    output_type operator()(const cell_container_types::host& cells, 
                           double* kernel_execution_duration,
                           cc_algorithm selected_algorithm) const;
};
}  // namespace traccc::stdpar
