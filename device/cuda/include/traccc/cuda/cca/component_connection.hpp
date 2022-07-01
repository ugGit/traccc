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

#include <benchmark/benchmark.h>

namespace traccc::cuda {
struct component_connection : algorithm<measurement_container_types::host(
                                  const cell_container_types::host& cells)> {
    // pass the benchmark state as nullptr to detect when none is passed, marks a non-breaking code extension
    output_type operator()(const cell_container_types::host& cells) const{
      return this->operator()(cells, nullptr);
    };
    output_type operator()(const cell_container_types::host& cells, benchmark::State *state) const;
};
}  // namespace traccc::cuda
