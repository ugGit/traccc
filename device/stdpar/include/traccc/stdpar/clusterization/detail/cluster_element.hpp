
#pragma once

#include "traccc/edm/cell.hpp"
#include "traccc/definitions/primitives.hpp"

// System include(s).
#include <cstddef>

struct cluster_id {
    traccc::event_id event = 0;
    std::size_t module_idx = 0;
    traccc::geometry_id module = 0;
    traccc::scalar threshold = 0.;
    traccc::transform3 placement = traccc::transform3{};
    traccc::pixel_data pixel;
};

struct cluster_element{
  cluster_id header;
  traccc::cell* items;
  unsigned int items_size;
};
