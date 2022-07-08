
#pragma once

#include "traccc/edm/cell.hpp"

// System include(s).
#include <cstddef>

struct cluster_id {
    event_id event = 0;
    std::size_t module_idx = 0;
    geometry_id module = 0;
    scalar threshold = 0.;
    transform3 placement = transform3{};
    pixel_data pixel;
};

struct cluster_element{
  cluster_id header;
  traccc::cell* items;
  unsigned int items_size;
};
