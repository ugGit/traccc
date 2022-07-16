
#pragma once

#include "traccc/edm/cell.hpp"
#include "traccc/definitions/primitives.hpp"

// System include(s).
#include <cstddef>

struct cluster_element{
  traccc::cell_module header;
  traccc::cell* items;
  unsigned int items_size;
};
