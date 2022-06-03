
#pragma once

#include "traccc/edm/cell.hpp"

// System include(s).
#include <cstddef>

struct cluster_element{
  std::size_t header;
  traccc::cell* items;
  unsigned int items_size;
};
