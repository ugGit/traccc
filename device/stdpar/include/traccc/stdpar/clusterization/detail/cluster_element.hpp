
#pragma once

#include "traccc/edm/cell.hpp"
#include "traccc/edm/cluster.hpp"

struct cluster_element{
  traccc::cluster_id header;
  traccc::cell* items;
  unsigned int items_size;
};
