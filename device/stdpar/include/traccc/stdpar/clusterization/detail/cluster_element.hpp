#include "traccc/edm/cell.hpp"
#include "traccc/edm/cluster.hpp"

struct cluster_element{
  cluster_id header;
  cell* items;
  unsigned int items_size;
};
