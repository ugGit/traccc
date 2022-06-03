/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "traccc/stdpar/clusterization/detail/sparse_ccl.hpp"
#include "traccc/stdpar/clusterization/detail/cluster_element.hpp"
#include "traccc/edm/cell.hpp"
#include "traccc/edm/cluster.hpp"
#include "traccc/utils/algorithm.hpp"

namespace traccc::stdpar {

/// Connected component labelling
///
/// Note that the separation between the public and private interface is
/// only there in the class because the compilers can't automatically figure
/// out the "vector type" of the templated implementation, without adding a
/// lot of "internal knowledge" about the vector types into this piece of
/// code. So instead the public operators are specifically implemented for
/// the host- and device versions of the EDM, making use of a single
/// implementation internally.
///
class component_connection
    : public algorithm<cluster_container_types::host(const cell_container_types::host&)> {
    //   public algorithm<cluster_container_types::host(const device_cell_collection&,
    //                                           const cell_module&)> {
    public:
    /// Constructor for component_connection
    ///
    /// @param mr is the memory resource
    component_connection(vecmem::memory_resource& mr) : m_mr(mr) {}

    // Placeholder since the calling interface for this algorithm differs.
    output_type operator()(
        const cell_container_types::host& cells_per_event) const override {
        printf("The stdpar version of the component_connection algorithm should not be called through this method!");
        output_type empty_result;
        return empty_result;
    }

    /// Implementation of an stdpar version for the component connection. It differs from the original algorith, but this is how it goes
    void operator()(const cell* cells,
                    unsigned int number_of_cells,
                    const cell_module& module,
                    cluster_element*& clusters,
                    unsigned int& num_clusters) {
        // Run the algorithm
        unsigned int* cluster_sizes = new unsigned int[number_of_cells]{}; // initialize values at 0
        unsigned int *connected_cells = new unsigned int[number_of_cells];
        
        detail::sparse_ccl(cells, connected_cells, number_of_cells,
                           num_clusters, cluster_sizes);

        clusters = new cluster_element[num_clusters];
        for(int i = 0; i < num_clusters; i++){
          // initialize the items arrays and store size information
          clusters[i].items = new cell[cluster_sizes[i]];
          clusters[i].items_size = 0; // use it as index when filling the items array later, will correspond at the end to cluster_sizes[i]
        }

        for(int i = 0; i < number_of_cells; i++){
          // get the cluster label info for the current cell
          unsigned int k = connected_cells[i]; 
          // assign add the cell to the cluster, and increase the current index for this cluster
          clusters[k].items[clusters[k].items_size++] = cells[i];
        }

        delete[] connected_cells;
        delete[] cluster_sizes;
    }

    private:
    std::reference_wrapper<vecmem::memory_resource> m_mr;

};  // class component_connection

}  // namespace traccc::stdpar
