/*
 * TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "traccc/edm/cell.hpp"
#include "traccc/edm/cluster.hpp"
#include "traccc/edm/internal_spacepoint.hpp"
#include "traccc/edm/measurement.hpp"
#include "traccc/edm/spacepoint.hpp"

// clusterization
#include "traccc/stdpar/clusterization/detail/cluster_element.hpp"
#include "traccc/stdpar/clusterization/component_connection_sparseccl.hpp"
#include "traccc/stdpar/clusterization/measurement_creation.hpp"
#include "traccc/stdpar/clusterization/spacepoint_formation.hpp"
#include "traccc/stdpar/clusterization/test.hpp"

#include <iostream>

namespace traccc::stdpar {

class clusterization_algorithm
    : public algorithm<measurement_container_types::host(const cell_container_types::host&)> {

    public:
    /// Constructor for clusterization algorithm
    ///
    /// @param mr is the memory resource
    clusterization_algorithm(vecmem::memory_resource& mr) : m_mr(mr) {
        cc = std::make_shared<traccc::stdpar::component_connection_sparseccl>(
            traccc::stdpar::component_connection_sparseccl(mr));
        mt = std::make_shared<traccc::stdpar::measurement_creation>(
            traccc::stdpar::measurement_creation(mr));
    }

    output_type operator()(
        const cell_container_types::host& cells_per_event) const override {

        // start data transformation
        unsigned int nbr_of_modules = cells_per_event.size();
        cell_module* data_header_array = new cell_module[nbr_of_modules];
        const vecmem::vector<cell>* vecmem_items = cells_per_event.get_items().data();
        cell** data_items_array = new cell*[nbr_of_modules];
        unsigned int* data_items_array_sizes = new unsigned int[nbr_of_modules];
        // copy the entire vector to a structure accessible from within stdpar loop
        for(int i=0; i < nbr_of_modules; i++) {
          data_header_array[i] = cells_per_event.get_headers().at(i);
          data_items_array_sizes[i] = vecmem_items[i].size();
          data_items_array[i] = new cell[vecmem_items[i].size()]; 
          for(int j=0; j < vecmem_items[i].size(); j++) {
            data_items_array[i][j] = vecmem_items[i][j];
          }
        }

        // TODO: parition the problem here, use the algo from CUDA, as this operates on traccc EDM an not flattened arrays

        // reserve as much space as there are modules
        cell_module* output_header_array = new cell_module[nbr_of_modules];
        measurement** output_items_array = new measurement*[nbr_of_modules];
        unsigned int* output_num_measurments_array = new unsigned int[nbr_of_modules];

        // init temporary structure in unified memory to prevent memory space problems when created within stdpar loop
        cluster_element** cluster_container = new cluster_element*[nbr_of_modules]; 
        
        // TODO: should not be necessary, but it is as the arrays can't be initialized within stdpar loop. Default constructs way to many measurements. 
        // init the output_items_array to welcome in the worst case as many measurements as there are activations in the module
        for(int i=0; i < nbr_of_modules; i++){
          const unsigned int nbr_of_activations = cells_per_event.at(i).items.size();
          output_items_array[i] = new measurement[nbr_of_activations];
          cluster_container[i] =  new cluster_element[nbr_of_activations]; 
          for(int j=0; j < nbr_of_activations; j++){
            cluster_container[i][j].items = new traccc::cell[nbr_of_activations]; // TODO: could also be set to `nbr_of_activations-j` to slightly reduce memory usage
          }
        }

        /*
         * Execute the CCA algorithm
         */
        std::for_each_n(std::execution::par, counting_iterator(0), nbr_of_modules, [=](unsigned int i){
          unsigned int num_clusters = 0;

          auto module = data_header_array[i];

          // The algorithmic code part: start
          cc->operator()(data_items_array[i], data_items_array_sizes[i], module, cluster_container[i], num_clusters);

          mt->operator()(cluster_container[i], module, num_clusters, output_items_array[i]);
          // The algorithmnic code part: end

          output_header_array[i] = module; // TODO: check if this is right, because we set placement and pixel to cluster container earlier
          output_num_measurments_array[i] = num_clusters;
        });

        /*
         * Convert data back to expected traccc EDM format
         */
        output_type measurements_per_event(&m_mr.get());
        measurements_per_event.reserve(nbr_of_modules); // reserve enough space
        for(int i=0; i < nbr_of_modules; i++){
          if(output_num_measurments_array[i] == 0) continue;
          // copy array to vector
          vecmem::vector<measurement> items (output_items_array[i], output_items_array[i] + output_num_measurments_array[i]);
          measurements_per_event.push_back(std::move(output_header_array[i]), std::move(items));
        }

        /*
         * Clean up used data structures
         */
        delete[] output_header_array;
        for(int i=0; i < nbr_of_modules; i++){
          delete[] output_items_array[i];
        }
        delete[] output_items_array;
        delete[] output_num_measurments_array;

        return measurements_per_event;
    }

    private:
    // algorithms
    std::shared_ptr<traccc::stdpar::component_connection_sparseccl> cc;
    std::shared_ptr<traccc::stdpar::measurement_creation> mt;
    std::reference_wrapper<vecmem::memory_resource> m_mr;
};

}  // namespace traccc::stdpar
