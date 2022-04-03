/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "traccc/edm/measurement.hpp"
#include "traccc/edm/spacepoint.hpp"
#include "traccc/utils/algorithm.hpp"

// STDPAR depencencies
#include "traccc/stdpar/utils/CountingIterator.h"
#include <algorithm>
#include <execution>

// Time measurements (temporary) // TODO: remove again (20220401)
#include <chrono>
#include <iostream>
#include <stdio.h>


#include <iostream>
#include <algorithm>
#include <execution>
#include <chrono>
#include <vector>

namespace traccc::stdpar {

/// Connected component labeling.
struct spacepoint_formation
    : public algorithm<host_spacepoint_collection(
          const cell_module&, const host_measurement_collection&)> {

    public:
    /// Constructor for spacepoint_formation
    ///
    /// @param mr is the memory resource
    spacepoint_formation(vecmem::memory_resource& mr) : m_mr(mr) {}

    /// Callable operator for the space point formation, based on one single
    /// module
    ///
    /// @param measurements are the input measurements, in this pixel
    /// demonstrator it one space
    ///    point per measurement
    ///
    /// C++20 piping interface
    ///
    /// @return a measurement collection - size of input/output container is
    /// identical
    output_type operator()(
        const cell_module& c,
        const host_measurement_collection& m) const override {        
        output_type spacepoints;
        this->operator()(c, m, spacepoints);
        return spacepoints;
    }

    /// Callable operator for the space point formation, based on one single
    /// module
    ///
    /// @param measurements are the input measurements, in this pixel
    /// demonstrator it one space
    ///    point per measurement
    ///
    /// void interface
    ///
    /// @return a measurement collection - size of input/output container is
    /// identical
    void operator()(const cell_module& module,
                    const host_measurement_collection& measurements,
                    output_type& spacepoints) const {
        // Run the algorithm

        // start crono
        const auto t1 = std::chrono::high_resolution_clock::now();

        spacepoints.reserve(measurements.size());
        measurement *measurements_array = new measurement[measurements.size()];
        spacepoint *spacepoints_array = new spacepoint[measurements.size()];


        const int DSIZE = 32; // 2*32*1048576;
  
  // Try to access elements first with vectors,
  // and in a second round by random access in arrays

  // Initialize vectors
  std::vector<float> a(DSIZE);
  std::vector<float> b(DSIZE);
  std::vector<float> c(DSIZE);
  for (int i = 0; i < DSIZE; i++){
    a.at(i) = rand()/(float)RAND_MAX;
    b.at(i) = rand()/(float)RAND_MAX;
  }

  // execute 
  std::transform(std::execution::par, a.begin(), a.end(), b.begin(), c.begin(), [](float x, float y) -> float {return x+y;});
  
  // stop crono
  const auto t2 = std::chrono::high_resolution_clock::now();
  const std::chrono::duration<double, std::milli> ms = t2 - t1;
 /*       
        // populate the array by copying the values
        for (int i = 0; i < measurements.size(); i++){
            measurements_array[i] = measurements.at(i);
        }

        // TODO: Problem does no come from the loop. It's the compilation of the file itself.
        // possibly from the counting iterator, probably from somewhere else.

        std::for_each_n(std::execution::par, counting_iterator(0), 2, //  measurements.size(), 
            [=](unsigned int i){
                printf("tst");
                // const auto m = measurements_array[i];
                /*
                point3 local_3d = {m.local[0], m.local[1], 0.};
                point3 global = module.placement.point_to_global(local_3d);
                variance3 variance = {0, 0, 0};
                spacepoint s({global, variance, m});

                // @todo add variance estimation
                spacepoints_array[i] = s;
            }
        ); 

        // stop crono
        const auto t2 = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double, std::milli> ms = t2 - t1;
        // std::cout << "Execution time [ms]: " << ms.count() << "\n";
        // std::cout << "-----------\n";
                              */ 
    }

    private:
    std::reference_wrapper<vecmem::memory_resource> m_mr;
};

}  // namespace traccc::stdpar
