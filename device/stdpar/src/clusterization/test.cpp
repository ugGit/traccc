/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#include <iostream>
#include <algorithm>
#include <execution>
#include <chrono>

#include <vector>

#include <stdio.h>


#include "traccc/stdpar/clusterization/test.hpp"
#include "traccc/stdpar/utils/CountingIterator.h"

namespace traccc {
namespace stdpar {

using std::cout;

void local_to_global(const cell_module& module, 
                     measurement *measurements_array, 
                     spacepoint *spacepoints_array, 
                     const int number_of_measurements){
    std::for_each_n(std::execution::par, counting_iterator(0), number_of_measurements, 
        [=](unsigned int i){
            const auto m = measurements_array[i];
            
            point3 local_3d = {m.local[0], m.local[1], 0.};
            
            point3 global = module.placement.point_to_global(local_3d);
            variance3 variance = {0, 0, 0};
            spacepoint s({global, variance, m});
            
            // @todo add variance estimation
            // spacepoints_array[i] = s;
            spacepoints_array[i] = {global, variance, m};
        }
    ); 
}

void execute(){
  const int DSIZE = 32; // 2*32*1048576;
  cout << "Start Vector Add Program\n";
  cout << "-----------\n";

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

  // start crono
  const auto t1 = std::chrono::high_resolution_clock::now();

  // execute 
  std::transform(std::execution::par, a.begin(), a.end(), b.begin(), c.begin(), [](float x, float y) -> float {return x+y;});
  
  // stop crono
  const auto t2 = std::chrono::high_resolution_clock::now();
  const std::chrono::duration<double, std::milli> ms = t2 - t1;
  
  cout << "Execution time [ms]: " << ms.count() << "\n";
  cout << "-----------\n";

  cout << a.at(0) << "\n";
  cout << b.at(0) << "\n";
  cout << c.at(0) << "\n";
  
}


}  // namespace stdpar
}  // namespace traccc
