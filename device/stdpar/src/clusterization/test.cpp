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


namespace traccc {
namespace stdpar {

using std::cout;

void execute(){
  const int DSIZE = 32; // 2*32*1048576;
  /*
  cout << "Start Vector Add Program\n";
  cout << "-----------\n";
  */

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
  std::transform(std::execution::par_unseq, a.begin(), a.end(), b.begin(), c.begin(), [](float x, float y) -> float {return x+y;});
  
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
