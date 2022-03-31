/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#include <iostream>

#include "traccc/stdpar/clusterization/test.hpp"

int main(){
  std::cout << "Start STDPAR Example\n";

  traccc::stdpar::execute();

  std::cout << "End STDPAR Example\n";
}
