/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <functional>

#include "tests/cca_test.hpp"
#include "traccc/stdpar/clusterization/clusterization_algorithm.hpp"

namespace {
vecmem::host_memory_resource resource;
traccc::stdpar::clusterization_algorithm cca(resource);

cca_function_t f = [](const traccc::cell_container_types::host &data) {
    std::map<traccc::geometry_id, vecmem::vector<traccc::measurement>> result;

    traccc::measurement_container_types::host mss = cca(data);

    for (std::size_t i = 0; i < mss.size(); ++i) {
        std::vector<traccc::measurement, std::pmr::polymorphic_allocator<traccc::measurement>> msv;

        for (std::size_t j = 0; j < mss.at(i).items.size(); ++j) {
            msv.push_back(mss.at(i).items.at(j));
        }

        result.emplace(mss.at(i).header.module, std::move(msv));
    }

    return result;
};
}  // namespace

TEST_P(ConnectedComponentAnalysisTests, Run) {
    test_connected_component_analysis(GetParam());
}

INSTANTIATE_TEST_SUITE_P(
    SimplifiedSvAlgorithmStdPar, ConnectedComponentAnalysisTests,
    ::testing::Combine(
        ::testing::Values(f),
        ::testing::ValuesIn(ConnectedComponentAnalysisTests::get_test_files())),
    ConnectedComponentAnalysisTests::get_test_name);
