/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <functional>

#include "tests/cca_test.hpp"
#include "traccc/stdpar/clusterization/component_connection_ssv.hpp"

namespace {

traccc::stdpar::component_connection_ssv ca;

cca_function_t f = [](const traccc::cell_container_types::host &data) {
    std::map<traccc::geometry_id, vecmem::vector<traccc::measurement>> result;

    auto measurements = ca(data);
    for (std::size_t i = 0; i < measurements.size(); i++) {
        result[measurements.at(i).header.module] = measurements.at(i).items;
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
