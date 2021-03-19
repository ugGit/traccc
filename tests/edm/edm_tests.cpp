/** TRACCC library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "definitions/primitives.hpp"
#include "definitions/algebra.hpp"
#include "edm/cell.hpp"
#include "edm/cluster.hpp"
#include "edm/measurement.hpp"
#include "edm/spacepoint.hpp"

#include <limits>
#include <cmath>

#include <gtest/gtest.h>

// Test of the cell EDM
TEST(edm, cell)
{

    using namespace traccc;

    cell invalidCell;

    ASSERT_EQ(invalidCell.channel0, std::numeric_limits<channel_id>::max());
    ASSERT_EQ(invalidCell.channel1, std::numeric_limits<channel_id>::max());
    ASSERT_TRUE(std::isnan(invalidCell.activation));
    ASSERT_TRUE(std::isnan(invalidCell.time));

    cell validCell{0, 1, 2., 3.};

    ASSERT_EQ(validCell.channel0, 0);
    ASSERT_EQ(validCell.channel1, 1);
    ASSERT_EQ(validCell.activation, 2.);
    ASSERT_EQ(validCell.time, 3.);
}

// Test of the cluster EDM
TEST(edm, cluster)
{

    using namespace traccc;

    cluster invalidCluser;
    ASSERT_TRUE(invalidCluser.cells.empty());

    cell aCell{0, 1, 2., 3.};
    cluster validCluster{{aCell}};
}

// Test the measurement EDM
TEST(edm, measurement)
{

    using namespace traccc;

    point2 local = {0., 0.};
    variance2 variance = {0., 0.};

    measurement copyM{local, variance};
    measurement moveM{std::move(local), std::move(variance)};
}

// Test the spacepoint EDM
TEST(edm, spaceppoint)
{

    using namespace traccc;

    point3 global = {33., 43., 10.};
    variance3 variance = {0.2, 0.2, 0.02};

    scalar phi = std::atan2(43., 33.);
    scalar rho = std::sqrt(33. * 33. + 43. * 43.);

    spacepoint validSpacePoint(std::move(global), std::move(variance));
}

// Google Test can be run manually from the main() function
// or, it can be linked to the gtest_main library for an already
// set-up main() function primed to accept Google Test test cases.
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}