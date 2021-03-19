/** TRACCC library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "definitions/primitives.hpp"
#include "definitions/algebra.hpp"
#include "edm/spacepoint.hpp"
#include "algorithms/tools/linear_projector.hpp"

#include <gtest/gtest.h>

// This defines the linear projection test
// coming from (0,0,z_v) - phi stays constant
TEST(algorithms, tools_linear_projector_zero)
{

    using namespace traccc;

    linear_projector zero_l_projector;

    spacepoint one{{33., 44., 5.}, {0., 0., 0.}};
    spacepoint two{{66., 88., 15.}, {0., 0., 0.}};
    scalar rho = 200.;

    // Test z @ rho projector
    scalar z_test = zero_l_projector.z_at_rho(one, two, rho);
    scalar theta = std::atan2(two.rho - one.rho, two.global[2] - one.global[2]);
    scalar z_ref = two.global[2] + (rho - two.rho) / std::tan(theta);
    ASSERT_EQ(z_test, z_ref);

    // Test (z,phi) @ rho projector
    auto [z_test_1, phi_test_1] = zero_l_projector.z_phi_at_rho(one, two, rho);
    ASSERT_EQ(z_test_1, z_ref);
    scalar phi_ref_0 = std::atan2(one.global[1], one.global[0]);
    scalar phi_ref_1 = std::atan2(one.global[1], one.global[0]);
    ASSERT_EQ(phi_ref_0, phi_test_1);
    ASSERT_EQ(phi_ref_1, phi_test_1);

    // What's the z at 0 ?
    scalar z_test_at_0 = zero_l_projector.z_at_rho(one, two, 0.);
    scalar z_ref_at_0 = one.global[2] + (-one.rho) / std::tan(theta);
    EXPECT_NEAR(z_test_at_0, z_ref_at_0, 1e-5);

    // Skew transform (to be tested later)
}

// This defines the linear projection test
// - not going through (0,0,z)
// This defines the linear projection test
// coming from (0,0,z_v) - phi stays constant
TEST(algorithms, tools_linear_projector_nonzero)
{

    using namespace traccc;

    linear_projector zero_l_projector;

    spacepoint one{{5., 44., 5.}, {0., 0., 0.}};
    spacepoint two{{66., 88., 15.}, {0., 0., 0.}};
    scalar rho = 200.;

    // Test z @ rho projector
    scalar z_test = zero_l_projector.z_at_rho(one, two, rho);
    scalar theta = std::atan2(two.rho - one.rho, two.global[2] - one.global[2]);
    scalar z_ref = two.global[2] + (rho - two.rho) / std::tan(theta);
    EXPECT_NEAR(z_test, z_ref, 1e-5);

    // Test (z,phi) @ rho projector
    auto [z_test_1, phi_test_1] = zero_l_projector.z_phi_at_rho(one, two, rho);
    EXPECT_NEAR(z_test_1, z_ref, 1e-5);
}

// Google Test can be run manually from the main() function
// or, it can be linked to the gtest_main library for an already
// set-up main() function primed to accept Google Test test cases.
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
