/** TRACCC library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "edm/spacepoint.hpp"
#include "definitions/primitives.hpp"
#include "definitions/algebra.hpp"
#include "algorithms/tools/quadratic_equation.hpp"

namespace traccc
{

  /// A purely linear projector
  struct linear_projector
  {

    /// Helix radius and vertex assumption
    scalar tz = 0.;

    linear_projector() = default;

    /// The linear projector, constructor for a dedicated
    /// vertex @param vtx, it ignores momentum and b field
    ///
    linear_projector(scalar vtxz, scalar /*pT*/, scalar /*B*/) : tz(vtxz)
    {
    }

    /// Helper method to estimate z at a given rho
    ///
    /// @param one first space point
    /// @param two second space point
    /// @param rho target rho
    ///
    /// @return z at given rho
    scalar z_at_rho(const spacepoint &one, const spacepoint &two, scalar rho) const
    {
      return one.global[2] - (two.global[2] - one.global[2]) / (two.rho - one.rho) * (one.rho - rho);
    }

    /// Helper method to estimate z at a given rho
    ///
    /// @param one first space point
    /// @param two second space point
    /// @param rho target rho
    ///
    /// @return z,phi at given rho
    std::array<scalar, 2> z_phi_at_rho(const spacepoint &one,
                                       const spacepoint &two,
                                       scalar rho) const
    {

      scalar z = one.global[2] - (two.global[2] - one.global[2]) / (two.rho - one.rho) * (one.rho - rho);

      // Swap coorinates x/y for numerical stability
      bool swap_x_y = std::abs(one.global[0] - two.global[0]) < 1e-3;

      unsigned int _x = swap_x_y ? 1 : 0;
      unsigned int _y = swap_x_y ? 0 : 1;
      scalar k = (one.global[_y] - two.global[_y]) / (one.global[_x] - two.global[_x]);
      scalar d = two.global[_y] - k * two.global[_x];
      scalar phi = std::numeric_limits<scalar>::quiet_NaN();

      quadratic_equation qe = {(1 + k * k), 2 * k * d, d * d - rho * rho};
      auto qe_solution = qe();

      if (std::get<0>(qe_solution) > 0)
      {
        vector3 rd = two.global - one.global;
        
        std::array<point3, 2> candidates;
        auto u01 = std::get<1>(qe_solution);
        std::array<scalar, 2> t01 = {0., 0.};

        candidates[0][_x] = u01[0];
        candidates[0][_y] = k * u01[0] + d;
        t01[0] = (candidates[0][_x] - one.global[_x]) / rd[_x];
        candidates[0][2] = one.global[2] + t01[0] * rd[2];

        candidates[1][_x] = u01[1];
        candidates[1][_y] = k * u01[1] + d;
        t01[1] = (candidates[1][_x] - one.global[_x]) / rd[_x];
        candidates[1][2] = one.global[2] + t01[1] * rd[2];

        // Chose the index, take the smaller positive one
        int cindex = (t01[0] < t01[1] and t01[0] > 0.) ? 0
                                                       : (t01[0] < 0. and t01[1] > 0. ? 1 : 0);
        if (t01[0] > 0. or t01[1] > 0.)
        {
          phi = getter::phi(candidates[cindex]);
        }
      }
      return {z, phi};
    }

    /// Helper method for the skew transform
    ///
    /// @param sp is the spacepoint to be transformed
    /// @param z is the z vertex estimate
    ///
    /// @return z,phi at given rho
    std::array<scalar, 2> skew_transform(const spacepoint &sp, scalar z) const
    {
      // The ciruclar component
      scalar phi = sp.phi;
      // The regular component
      scalar dz = sp.global[2] - tz;
      scalar zrho = dz / sp.rho;
      scalar eta = std::log(zrho + std::sqrt(1. + zrho * zrho));
      eta *= (dz < 0.) ? 1 : -1;

      return {phi, eta};
    }
  };

} // end of namespace