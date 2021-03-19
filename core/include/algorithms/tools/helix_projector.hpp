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

namespace traccc
{

  /// A purely linear projector
  struct helix_projector
  {

    // Helix radius and vertex assumption
    scalar R = 0.;
    scalar tz = 0.;

    // The approximate helixal projector given pT and B
    helix_projector(scalar vtxz, scalar pT, scalar B)
        : R(pT / B), tz(vtxz)
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

    std::array<scalar, 2> z_phi_at_rho(const spacepoint &one,
                                       const spacepoint &two,
                                       scalar rho) const
    {

      scalar z = one.global[2] - (two.global[2] - one.global[2]) / (two.rho - one.rho) * (one.rho - rho);

      //scalar zphi = one.phi - (two.global[2]-one.global[2])/(two.rho-one.rho) * (one.rho - rho)

      return {z, one.phi};
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
      scalar phi = sp.phi - 0.5 * sp.rho / R;
      // The regular component
      scalar dz = sp.global[2] - tz;
      //scalar zrho = dz/hit.rho;
      scalar theta = std::atan2(sp.rho, dz);
      scalar eta = -std::log(tan(0.5 * theta));

      //scalar eta = log(zrho+sqrt(1.+zrho*zrho));
      //eta *= (dz < 0.) ? 1 : -1;
      //
      return {phi, eta};
    }
  };

} // end of namespace