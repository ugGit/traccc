/** TRACCC library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <tuple>
#include <array>

namespace traccc
{
    /** Struct to solve a quadratic equation of type p[0] * x^2 + p[1] * x + p[2] = 0 
     *
     * Note: this is algorithmic identical to detray/quadratic_equation
     */
    struct quadratic_equation
    {
        std::array<scalar, 3> _params = {0., 0., 0.};

        /** Solve the quadratic equation 
         **/
        std::tuple<int, std::array<scalar, 2> >
        operator()() const
        {
            scalar discriminant = _params[1] * _params[1] - 4 * _params[0] * _params[2];
            if (discriminant < 0.)
            {
                return {0, {std::numeric_limits<scalar>::infinity(), std::numeric_limits<scalar>::infinity()}};
            }
            else
            {
                int solutions = (discriminant == 0.) ? 1 : 2;
                double q = -0.5 * (_params[1] + (_params[1] > 0 ? std::sqrt(discriminant)
                                                                : -std::sqrt(discriminant)));
                scalar first = q / _params[0];
                scalar second = _params[2] / q;
                std::array<scalar, 2> poles = {first, second};
                std::sort(poles.begin(), poles.end());
                return {solutions, poles};
            }
        }
    };

} // namespace traccc
