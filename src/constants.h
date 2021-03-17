#pragma once

#include <cmath>
#include <Eigen/Dense>

constexpr double G = 6.674e-11;      // Gravitational constant
constexpr double r2d = 180 / M_PI;
constexpr double d2r = M_PI / 180;

namespace EARTH{

    constexpr double Me = 5.972e24;      // Mass of Earth [kg]
    constexpr double Re = 6378.1363e3;        // Earth's radius [m]
    constexpr double mu = G * Me;        // Gravitational parameter of Earth
    // http://adsabs.harvard.edu/full/1964PASJ...16..263K
    constexpr double J2 = 1.0826e-3; // Earth J2 coefficient
    constexpr double J3 = -2.546e-6;     // Earth J3 coefficient
    constexpr double J4 = -1.649e-6;     // Earth J4 coefficient
    constexpr double J5 = -0.210e-6;     // Earth J5 coefficient
    constexpr double J6 = 0.646e-6;      // Earth J6 coefficient
    constexpr double EARTH_ROT_SPEED = 2 * M_PI / 24 / 60 / 60;
    inline const Eigen::Vector3d EARTH_ROT_VECTOR(0,0,EARTH_ROT_SPEED);

    // http://www.braeunig.us/space/atmmodel.htm
    constexpr double rho_coeff[4][5] = {{1.199282E-09, -1.451051E-06, 6.910474E-04, -0.1736220, -5.321644},
                                        {1.140564E-10, -2.130756E-07, 1.570762E-04, -0.07029296, -12.89844},
                                        {8.105631E-12, -2.358417E-09, -2.635110E-06, -0.01562608, -20.02246},
                                        {-3.701195E-12, -8.608611E-09, 5.118829E-05, -0.06600998, -6.137674}};
}




