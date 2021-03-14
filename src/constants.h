#pragma once

#include <cmath>

constexpr double G = 6.674e-11;      // Gravitational constant
constexpr double Me = 5.972e24;      // Mass of Earth [kg]
constexpr double Re = 6378.1363e3;        // Earth's radius [m]
constexpr double mu = G * Me;        // Gravitational parameter of Earth
// http://adsabs.harvard.edu/full/1964PASJ...16..263K
constexpr double J2 = 1.0826e-3; // Earth J2 coefficient
constexpr double J3 = -2.546e-6;     // Earth J3 coefficient
constexpr double J4 = -1.649e-6;     // Earth J4 coefficient
constexpr double J5 = -0.210e-6;     // Earth J5 coefficient
constexpr double J6 = 0.646e-6;      // Earth J6 coefficient

constexpr double r2d = 180 / M_PI;
constexpr double d2r = M_PI / 180;
