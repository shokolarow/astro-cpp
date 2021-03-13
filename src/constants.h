#pragma once

#include <cmath>

constexpr double G = 6.674e-11;      // Gravitational constant
constexpr double Me = 5.972e24;      // Mass of Earth [kg]
constexpr double mu = G * Me;        // Gravitational parameter of Earth
constexpr double J2 = 1.08262668e-3; // Earth J2 coefficient
constexpr double r2d = 180 / M_PI;
constexpr double d2r = M_PI / 180;
