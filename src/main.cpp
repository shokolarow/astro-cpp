#include "Propagator.h"
#include <iostream>
#include <string>
#include "tools.h"
#include <Eigen/Core>

int main() {

    std::string s1 = "1 27424U 02022A   21065.93350969  .00000101  00000-0  32458-4 0  9991";
    std::string s2 = "2 27424  98.2119   8.7700 0000212 129.0710 304.4809 14.57112576  2108";

    auto koe_date = parse_tle(s1, s2);
    auto koe = koe_date.first;
    auto date = koe_date.second;
    auto csv = koe_to_csv(koe(0), koe(1), koe(2), koe(3), koe(4), koe(5));

    Body sat(100, 0, 50, 2.2, "Satellite");
    std::vector<Body> bodies = {sat};

    double T = 60 * 60;
    double dt = 0.1;

    Propagator sim(csv, bodies, T, dt);

    auto state = sim.propagate();

    return 0;

}