#include "Propagator.h"
#include <Eigen/Core>
#include <iostream>
#include <string>
#include "tools.h"

int main() {
    std::string s1 = "1 27424U 02022A   21065.93350969  .00000101  00000-0  32458-4 0  9991";
    std::string s2 = "2 27424  98.2119   8.7700 0000212 129.0710 304.4809 14.57112576  2108";

    auto koe_date = parse_tle(s1, s2);
    auto koe = koe_date.first;
    auto date = koe_date.second;
    auto csv = koe_to_csv(koe(0), koe(1), koe(2), koe(3), koe(4), koe(5));

    std::cout << date << std::endl;
    std::cout << csv << std::endl;
    return 0;
}