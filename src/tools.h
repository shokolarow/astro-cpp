#pragma once

#include <cmath>
#include <vector>
#include <iostream>
#include <iterator>
#include <ctime>
#include <utility>
#include <string>

#include <Eigen/Dense>

#include "constants.h"
#include "tools.h"


struct datetime{
    double year,month,day,hour,minute,second;
    datetime(double yr, double mn, double d, double h, double m, double s);
};

std::ostream& operator<<(std::ostream &o, const datetime& date);

Eigen::RowVectorXd koe_to_csv(double a, double e, double i, double om, double w, double nu);
double MA2v(double MA, double e);
std::pair<Eigen::Matrix<double,1,6>, datetime> parse_tle(const std::string& s1, const std::string& s2);

constexpr double compute_rho(double alt){
    using namespace EARTH;
    // http://www.braeunig.us/space/atmmodel.htm
    alt /= 1E3;
    double rho = 0;
    if ((alt >= 200) && (alt < 300)){
        rho = std::exp(rho_coeff[0][0]*std::pow(alt,4) +
                       rho_coeff[0][1]*std::pow(alt,3) +
                       rho_coeff[0][2]*std::pow(alt,2) +
                       rho_coeff[0][3]*alt + rho_coeff[0][4]);
    }
    else if ((alt >= 300) && (alt < 500)){
        rho = std::exp(rho_coeff[1][0]*std::pow(alt,4) +
                       rho_coeff[1][1]*std::pow(alt,3) +
                       rho_coeff[1][2]*std::pow(alt,2) +
                       rho_coeff[1][3]*alt + rho_coeff[1][4]);
    }
    else if ((alt >= 500) && (alt < 750)){
        rho = std::exp(rho_coeff[2][0]*std::pow(alt,4) +
                       rho_coeff[2][1]*std::pow(alt,3) +
                       rho_coeff[2][2]*std::pow(alt,2) +
                       rho_coeff[2][3]*alt + rho_coeff[2][4]);
    }
    else if ((alt >= 750) && (alt < 1000)){
        rho = std::exp(rho_coeff[3][0]*std::pow(alt,4) +
                       rho_coeff[3][1]*std::pow(alt,3) +
                       rho_coeff[3][2]*std::pow(alt,2) +
                       rho_coeff[3][3]*alt + rho_coeff[3][4]);
    }
    return rho;
}

