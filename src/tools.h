#pragma once

#include <string>
#include <Eigen/Dense>
#include <utility>

struct datetime{
    double year,month,day,hour,minute,second;
    datetime(double yr, double mn, double d, double h, double m, double s);
};

std::ostream& operator <<(std::ostream &o, const datetime& date);

Eigen::RowVectorXd koe_to_csv(double a, double e, double i, double om, double w, double nu);
double MA2v(double MA, double e);
std::pair<Eigen::Matrix<double,1,6>, datetime> parse_tle(const std::string& s1, const std::string& s2);
