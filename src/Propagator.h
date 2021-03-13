#pragma once

#include <Eigen/Dense>
#include <string>
#include <vector>

class Body {
public:
    double m;
    double T;
    std::string name;

    Body(double mass, double thrust, std::string name);
};

class Propagator {
public:
    Eigen::MatrixXd u;
    Eigen::RowVectorXd ut;

    Propagator(const Eigen::RowVectorXd &init, const std::vector <Body> &bodies, double T, double dt);
    Eigen::RowVectorXd propagate();

private:
    static constexpr int dof = 6;
    const int Nbody;
    const int M_Multibody;
    const double T_final;
    double dt;
    int N;

    Eigen::RowVectorXd dudt;
    Eigen::RowVectorXd k1;
    Eigen::RowVectorXd k2;
    Eigen::RowVectorXd k3;
    Eigen::RowVectorXd k4;

    void jacob(const Eigen::RowVectorXd &u, Eigen::RowVectorXd &dudt);
    void rk4(Eigen::RowVectorXd &u);

};
