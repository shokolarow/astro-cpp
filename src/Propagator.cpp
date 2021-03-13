#include "Propagator.h"
#include <cmath>
#include "constants.h"

Body::Body(double mass, double thrust, std::string name)
        : m(mass), T(thrust), name(std::move(name)) {}

Propagator::Propagator(const Eigen::RowVectorXd &init, const std::vector <Body> &bodies, double T, double dt)
        : Nbody(bodies.size()), M_Multibody(dof * Nbody), T_final(T), dt(dt) {
    N = static_cast<int>(T / dt);
    u = Eigen::MatrixXd::Zero(N, M_Multibody);
    ut = Eigen::RowVectorXd::Zero(M_Multibody);
    dudt = Eigen::RowVectorXd::Zero(M_Multibody);
    k1 = Eigen::RowVectorXd::Zero(M_Multibody);
    k2 = Eigen::RowVectorXd::Zero(M_Multibody);
    k3 = Eigen::RowVectorXd::Zero(M_Multibody);
    k4 = Eigen::RowVectorXd::Zero(M_Multibody);

    u.row(0) = init;
    ut = init;
}

Eigen::RowVectorXd Propagator::propagate() {
    for (int i = 0; i < N - 1; i++) {
        rk4(ut);
        u.row(i + 1) = ut;
    }

    return ut;
}

void Propagator::jacob(const Eigen::RowVectorXd &u, Eigen::RowVectorXd &dudt) {
    for (int i = 0; i < Nbody; i++) {
        dudt(i * dof) = u(i * dof + 3);
        dudt(i * dof + 1) = u(i * dof + 4);
        dudt(i * dof + 2) = u(i * dof + 5);

        double R = std::sqrt(u(i * dof) * u(i * dof) +
                             u(i * dof + 1) * u(i * dof + 1) +
                             u(i * dof + 2) * u(i * dof + 2));
        double R3 = R * R * R;

        dudt(i * dof + 3) = -mu * u(i * dof) / R3;
        dudt(i * dof + 4) = -mu * u(i * dof + 1) / R3;
        dudt(i * dof + 5) = -mu * u(i * dof + 2) / R3;
    }
}

void Propagator::rk4(Eigen::RowVectorXd &u) {
    jacob(u, k1);
    jacob(u + k1 * dt / 2, k2);
    jacob(u + k2 * dt / 2, k3);
    jacob(u + k3 * dt, k4);
    u.noalias() += (k1 + 2 * k2 + 2 * k3 + k4) * dt / 6;
}