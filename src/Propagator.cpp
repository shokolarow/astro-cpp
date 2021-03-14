#include "Propagator.h"
#include <cmath>
#include "constants.h"

#define USE_J2
#define NBODY

Body::Body(double mass, double thrust, std::string name)
        : m(mass), T(thrust), name(std::move(name)) {}

Propagator::Propagator(const Eigen::RowVectorXd &init, const std::vector <Body> &bodies, double T, double dt)
        : bodies(bodies), Nbody(bodies.size()), M_Multibody(dof * Nbody), T_final(T), dt(dt) {
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
        double X = u(i*dof);
        double Y = u(i*dof + 1);
        double Z = u(i*dof + 2);

        dudt(i * dof) = u(i * dof + 3);
        dudt(i * dof + 1) = u(i * dof + 4);
        dudt(i * dof + 2) = u(i * dof + 5);

        double R = std::sqrt(X*X + Y*Y + Z*Z);
        double R3 = R * R * R;

        // Gravitational term
        dudt(i * dof + 3) = -mu * X / R3;
        dudt(i * dof + 4) = -mu * Y / R3;
        dudt(i * dof + 5) = -mu * Z / R3;

#ifdef USE_J2
        // J2 effect
        double coeff = 3 * J2 * mu * Re * Re / 2 / std::pow(R,5);
        dudt(i * dof + 3) += coeff * (5 * (Z * Z / R / R) - 1) * X;
        dudt(i * dof + 4) += coeff * (5 * (Z * Z / R / R) - 1) * Y;
        dudt(i * dof + 5) += coeff * (5 * (Z * Z / R / R) - 3) * Z;
#endif

#ifdef NBODY
      for (int j = 0; j < Nbody; j++){
          if (i != j) {
              double Xo = u(j*dof);
              double Yo = u(j*dof + 1);
              double Zo = u(j*dof + 2);
              double Rr = std::sqrt(std::pow(Xo-X,2) +
                                    std::pow(Yo-Y,2) +
                                    std::pow(Zo-Z,2));
              dudt(i * dof + 3) += G * bodies[j].m * (Xo - X) / std::pow(Rr,3);
              dudt(i * dof + 4) += G * bodies[j].m * (Yo - Y) / std::pow(Rr,3);
              dudt(i * dof + 5) += G * bodies[j].m * (Zo - Z) / std::pow(Rr,3);
          }
      }
#endif
    }
}

inline void Propagator::rk4(Eigen::RowVectorXd &u) {
    jacob(u, k1);
    jacob(u + k1 * dt / 2, k2);
    jacob(u + k2 * dt / 2, k3);
    jacob(u + k3 * dt, k4);
    u.noalias() += (k1 + 2 * k2 + 2 * k3 + k4) * dt / 6;
}