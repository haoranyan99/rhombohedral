// LinearAlgebra/MathFunctions.h
#pragma once
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
#include "Constants.h"

namespace la {

inline int mod_pos_int(int a, int m) {
    int r = a % m;
    return (r < 0) ? (r + m) : r;
}

inline double fermi(double E, double mu, double T_K)
{
    if (T_K <= 0.0) {
        return (E < mu) ? 1.0 : 0.0;
    }

    const double beta = 1.0 / (la::kB_eV * T_K);
    const double x = (E - mu) * beta;

    if (x >  50.0) return 0.0;
    if (x < -50.0) return 1.0;

    return 1.0 / (std::exp(x) + 1.0);
}

inline std::vector<double> linspace(double a, double b, int n) {
    if (n < 1) throw std::runtime_error("linspace: n must be >= 1");
    std::vector<double> x;
    x.reserve((size_t)n);
    if (n == 1) { x.push_back(a); return x; }
    const double dx = (b - a) / double(n - 1);
    for (int i = 0; i < n; ++i) x.push_back(a + dx * double(i));
    return x;
}

inline double doping_to_filling(double doping, double a = 2.46)
{
    // Area of hexagonal unit cell: sqrt(3)/2 * a^2  (Angstrom^2)
    // 1 Angstrom^2 = 1e-16 cm^2
    // doping unit = 1e12 cm^-2
    //
    // delta_filling = doping * A_cell * 1e12 * 1e-16
    //               = doping * (sqrt(3)/2 * a^2) * 1e-4
    const double unit_conv = std::sqrt(3.0) * 0.5 * a * a * 1.0e-4;

    return 0.5 + doping * unit_conv;
}

inline double filling_to_doping(double filling, double a = 2.46)
{
    const double unit_conv = std::sqrt(3.0) * 0.5 * a * a * 1.0e-4;

    return (filling - 0.5) / unit_conv;
}


} // namespace la

