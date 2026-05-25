// File: LinearAlgebra/Constants.h
#pragma once
#include <limits>

namespace la {

constexpr double pi     = 3.141592653589793238462643383279502884;

// Boltzmann constant in eV / K
constexpr double kB_eV = 8.617333262e-5;

// Bohr magneton in eV / Tesla
constexpr double muB_eV_per_T = 5.7883818060e-5;

// 2 m_e / hbar^2 in units of 1 / (eV Angstrom^2).
// Converts the eV*Angstrom^2 orbital-moment kernel in Eq. (17) to mu_B.
constexpr double orbital_moment_muB_prefactor = 0.262464239045862;

constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

} // namespace la
