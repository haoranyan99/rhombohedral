// File: LinearAlgebra/Constants.h
#pragma once
#include <limits>

namespace la {

constexpr double pi     = 3.141592653589793238462643383279502884;

// Boltzmann constant in eV / K
constexpr double kB_eV = 8.617333262e-5;

constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

} // namespace la
