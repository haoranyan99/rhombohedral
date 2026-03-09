// File: Util/IO/rg_fermiPatch.h
#pragma once

#include <Eigen/Dense>

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <limits>
#include <stdexcept>
#include <array>

#include "Common/DataContainers.h"
#include "Util/IO/rg_io.h"

namespace rgio {

using Vec2 = Eigen::Vector2d;

// ============================================================
// Parameter for cal_chi
// ============================================================
struct cal_chi_param {
    core::GridData qgrid;
    std::array<int,4> q_range;   // {iq_min, iq_max, jq_min, jq_max}

    int Nk = 0;
    int dim = 0;

    double T_K = 0.0;
    double Ef = 0.0;
    double doping = 0.0;
    double filling = 0.0;

    double polar_mu = 0.0;
    bool boundary_periodic = true;
    double eta = 1e-6;
};

// ============================================================
// Output for cal_chi
// ============================================================
struct cal_chi_result {
    core::GridData qgrid;
    double occTot_up = 0.0;
    double occTot_dn = 0.0;
    double doping = 0.0;
    double filling = 0.0;
};

// ============================================================
// Read fermiPatch file into full rectangular GridData
// ============================================================
inline core::GridData
read_fermiPatch_grid_with_bands(const std::string& file_path,
                                rgio::cal_chi_param& param)
{
    std::ifstream fin(file_path);
    if (!fin)
        throw std::runtime_error("Cannot open fermiPatch file: " + file_path);

    int& Nk_       = param.Nk;
    int& dim_      = param.dim;
    double& EF_    = param.Ef;
    double& T_     = param.T_K;
    double& doping_  = param.doping;
    double& filling_ = param.filling;

    int nLayer = -1;
    int Nk = -1;
    double dk_frac = std::numeric_limits<double>::quiet_NaN();
    double EF_eV   = std::numeric_limits<double>::quiet_NaN();
    double T_K     = std::numeric_limits<double>::quiet_NaN();
    double doping  = std::numeric_limits<double>::quiet_NaN();
    double filling = std::numeric_limits<double>::quiet_NaN();
    std::string kmesh_type;

    std::string line;

    // ---------------- read header ----------------
    while (std::getline(fin, line)) {
        if (!is_hash_line(line)) break;

        std::string nohash = trim(line);
        if (!nohash.empty() && nohash[0] == '#')
            nohash = trim(nohash.substr(1));

        parse_header_int(nohash, "nLayer", nLayer);
        parse_kmesh_line(nohash, kmesh_type, Nk, dk_frac);
        parse_header_double(nohash, "EF", EF_eV);
        parse_header_double(nohash, "T", T_K);
        parse_header_double(nohash, "doping", doping);
        parse_header_double(nohash, "filling", filling);
    }

    if (nLayer <= 0)
        throw std::runtime_error("Cannot parse nLayer from header: " + file_path);

    if (Nk < 0)
        throw std::runtime_error("Cannot parse Nk from header: " + file_path);

    const int dim = 2 * nLayer;

    Nk_      = Nk;
    dim_     = dim;
    EF_      = EF_eV;
    T_       = T_K;
    doping_  = doping;
    filling_ = filling;

    // rewind
    fin.clear();
    fin.seekg(0, std::ios::beg);

    // ---------------- allocate full rect grid ----------------
    const int P = 2 * Nk + 1;
    const size_t Ntot = (size_t)P * (size_t)P;

    core::GridData g;
    g.resize(Ntot);
    g.dx = 1.0 / double(P);
    g.dy = 1.0 / double(P);

    for (int iq = -Nk; iq <= Nk; ++iq) {
        for (int jq = -Nk; jq <= Nk; ++jq) {
            const size_t idx =
                (size_t)(iq + Nk) * (size_t)P + (size_t)(jq + Nk);
            g.iq[idx] = iq;
            g.jq[idx] = jq;
        }
    }

    auto& kvecF   = g.add<Vec2>("kvec").v;
    auto& occAvgF = g.add<double>("occ_k_avg").v;
    auto& insideF = g.add<unsigned char>("inside").v;
    auto& evalsF  = g.add<std::vector<double>>("evals").v;
    auto& occbF   = g.add<std::vector<double>>("occ_band").v;

    const double NaN = std::numeric_limits<double>::quiet_NaN();
    for (size_t i = 0; i < Ntot; ++i) {
        kvecF[i]   = Vec2(NaN, NaN);
        occAvgF[i] = NaN;
        insideF[i] = (unsigned char)0;
        evalsF[i].assign(dim, NaN);
        occbF[i].assign(dim, NaN);
    }

    // ---------------- read data rows ----------------
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        if (is_hash_line(line)) continue;

        std::istringstream iss(line);

        long long idx1;
        int iq, jq;
        double kx, ky, occ;

        if (!(iss >> idx1 >> iq >> jq >> kx >> ky >> occ))
            continue;

        if (iq < -Nk || iq > Nk || jq < -Nk || jq > Nk)
            continue;

        const size_t flat =
            (size_t)(iq + Nk) * (size_t)P + (size_t)(jq + Nk);

        insideF[flat] = (unsigned char)1;
        kvecF[flat]   = Vec2(kx, ky);
        occAvgF[flat] = occ;

        for (int b = 0; b < dim; ++b) {
            double Eb;
            if (!(iss >> Eb)) break;
            evalsF[flat][b] = Eb;
        }

        for (int b = 0; b < dim; ++b) {
            double fb;
            if (!(iss >> fb)) break;
            occbF[flat][b] = fb;
        }
    }

    g.assert_consistent();
    return g;
}

} // namespace rgio