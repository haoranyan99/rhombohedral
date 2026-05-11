// File: MeasureEngines/cal_susceptibility.h
#pragma once

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "Common/DataContainers.h"
#include "LinearAlgebra/MathFunctions.h"
#include "Util/MPI/mpi.h"

namespace rgio {

using Vec2 = Eigen::Vector2d;
using cd   = std::complex<double>;

struct cal_chi_param {
    int iq = 0;
    int jq = 0;

    int Nk  = 0;
    int dim = 0;

    double T_K = 0.0;
    double Ef  = 0.0;

    double doping  = 0.0;
    double filling = 0.0;

    double polar_mu = 0.0;
    double eta = 1e-6;

    double lattice_a = 2.46;

    // Final normalization:
    // chi = sum_k(...) * area_density.
    // If <=0, kernel uses fgrid.dx * fgrid.dy.
    double area_density = 0.0;

    bool boundary_periodic = false;
    bool use_form_factor = true;

    std::string mesh_type;
};

struct cal_chi_result {
    core::GridData qgrid;

    double occTot_up = 0.0;
    double occTot_dn = 0.0;

    double doping  = 0.0;
    double filling = 0.0;
};

inline core::GridData read_fermiPatch_bin(
    const std::string& file_path,
    cal_chi_param& param
) {
    std::ifstream ifs(file_path, std::ios::binary);
    if (!ifs) {
        throw std::runtime_error("read_fermiPatch_bin: cannot open file: " + file_path);
    }

    int32_t magic = 0;
    int32_t version = 0;
    int32_t NkTot_i32 = 0;
    int32_t dim_i32 = 0;

    double EF = 0.0;
    double T_K = 0.0;
    double doping = 0.0;
    double filling = 0.0;
    double dx = 0.0;
    double dy = 0.0;

    char mesh_type_buf[32] = {};

    ifs.read(reinterpret_cast<char*>(&magic), sizeof(magic));
    ifs.read(reinterpret_cast<char*>(&version), sizeof(version));
    ifs.read(reinterpret_cast<char*>(&NkTot_i32), sizeof(NkTot_i32));
    ifs.read(reinterpret_cast<char*>(&dim_i32), sizeof(dim_i32));

    ifs.read(reinterpret_cast<char*>(&EF), sizeof(EF));
    ifs.read(reinterpret_cast<char*>(&T_K), sizeof(T_K));
    ifs.read(reinterpret_cast<char*>(&doping), sizeof(doping));
    ifs.read(reinterpret_cast<char*>(&filling), sizeof(filling));

    ifs.read(reinterpret_cast<char*>(&dx), sizeof(dx));
    ifs.read(reinterpret_cast<char*>(&dy), sizeof(dy));
    ifs.read(reinterpret_cast<char*>(mesh_type_buf), sizeof(mesh_type_buf));

    if (!ifs) {
        throw std::runtime_error("read_fermiPatch_bin: failed to read header: " + file_path);
    }

    if (magic != 20260510) {
        throw std::runtime_error("read_fermiPatch_bin: wrong magic number");
    }

    if (version != 3) {
        throw std::runtime_error("read_fermiPatch_bin: unsupported version");
    }

    if (NkTot_i32 <= 0 || dim_i32 <= 0) {
        throw std::runtime_error("read_fermiPatch_bin: invalid NkTot or dim");
    }

    const size_t NkTot = static_cast<size_t>(NkTot_i32);
    const int dim = static_cast<int>(dim_i32);

    core::GridData g;
    g.resize(NkTot);
    g.dx = dx;
    g.dy = dy;

    auto& kvecF   = g.add<Vec2>("kvec").v;
    auto& occAvgF = g.add<double>("occ_k_avg").v;
    auto& evalsF  = g.add<std::vector<double>>("evals").v;
    auto& occbF   = g.add<std::vector<double>>("occ_band").v;
    auto& evecReF = g.add<std::vector<double>>("evec_re").v;
    auto& evecImF = g.add<std::vector<double>>("evec_im").v;

    for (size_t ik = 0; ik < NkTot; ++ik) {
        int32_t iq = 0;
        int32_t jq = 0;

        double kx = 0.0;
        double ky = 0.0;
        double occ_avg = 0.0;

        ifs.read(reinterpret_cast<char*>(&iq), sizeof(iq));
        ifs.read(reinterpret_cast<char*>(&jq), sizeof(jq));
        ifs.read(reinterpret_cast<char*>(&kx), sizeof(kx));
        ifs.read(reinterpret_cast<char*>(&ky), sizeof(ky));
        ifs.read(reinterpret_cast<char*>(&occ_avg), sizeof(occ_avg));

        g.iq[ik] = static_cast<int>(iq);
        g.jq[ik] = static_cast<int>(jq);
        kvecF[ik] = Vec2(kx, ky);
        occAvgF[ik] = occ_avg;

        evalsF[ik].resize(static_cast<size_t>(dim));
        occbF[ik].resize(static_cast<size_t>(dim));
        evecReF[ik].resize(static_cast<size_t>(dim) * static_cast<size_t>(dim));
        evecImF[ik].resize(static_cast<size_t>(dim) * static_cast<size_t>(dim));

        ifs.read(
            reinterpret_cast<char*>(evalsF[ik].data()),
            sizeof(double) * static_cast<size_t>(dim)
        );

        ifs.read(
            reinterpret_cast<char*>(occbF[ik].data()),
            sizeof(double) * static_cast<size_t>(dim)
        );

        ifs.read(
            reinterpret_cast<char*>(evecReF[ik].data()),
            sizeof(double) * static_cast<size_t>(dim) * static_cast<size_t>(dim)
        );

        ifs.read(
            reinterpret_cast<char*>(evecImF[ik].data()),
            sizeof(double) * static_cast<size_t>(dim) * static_cast<size_t>(dim)
        );

        if (!ifs) {
            throw std::runtime_error("read_fermiPatch_bin: failed reading body");
        }
    }

    int Nk_infer = 0;
    for (size_t ik = 0; ik < NkTot; ++ik) {
        Nk_infer = std::max(Nk_infer, std::abs(g.iq[ik]));
        Nk_infer = std::max(Nk_infer, std::abs(g.jq[ik]));
        Nk_infer = std::max(Nk_infer, std::abs(g.iq[ik] + g.jq[ik]));
    }

    param.Nk = Nk_infer;
    param.dim = dim;
    param.Ef = EF;
    param.T_K = T_K;
    param.doping = doping;
    param.filling = filling;
    param.mesh_type = std::string(mesh_type_buf);
    param.area_density = dx * dy;

    g.assert_consistent();
    return g;
}

} // namespace rgio


namespace rg {

using Vec2 = Eigen::Vector2d;
using cd   = std::complex<double>;

inline void infer_grid_range_(
    const core::GridData& g,
    int& iq_min,
    int& iq_max,
    int& jq_min,
    int& jq_max
) {
    iq_min =  1000000000;
    iq_max = -1000000000;
    jq_min =  1000000000;
    jq_max = -1000000000;

    for (size_t i = 0; i < g.size(); ++i) {
        iq_min = std::min(iq_min, g.iq[i]);
        iq_max = std::max(iq_max, g.iq[i]);
        jq_min = std::min(jq_min, g.jq[i]);
        jq_max = std::max(jq_max, g.jq[i]);
    }
}

inline bool find_shifted_index_(
    const core::GridData& g,
    const std::string& mesh_type,
    int iq,
    int jq,
    bool periodic,
    size_t& idx
) {
    if (!periodic || mesh_type == "hex") {
        return g.ij_to_idx(iq, jq, idx);
    }

    int iq_min = 0;
    int iq_max = 0;
    int jq_min = 0;
    int jq_max = 0;

    infer_grid_range_(g, iq_min, iq_max, jq_min, jq_max);

    const int N1 = iq_max - iq_min + 1;
    const int N2 = jq_max - jq_min + 1;

    const int wiq = iq_min + la::mod_pos_int(iq - iq_min, N1);
    const int wjq = jq_min + la::mod_pos_int(jq - jq_min, N2);

    return g.ij_to_idx(wiq, wjq, idx);
}

inline bool infer_dq_vectors_(
    const core::GridData& g,
    Vec2& dq1,
    Vec2& dq2
) {
    const auto& kvec = g.get<Vec2>("kvec").v;

    size_t i00 = 0;
    size_t i10 = 0;
    size_t i01 = 0;

    if (
        g.ij_to_idx(0, 0, i00)
     && g.ij_to_idx(1, 0, i10)
     && g.ij_to_idx(0, 1, i01)
    ) {
        dq1 = kvec[i10] - kvec[i00];
        dq2 = kvec[i01] - kvec[i00];
        return true;
    }

    return false;
}

inline double form_factor_(
    const std::vector<double>& u1_re,
    const std::vector<double>& u1_im,
    const std::vector<double>& u2_re,
    const std::vector<double>& u2_im,
    int dim,
    int band1,
    int band2
) {
    cd ov(0.0, 0.0);

    for (int a = 0; a < dim; ++a) {
        const size_t p1 =
            static_cast<size_t>(a) * static_cast<size_t>(dim)
          + static_cast<size_t>(band1);

        const size_t p2 =
            static_cast<size_t>(a) * static_cast<size_t>(dim)
          + static_cast<size_t>(band2);

        const cd u1(u1_re[p1], u1_im[p1]);
        const cd u2(u2_re[p2], u2_im[p2]);

        ov += std::conj(u1) * u2;
    }

    return std::norm(ov);
}

inline rgio::cal_chi_result cal_chi_grid_from_fermiPatch(
    const core::GridData& fgrid,
    rgio::cal_chi_param param
) {
    int rank = 0;
    int nprocs = 1;
    rgmpi::rank_size(rank, nprocs);

    fgrid.assert_consistent();

    const size_t NkTot = fgrid.size();
    if (NkTot == 0) {
        rgmpi::abort_all("cal_chi_grid_from_fermiPatch: empty fgrid");
    }

    const int dim = param.dim;
    if (dim <= 0) {
        rgmpi::abort_all("cal_chi_grid_from_fermiPatch: dim <= 0");
    }

    if (!(param.eta > 0.0)) {
        rgmpi::abort_all("cal_chi_grid_from_fermiPatch: eta must be > 0");
    }

    const auto& evals  = fgrid.get<std::vector<double>>("evals").v;
    const auto& occ    = fgrid.get<std::vector<double>>("occ_band").v;
    const auto& evecRe = fgrid.get<std::vector<double>>("evec_re").v;
    const auto& evecIm = fgrid.get<std::vector<double>>("evec_im").v;

    rgio::cal_chi_result result;
    result.doping = param.doping;
    result.filling = param.filling;

    core::GridData& qgrid = result.qgrid;
    qgrid.resize(1);
    qgrid.dx = fgrid.dx;
    qgrid.dy = fgrid.dy;

    qgrid.iq[0] = param.iq;
    qgrid.jq[0] = param.jq;

    auto& chiF   = qgrid.add<cd>("chi").v;
    auto& nPairF = qgrid.add<long long>("nKpair").v;
    auto& qvecF  = qgrid.add<Vec2>("qvec").v;
    auto& qxF    = qgrid.add<double>("qx").v;
    auto& qyF    = qgrid.add<double>("qy").v;

    Vec2 dq1(0.0, 0.0);
    Vec2 dq2(0.0, 0.0);

    if (infer_dq_vectors_(fgrid, dq1, dq2)) {
        const Vec2 qvec =
            static_cast<double>(param.iq) * dq1
          + static_cast<double>(param.jq) * dq2;

        qvecF[0] = qvec;
        qxF[0] = qvec.x();
        qyF[0] = qvec.y();
    } else {
        const double NaN = std::numeric_limits<double>::quiet_NaN();
        qvecF[0] = Vec2(NaN, NaN);
        qxF[0] = NaN;
        qyF[0] = NaN;
    }

    chiF[0] = cd(0.0, 0.0);
    nPairF[0] = 0;

    const bool use_polar = std::abs(param.polar_mu) > 1e-15;

    std::vector<double> occ_up;
    std::vector<double> occ_dn;

    if (use_polar) {
        occ_up.assign(NkTot * static_cast<size_t>(dim), 0.0);
        occ_dn.assign(NkTot * static_cast<size_t>(dim), 0.0);

        result.occTot_up = 0.0;
        result.occTot_dn = 0.0;

        for (size_t ik = 0; ik < NkTot; ++ik) {
            for (int b = 0; b < dim; ++b) {
                const double E = evals[ik][static_cast<size_t>(b)];

                const double E_up = E - 0.5 * param.polar_mu;
                const double E_dn = E + 0.5 * param.polar_mu;

                const double fu = la::fermi(E_up, param.Ef, param.T_K);
                const double fd = la::fermi(E_dn, param.Ef, param.T_K);

                const size_t p =
                    ik * static_cast<size_t>(dim)
                  + static_cast<size_t>(b);

                occ_up[p] = fu;
                occ_dn[p] = fd;

                result.occTot_up += fu;
                result.occTot_dn += fd;
            }
        }

        result.filling =
            (result.occTot_up + result.occTot_dn)
          / (2.0 * static_cast<double>(NkTot) * static_cast<double>(dim));

        result.doping =
            la::filling_to_doping(result.filling, param.lattice_a);
    }

    const auto [k_start, k_end] =
        rgmpi::block_1d_int(static_cast<int>(NkTot), rank, nprocs);

    cd chi_local(0.0, 0.0);
    long long nPair_local = 0;

    for (int ii = k_start; ii < k_end; ++ii) {
        const size_t ik1 = static_cast<size_t>(ii);

        const int iq1 = fgrid.iq[ik1];
        const int jq1 = fgrid.jq[ik1];

        size_t ik2 = 0;
        if (!find_shifted_index_(
                fgrid,
                param.mesh_type,
                iq1 + param.iq,
                jq1 + param.jq,
                param.boundary_periodic,
                ik2
            )) {
            continue;
        }

        for (int b = 0; b < dim; ++b) {
            const double Eb = evals[ik1][static_cast<size_t>(b)];

            for (int m = 0; m < dim; ++m) {
                const double Em = evals[ik2][static_cast<size_t>(m)];

                double FF = 1.0;
                if (param.use_form_factor) {
                    FF = form_factor_(
                        evecRe[ik1],
                        evecIm[ik1],
                        evecRe[ik2],
                        evecIm[ik2],
                        dim,
                        b,
                        m
                    );
                }

                if (!use_polar) {
                    const double fb = occ[ik1][static_cast<size_t>(b)];
                    const double fm = occ[ik2][static_cast<size_t>(m)];

                    chi_local +=
                        FF * (fm - fb) / cd(Eb - Em, param.eta);

                    ++nPair_local;
                } else {
                    const size_t p1 =
                        ik1 * static_cast<size_t>(dim)
                      + static_cast<size_t>(b);

                    const size_t p2 =
                        ik2 * static_cast<size_t>(dim)
                      + static_cast<size_t>(m);

                    const double Ebu = Eb - 0.5 * param.polar_mu;
                    const double Ebd = Eb + 0.5 * param.polar_mu;
                    const double Emu = Em - 0.5 * param.polar_mu;
                    const double Emd = Em + 0.5 * param.polar_mu;

                    const double fbu = occ_up[p1];
                    const double fbd = occ_dn[p1];
                    const double fmu = occ_up[p2];
                    const double fmd = occ_dn[p2];

                    chi_local +=
                        0.25 * FF * (fmu - fbu) / cd(Ebu - Emu, param.eta);

                    chi_local +=
                        0.25 * FF * (fmd - fbd) / cd(Ebd - Emd, param.eta);

                    chi_local +=
                        0.25 * FF * (fmu - fbd) / cd(Ebd - Emu, param.eta);

                    chi_local +=
                        0.25 * FF * (fmd - fbu) / cd(Ebu - Emd, param.eta);

                    nPair_local += 4;
                }
            }
        }
    }

    double chi_re_sum = 0.0;
    double chi_im_sum = 0.0;
    long long nPair_sum = 0;

    rgmpi::allreduce_sum(chi_local.real(), chi_re_sum);
    rgmpi::allreduce_sum(chi_local.imag(), chi_im_sum);
    rgmpi::allreduce_sum(nPair_local, nPair_sum);

    double area_density = param.area_density;
    if (!(area_density > 0.0)) {
        area_density = fgrid.dx * fgrid.dy;
    }
    if (!(area_density > 0.0)) {
        area_density = 1.0;
    }

    chiF[0] = cd(chi_re_sum, chi_im_sum) * area_density;
    nPairF[0] = nPair_sum;

    qgrid.assert_consistent();
    return result;
}

inline rgio::cal_chi_result cal_chi_grid_from_fermiPatch(
    const std::string& file_path,
    rgio::cal_chi_param& param
) {
    core::GridData fgrid =
        rgio::read_fermiPatch_bin(file_path, param);

    return cal_chi_grid_from_fermiPatch(fgrid, param);
}

} // namespace rg