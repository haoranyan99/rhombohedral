// File: MeasureEngines/cal_susceptibility.h
#pragma once
#include <Eigen/Dense>
#include <complex>
#include <limits>
#include <vector>
#include <cmath>
#include <algorithm>
#include <array>
#include <filesystem>
#include <string>
#include <iostream>

#ifdef USE_MPI
  #include <mpi.h>
#endif

#include "Common/DataContainers.h"
#include "LinearAlgebra/Constants.h"
#include "LinearAlgebra/MathFunctions.h"
#include "Util/MPI/mpi.h"
#include "Util/IO/rg_fermiPatch.h"
#include "Util/IO/rg_io.h"

namespace rg {

using Vec2 = Eigen::Vector2d;
using cd   = std::complex<double>;

inline int wrap_q_index(int i, int Nk) {
    const int P = 2 * Nk;
    return ((i + Nk) % P + P) % P - Nk;
}

// ============================================================
// Core kernel: compute chi from an already-loaded fgrid.
// - If polar_mu==0: use occ_band from fgrid.
// - If polar_mu!=0:
//     if fgrid already has occ_up/occ_dn (injected by path overload), use them;
//     else fallback: compute occ_up/occ_dn from SAME fgrid by Fermi function.
// ============================================================
inline rgio::cal_chi_result cal_chi_grid_from_fermiPatch(
    const core::GridData& fgrid_in,
    const rgio::cal_chi_param& param)
{
    int rank = 0, nprocs = 1;
    rgmpi::rank_size(rank, nprocs);

    rgio::cal_chi_result result;
    double& occTot_up = result.occTot_up;
    double& occTot_dn = result.occTot_dn;
    double& doping    = result.doping;
    double& filling   = result.filling;
    doping  = param.doping;
    filling = param.filling;

    core::GridData fgrid = fgrid_in;

    const int Nk = param.Nk;
    const int dim = param.dim;
    const auto& q_range = param.q_range;
    const double T_K = param.T_K;
    const double Ef = param.Ef;
    const double polar_mu = param.polar_mu;
    const bool boundary_periodic = param.boundary_periodic;
    const double eta = param.eta;

    const int P = 2 * Nk + 1;
    const size_t Nrect = (size_t)P * (size_t)P;

    const int q_iq_min = q_range[0];
    const int q_iq_max = q_range[1];
    const int q_jq_min = q_range[2];
    const int q_jq_max = q_range[3];
    const int Nq1 = q_iq_max - q_iq_min + 1;
    const int Nq2 = q_jq_max - q_jq_min + 1;
    const size_t Nout = (size_t)Nq1 * (size_t)Nq2;

    const auto& evalsV = fgrid.get<std::vector<double>>("evals").v;
    const auto& kvecF  = fgrid.get<Vec2>("kvec").v;

    const std::vector<unsigned char>* insidePtr = nullptr;
    try { insidePtr = &fgrid.get<unsigned char>("inside").v; } catch (...) { insidePtr = nullptr; }

    // -------- polar branch: ensure occ_up/occ_dn exist, else fallback compute --------
    if (polar_mu != 0.0) {
        bool has_updn = false;
        try {
            (void)fgrid.get<std::vector<double>>("occ_up").v;
            (void)fgrid.get<std::vector<double>>("occ_dn").v;
            has_updn = true; 
        } catch (...) {
            has_updn = false;
        }

        if (!has_updn) {
            if (rank == 0) {
                std::cout << "[cal_chi] FAIL warmup\n";
            }

            auto& occUp_b = fgrid.add<std::vector<double>>("occ_up").v;
            auto& occDn_b = fgrid.add<std::vector<double>>("occ_dn").v;
            occUp_b.assign(Nrect, std::vector<double>(dim, 0.0));
            occDn_b.assign(Nrect, std::vector<double>(dim, 0.0));

            occTot_up = 0.0;
            occTot_dn = 0.0;

            const double Ef_up = Ef - 0.5 * polar_mu;
            const double Ef_dn = Ef + 0.5 * polar_mu;

            for (size_t idxk = 0; idxk < Nrect; ++idxk) {
                if (insidePtr && !(*insidePtr)[idxk]) continue;
                const auto& E = evalsV[idxk];

                std::vector<double> fu(dim, 0.0), fd(dim, 0.0);
                for (int b = 0; b < dim; ++b) {
                    fu[(size_t)b] = la::fermi(E[(size_t)b] - 0.5 * polar_mu, Ef_up, T_K);
                    fd[(size_t)b] = la::fermi(E[(size_t)b] + 0.5 * polar_mu, Ef_dn, T_K);
                    occTot_up += fu[(size_t)b];
                    occTot_dn += fd[(size_t)b];
                }
                occUp_b[idxk] = std::move(fu);
                occDn_b[idxk] = std::move(fd);
            }

            filling = (occTot_up + occTot_dn) / (2.0 * double(Nrect) * double(dim));
            doping  = la::filling_to_doping(filling);

        } else {
            if (rank == 0) {
                std::cout << "[cal_chi] SUCCESS warmup\n";
            }

            // injected case: recompute totals for metadata
            const auto& occUp_b = fgrid.get<std::vector<double>>("occ_up").v;
            const auto& occDn_b = fgrid.get<std::vector<double>>("occ_dn").v;

            occTot_up = 0.0;
            occTot_dn = 0.0;
            for (size_t idxk = 0; idxk < Nrect; ++idxk) {
                if (insidePtr && !(*insidePtr)[idxk]) continue;
                for (int b = 0; b < dim; ++b) {
                    occTot_up += occUp_b[idxk][(size_t)b];
                    occTot_dn += occDn_b[idxk][(size_t)b];
                }
            }

            filling = (occTot_up + occTot_dn) / (2.0 * double(Nrect) * double(dim));
            doping  = la::filling_to_doping(filling);
        }
    } else {
        // polar=0: use param filling metadata if needed
        occTot_up = filling * (2.0 * double(Nrect) * double(dim));
        occTot_dn = filling * (2.0 * double(Nrect) * double(dim));
    }

    // ---- output qgrid ----
    core::GridData& qgrid = result.qgrid;
    qgrid.resize(Nout); 
    qgrid.dx = fgrid.dx;
    qgrid.dy = fgrid.dy;

    auto& chiF   = qgrid.add<cd>("chi").v;
    auto& nPairF = qgrid.add<long long>("nKpair").v;
    auto& qvecF  = qgrid.add<Vec2>("qvec").v;
    auto& qxF    = qgrid.add<double>("qx").v;
    auto& qyF    = qgrid.add<double>("qy").v;

    // fill q-grid
    {
        size_t out = 0;
        for (int jq = q_jq_min; jq <= q_jq_max; ++jq) {
            for (int iq = q_iq_min; iq <= q_iq_max; ++iq, ++out) {
                qgrid.iq[out] = iq;
                qgrid.jq[out] = jq;
                Vec2 qvec = kvecF[fgrid.ij_to_idx_or_throw(iq, jq)];
                qvecF[out] = qvec;
                qxF[out] = qvec.x();
                qyF[out] = qvec.y();
                chiF[out] = cd(0.0, 0.0);
                nPairF[out] = 0;
            }
        }
    }

    // MPI split
    const auto [lin_s, lin_e] = rgmpi::block_1d_int((int)Nrect, rank, nprocs);

    // main loop
    for (size_t out = 0; out < Nout; ++out) {
        const int dq_iq = qgrid.iq[out];
        const int dq_jq = qgrid.jq[out];

        cd chi_local(0.0, 0.0);
        long long nPair_local = 0;

        for (size_t lin = (size_t)lin_s; lin < (size_t)lin_e; ++lin) {
            const size_t idxk1 = lin;
            if (insidePtr && !(*insidePtr)[idxk1]) continue;

            const auto [iq1, jq1] = fgrid.idx_to_ij(idxk1);
            int iq2 = iq1 + dq_iq;
            int jq2 = jq1 + dq_jq;

            if (boundary_periodic) {
                iq2 = wrap_q_index(iq2, Nk);
                jq2 = wrap_q_index(jq2, Nk);
            } else {
                if (iq2 < -Nk || iq2 > Nk) continue;
                if (jq2 < -Nk || jq2 > Nk) continue;
            }

            size_t idxk2 = 0;
            if (!fgrid.ij_to_idx(iq2, jq2, idxk2)) continue;
            if (insidePtr && !(*insidePtr)[idxk2]) continue;

            if (polar_mu == 0.0) {
                const auto& occTot_b = fgrid.get<std::vector<double>>("occ_band").v;
                const auto& E1 = evalsV[idxk1];
                const auto& E2 = evalsV[idxk2];
                const auto& f1 = occTot_b[idxk1];
                const auto& f2 = occTot_b[idxk2];

                for (int b = 0; b < dim; ++b) {
                    const double Eb = E1[(size_t)b];
                    const double fb = f1[(size_t)b];
                    for (int m = 0; m < dim; ++m) {
                        const double Em = E2[(size_t)m];
                        const double fm = f2[(size_t)m];
                        chi_local += (fm - fb) / cd(Eb - Em, eta);
                    }
                }
                nPair_local += (long long)dim * (long long)dim;

            } else {
                const auto& E1 = evalsV[idxk1];
                const auto& E2 = evalsV[idxk2];

                const auto& occUp_b = fgrid.get<std::vector<double>>("occ_up").v;
                const auto& occDn_b = fgrid.get<std::vector<double>>("occ_dn").v;
                const auto& f1u = occUp_b[idxk1];
                const auto& f1d = occDn_b[idxk1];
                const auto& f2u = occUp_b[idxk2];
                const auto& f2d = occDn_b[idxk2];

                for (int b = 0; b < dim; ++b) {
                    const double E1u = E1[(size_t)b] - 0.5 * polar_mu;
                    const double E1d = E1[(size_t)b] + 0.5 * polar_mu;
                    for (int m = 0; m < dim; ++m) {
                        const double E2u = E2[(size_t)m] - 0.5 * polar_mu;
                        const double E2d = E2[(size_t)m] + 0.5 * polar_mu;

                        chi_local += (f1u[(size_t)b] - f2u[(size_t)m]) / 4.0 / cd(-E1u + E2u, eta);
                        chi_local += (f1d[(size_t)b] - f2u[(size_t)m]) / 4.0 / cd(-E1d + E2u, eta);
                        chi_local += (f1u[(size_t)b] - f2d[(size_t)m]) / 4.0 / cd(-E1u + E2d, eta);
                        chi_local += (f1d[(size_t)b] - f2d[(size_t)m]) / 4.0 / cd(-E1d + E2d, eta);
                    }
                }
                nPair_local += (long long)dim * (long long)dim;
            }
        }

        double re = chi_local.real();
        double im = chi_local.imag();
        double re_sum = 0.0, im_sum = 0.0;
        long long np_sum = 0;

        rgmpi::reduce_sum(re, re_sum);
        rgmpi::reduce_sum(im, im_sum);
        rgmpi::reduce_sum(nPair_local, np_sum);

        if (rank == 0) {
            chiF[out]   = cd(re_sum, im_sum) / double(Nrect);
            nPairF[out] = np_sum;
        }
    }

#ifdef USE_MPI
    if (rgmpi::inited() && nprocs > 1) {
        std::vector<double> re_all(Nout, 0.0), im_all(Nout, 0.0);
        if (rank == 0) {
            for (size_t out = 0; out < Nout; ++out) {
                re_all[out] = chiF[out].real();
                im_all[out] = chiF[out].imag();
            }
        }
        MPI_Bcast(re_all.data(), (int)Nout, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(im_all.data(), (int)Nout, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(nPairF.data(), (int)Nout, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
        if (rank != 0) {
            for (size_t out = 0; out < Nout; ++out) chiF[out] = cd(re_all[out], im_all[out]);
        }
    }
#endif

    return result;
}


// ============================================================
// Path overload (simplified):
// 1) read base file -> fgrid
// 2) if polar_mu!=0: try locate mu_up/mu_dn fermiPatch and inject occ_up/occ_dn
//    if fail -> debug + fallback (kernel will self-compute occ_up/occ_dn)
// 3) call core kernel cal_chi_grid_from_fermiPatch(fgrid,param)
// ============================================================
// ============================================================
// Path overload (minimal print):
// Only print UP/DN success or failure
// ============================================================
inline rgio::cal_chi_result cal_chi_grid_from_fermiPatch(
    const std::string& file_path,
    rgio::cal_chi_param& param)
{
    int rank = 0, nprocs = 1;
    rgmpi::rank_size(rank, nprocs);

    namespace fs = std::filesystem;

    // (1) read base
    core::GridData fgrid = rgio::read_fermiPatch_grid_with_bands(file_path, param);

    if (param.polar_mu == 0.0) {
        return cal_chi_grid_from_fermiPatch(fgrid, param);
    }

    auto FAIL = [&]() -> rgio::cal_chi_result {
        return cal_chi_grid_from_fermiPatch(fgrid, param);
    };

    auto newest_fermiPatch_in_dir = [&](const fs::path& dir) -> std::string {
        if (!fs::exists(dir) || !fs::is_directory(dir)) return {};
        fs::file_time_type best_t{};
        fs::path best_p;
        bool has = false;
        for (const auto& ent : fs::directory_iterator(dir)) {
            if (!ent.is_regular_file()) continue;
            const fs::path p = ent.path();
            if (p.extension() != ".txt") continue;
            const std::string fn = p.filename().string();
            if (fn.rfind("fermiPatch_", 0) != 0) continue;
            const auto t = fs::last_write_time(p);
            if (!has || t > best_t) { best_t = t; best_p = p; has = true; }
        }
        return has ? best_p.string() : std::string{};
    };

    // locate shifted folders
    const fs::path p(file_path);
    const fs::path mu_dir = p.parent_path();
    const fs::path T_dir  = mu_dir.parent_path();
    const std::string mu_tag = mu_dir.filename().string();

    double muX = std::numeric_limits<double>::quiet_NaN();
    if (!rgio::parse_tag_value(mu_tag, "mu", muX)) {
        if (rank == 0) {
            std::cout << "[shift] UP MISS (cannot parse mu tag=" << mu_tag << ")\n";
            std::cout << "[shift] DN MISS (cannot parse mu tag=" << mu_tag << ")\n";
            std::cout << "[shift] UP FAIL\n";
            std::cout << "[shift] DN FAIL\n";
        }
        return FAIL();
    }

    const double mu_up = muX - 0.5 * param.polar_mu;
    const double mu_dn = muX + 0.5 * param.polar_mu;

    const fs::path up_dir = T_dir / ("mu" + rgio::tag6(mu_up));
    const fs::path dn_dir = T_dir / ("mu" + rgio::tag6(mu_dn));

    const std::string up_file = newest_fermiPatch_in_dir(up_dir);
    const std::string dn_file = newest_fermiPatch_in_dir(dn_dir);

    const bool up_found = !up_file.empty();
    const bool dn_found = !dn_file.empty();

    if (rank == 0) {
        std::cout << "[shift] UP " << (up_found ? "FOUND " : "MISS  ")
                  << (up_found ? up_file : up_dir.string()) << "\n";
        std::cout << "[shift] DN " << (dn_found ? "FOUND " : "MISS  ")
                  << (dn_found ? dn_file : dn_dir.string()) << "\n";
    }

    if (!up_found || !dn_found) {
        if (rank == 0) {
            std::cout << "[shift] UP " << (up_found ? "OK" : "FAIL") << "\n";
            std::cout << "[shift] DN " << (dn_found ? "OK" : "FAIL") << "\n";
        }
        return FAIL();
    }

    // read shifted grids
    bool up_ok = false, dn_ok = false;
    rgio::cal_chi_param pu = param, pd = param;
    core::GridData f_up, f_dn;

    try { f_up = rgio::read_fermiPatch_grid_with_bands(up_file, pu); up_ok = true; }
    catch (...) { up_ok = false; }

    try { f_dn = rgio::read_fermiPatch_grid_with_bands(dn_file, pd); dn_ok = true; }
    catch (...) { dn_ok = false; }

    if (rank == 0) {
        std::cout << "[shift] UP " << (up_ok ? "OK" : "FAIL") << "\n";
        std::cout << "[shift] DN " << (dn_ok ? "OK" : "FAIL") << "\n";
    }

    if (!up_ok || !dn_ok) return FAIL();

    // minimal consistency check
    if (pu.Nk != param.Nk || pd.Nk != param.Nk ||
        pu.dim != param.dim || pd.dim != param.dim ||
        f_up.size() != fgrid.size() || f_dn.size() != fgrid.size())
    {
        return FAIL();
    }

    // inject occ_up/occ_dn from shifted occ_band
    const auto& occb_up = f_up.get<std::vector<double>>("occ_band").v;
    const auto& occb_dn = f_dn.get<std::vector<double>>("occ_band").v;

    fgrid.add<std::vector<double>>("occ_up").v = occb_up;
    fgrid.add<std::vector<double>>("occ_dn").v = occb_dn;

    return cal_chi_grid_from_fermiPatch(fgrid, param);
}



} // namespace rg