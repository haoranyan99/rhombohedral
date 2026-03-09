// File: Source/Models/RG_ModelBase.h
#pragma once

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <complex>
#include <stdexcept>
#include <vector>
#include <ostream>
#include <fstream>


#ifdef USE_MPI
  #include <mpi.h>
#endif

#include "Common/DataContainers.h"
#include "PhysStruct/RG_Structure.h"
#include "LinearAlgebra/Constants.h"
#include "LinearAlgebra/MathFunctions.h"
#include "Util/MPI/mpi.h"

namespace rg {

class RG_ModelBase {
public:
    using Vec2 = Eigen::Vector2d;
    using Vec3 = Eigen::Vector3d;
    using cd   = std::complex<double>;

    RG_ModelBase() = default;
    explicit RG_ModelBase(const RG_Structure& st) : st_(st) {}
    virtual ~RG_ModelBase() = default;

    // ---- common members ----
    void   set_Dfield(double D) { Dfield_ = D; }
    double get_Dfield() const   { return Dfield_; }

    const RG_Structure& get_structure() const { return st_; }

    // ---- the one thing derived models must provide ----
    virtual Eigen::MatrixXcd build_Hk(const Vec2& k,
                                      bool enforce_hermitian = true) const = 0;

    // ---- common utilities (implemented once) ----
    virtual Eigen::VectorXd bands_at_k(const Vec2& k,
                                       bool enforce_hermitian = true) const;

    virtual Eigen::MatrixXd bands_along_path(const std::vector<Vec2>& klist,
                                             bool enforce_hermitian = true) const;

    virtual Eigen::MatrixXd bands_along_path(const core::PathData& path,
                                             bool enforce_hermitian = true) const;

    virtual core::SeriesData cal_dos_gaussian(const core::GridData& kpatch,
                                            double e_low, double e_high,
                                            int num_e,
                                            double eta,
                                            double T_K,
                                            bool enforce_hermitian = true) const;

    virtual core::GridData cal_chi_grid_Ef(const core::GridData& kpatch,
                                        int Nq,
                                        int q_center_iq,
                                        int q_center_jq,
                                        double Ef,
                                        double T_K,
                                        double polar_mu,
                                        double eta,
                                        bool enforce_hermitian,
                                        bool boundary_periodic,
                                        double nest_thr = 50) const;

    virtual core::GridData cal_chi_grid_from_fermiPatch(const std::string& in_path,
                                        int Nq,
                                        int q_center_iq,
                                        int q_center_jq,
                                        double eta,
                                        bool enforce_hermitian,
                                        bool boundary_periodic) const;

    virtual double find_EF_from_filling(double filling_target, double T_K, const core::GridData& kpatch) const;

    virtual double find_EF_from_filling_T0(double filling_target, const core::GridData& kpatch) const;

    virtual double find_filling_from_EF(double EF, double T_K, const core::GridData& kpatch) const;
    
    virtual double find_filling_from_EF_T0(double EF, const core::GridData& kpatch) const;

    virtual core::GridData cal_fermi_patch_from_mu(double EF, double T_K, const core::GridData& kpatch) const;

    virtual core::GridData cal_fermi_patch_from_filling(double filling, double T_K, const core::GridData& kpatch) const;
    
    virtual void write_info(std::ostream& os, const std::string& prefix = "# ") const
    {
        (void)os;
        (void)prefix;
    }

protected:
    static Eigen::MatrixXcd add_Dfield_(int norb, int Nl, double D);
    void write_info_common_(std::ostream& os, const std::string& prefix) const;

protected:
    RG_Structure st_;
    double Dfield_ = 0.0;
};

// ---------------- inline ----------------

inline Eigen::MatrixXcd RG_ModelBase::add_Dfield_(int norb, int Nl, double D)
{
    Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(norb, norb);

    std::vector<int> diag;
    if (Nl == 5)      diag = {-2,-2,-1,-1,0,0,1,1,2,2};
    else if (Nl == 6) diag = {-2,-2,-1,-1,0,0,1,1,2,2,3,3};
    else if (Nl == 4) diag = {-1,-1,0,0,1,1,2,2};
    else throw std::runtime_error("RG_SKModel::add_Dfield_: unsupported layer number (only 4/5/6 implemented)");

    if ((int)diag.size() != norb)
        throw std::runtime_error("RG_SKModel::add_Dfield_: diag size != norb (check layer number vs orbit count)");

    for (int i = 0; i < norb; ++i) M(i,i) = double(diag[i]) * D;
    return M;
}

inline Eigen::VectorXd RG_ModelBase::bands_at_k(const Vec2& k, bool enforce_hermitian) const {
    Eigen::MatrixXcd H = build_Hk(k, enforce_hermitian);

    // enforce Hermitian if requested (safety)
    if (enforce_hermitian) {
        H = 0.5 * (H + H.adjoint());
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(H);
    if (es.info() != Eigen::Success) {
        throw std::runtime_error("RG_ModelBase::bands_at_k: eigen solve failed");
    }
    return es.eigenvalues();
}

inline Eigen::MatrixXd RG_ModelBase::bands_along_path(const std::vector<Vec2>& klist,
                                                      bool enforce_hermitian) const {
    if (klist.empty()) return Eigen::MatrixXd();

    Eigen::VectorXd e0 = bands_at_k(klist.front(), enforce_hermitian);
    const int norb = (int)e0.size();

    Eigen::MatrixXd E((int)klist.size(), norb);
    E.row(0) = e0.transpose();

    for (int i = 1; i < (int)klist.size(); ++i) {
        E.row(i) = bands_at_k(klist[(size_t)i], enforce_hermitian).transpose();
    }
    return E;
}

inline Eigen::MatrixXd RG_ModelBase::bands_along_path(const core::PathData& path,
                                                      bool enforce_hermitian) const {
    return bands_along_path(path.k_list, enforce_hermitian);
}


inline core::SeriesData RG_ModelBase::cal_dos_gaussian(
    const core::GridData& kpatch,
    double e_low, double e_high,
    int num_e,
    double eta,
    double T_K,
    bool enforce_hermitian) const
{
    if (num_e < 2) rgmpi::abort_all("cal_dos_gaussian: num_e must be >= 2");
    if (!(e_high > e_low)) rgmpi::abort_all("cal_dos_gaussian: e_high must be > e_low");
    kpatch.assert_consistent();
    const size_t NkTot = kpatch.size();
    if (NkTot == 0) rgmpi::abort_all("cal_dos_gaussian: empty kpatch");

    int rank = 0, nprocs = 1;
    rgmpi::rank_size(rank, nprocs);

    const auto& kvec = kpatch.get<Vec2>("kvec").v;
    const unsigned char* inPtr = nullptr;
    try { inPtr = kpatch.get<unsigned char>("inside").v.data(); }
    catch (...) { inPtr = nullptr; }

    long long nK = 0;
    for (size_t i = 0; i < NkTot; ++i) if (!inPtr || inPtr[i]) ++nK;
    if (nK <= 0) rgmpi::abort_all("cal_dos_gaussian: no valid k points");
    const int norb = static_cast<int>(st_.atoms().size());
    if (norb <= 0) rgmpi::abort_all("cal_dos_gaussian: norb <= 0");
    const double dE = (e_high - e_low) / static_cast<double>(num_e - 1);
    if (!(eta > 0.0)) eta = 0.6 * dE;
    if (!(eta > 0.0)) rgmpi::abort_all("cal_dos_gaussian: eta must be > 0");

    core::SeriesData out;
    out.resize((size_t)num_e);
    out.dx = dE;

    // output fields
    auto& dos_field    = out.add<double>("dos").v;
    auto& fill_field   = out.add<double>("filling").v;
    auto& doping_field = out.add<double>("doping").v;

    for (int ie = 0; ie < num_e; ++ie) {
        const double E = e_low + dE * static_cast<double>(ie);
        out.x[(size_t)ie] = E;
        dos_field[(size_t)ie]    = 0.0;
        fill_field[(size_t)ie]   = 0.0;
        doping_field[(size_t)ie] = 0.0;
    }

    // distribute valid-k index rr_global in [0,nK) among ranks
    const long long base  = nK / nprocs;
    const long long rem   = nK % nprocs;
    const long long my_nK = base + (rank < rem ? 1 : 0);
    const long long my_s  = rank * base + std::min<long long>(rank, rem);
    const long long my_e  = my_s + my_nK;

    // local accumulators over energy grid
    std::vector<double> dos_local((size_t)num_e, 0.0);
    std::vector<double> fill_local((size_t)num_e, 0.0);

    // gaussian prefactors
    const double inv_sqrt2pi = 1.0 / std::sqrt(2.0 * la::pi);
    const double pref = inv_sqrt2pi / eta;
    const double inv2eta2 = 1.0 / (2.0 * eta * eta);

    long long rr_global = 0;
    for (size_t i = 0; i < NkTot; ++i) {
        if (inPtr && !inPtr[i]) continue;

        if (rr_global >= my_s && rr_global < my_e) {
            Eigen::MatrixXcd Hk = build_Hk(kvec[i], enforce_hermitian);
            if (enforce_hermitian) Hk = 0.5 * (Hk + Hk.adjoint());

            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(Hk);
            if (es.info() != Eigen::Success) {
                rgmpi::abort_all("cal_dos_gaussian: eigen decomposition failed");
            }

            const auto ev = es.eigenvalues(); // size = norb

            for (int ie = 0; ie < num_e; ++ie) {
                const double E = out.x[(size_t)ie];
                double acc_dos  = 0.0;
                double acc_fill = 0.0;

                for (int b = 0; b < norb; ++b) {
                    const double eb = ev(b);

                    // DOS Gaussian
                    const double de = E - eb;
                    acc_dos += std::exp(-(de * de) * inv2eta2) * pref;

                    // filling(EF=E)
                    acc_fill += la::fermi(eb, /*mu=*/E, T_K);
                }

                dos_local[(size_t)ie]  += acc_dos;
                fill_local[(size_t)ie] += acc_fill;
            }
        }

        ++rr_global;
    }

    // reduce + normalize
    const double invNk      = 1.0 / double(nK);
    const double invNkNorb  = 1.0 / (double(nK) * double(norb));

    const double a_Ang = st_.a(); // <<< use structure a(P) in Angstrom

    for (int ie = 0; ie < num_e; ++ie) {
        double dos_sum  = 0.0;
        double fill_sum = 0.0;

        rgmpi::allreduce_sum(dos_local[(size_t)ie],  dos_sum);
        rgmpi::allreduce_sum(fill_local[(size_t)ie], fill_sum);

        const double dos     = dos_sum  * invNk;
        const double filling = fill_sum * invNkNorb;
        const double doping  = la::filling_to_doping(filling, a_Ang);

        dos_field[(size_t)ie]    = dos;
        fill_field[(size_t)ie]   = filling;
        doping_field[(size_t)ie] = doping;
    }

    out.assert_consistent();
    return out;
}



inline core::GridData RG_ModelBase::cal_chi_grid_Ef(
    const core::GridData& kpatch,
    int Nq,
    int q_center_iq,
    int q_center_jq,
    double EF,
    double T_K,
    double polar_mu,
    double eta,
    bool enforce_hermitian,
    bool boundary_periodic,
    double nest_thr) const
{
    int rank = 0, nprocs = 1;
    rgmpi::rank_size(rank, nprocs);

    using Vec2 = Eigen::Vector2d;
    using cd   = std::complex<double>;

    if (Nq < 0)            rgmpi::abort_all("cal_chi_grid_Ef: Nq must be >= 0");
    if (!(eta > 0.0))      rgmpi::abort_all("cal_chi_grid_Ef: eta must be > 0");
    if (!(T_K >= 0.0))     rgmpi::abort_all("cal_chi_grid_Ef: T_K must be >= 0");
    if (!(nest_thr > 0.0)) rgmpi::abort_all("cal_chi_grid_Ef: nest_thr must be > 0");

    kpatch.assert_consistent();
    const size_t NkTot_sz = kpatch.size();
    if (NkTot_sz == 0) rgmpi::abort_all("cal_chi_grid_Ef: empty kpatch");

    // infer N1,N2 from iq/jq ranges (assumes full rectangle)
    auto infer_N1N2 = [&]() -> std::pair<int,int> {
        int iq_min =  1000000000, iq_max = -1000000000;
        int jq_min =  1000000000, jq_max = -1000000000;
        for (size_t t = 0; t < NkTot_sz; ++t) {
            iq_min = std::min(iq_min, kpatch.iq[t]);
            iq_max = std::max(iq_max, kpatch.iq[t]);
            jq_min = std::min(jq_min, kpatch.jq[t]);
            jq_max = std::max(jq_max, kpatch.jq[t]);
        }
        const int N1 = (iq_max - iq_min + 1);
        const int N2 = (jq_max - jq_min + 1);
        if ((long long)N1 * (long long)N2 != (long long)NkTot_sz) {
            rgmpi::abort_all("cal_chi_grid_Ef: cannot infer N1,N2 from iq/jq (size mismatch)");
        }
        return {N1, N2};
    };

    const auto [N1, N2] = infer_N1N2();

    int iq_min =  1000000000, jq_min =  1000000000;
    for (size_t t = 0; t < NkTot_sz; ++t) {
        iq_min = std::min(iq_min, kpatch.iq[t]);
        jq_min = std::min(jq_min, kpatch.jq[t]);
    }

    const auto& kvec = kpatch.get<Vec2>("kvec").v;

    // inside mask optional (used only for truncated)
    const std::vector<unsigned char>* inHexPtr = nullptr;
    std::vector<unsigned char> inAll;
    try {
        inHexPtr = &kpatch.get<unsigned char>("inside").v;
        if (inHexPtr->size() != NkTot_sz)
            rgmpi::abort_all("cal_chi_grid_Ef: inside size mismatch");
    } catch (...) {
        inAll.assign(NkTot_sz, (unsigned char)1);
        inHexPtr = &inAll;
    }
    const auto& inHex = *inHexPtr;

    // q basis vectors
    const Vec2 dq1 = kpatch.dx * st_.b1();
    const Vec2 dq2 = kpatch.dy * st_.b2();

    // lambdas
    auto idx_ij = [N2](int i, int j) -> int { return i * N2 + j; };

    auto mod_pos = [](int a, int m) -> int {
        int r = a % m;
        return (r < 0) ? (r + m) : r;
    };

    auto block_1d = [](int n, int r, int p) -> std::pair<int,int> {
        if (p <= 1) return {0, n};
        const int base = n / p;
        const int rem  = n % p;
        const int my_n = base + (r < rem ? 1 : 0);
        const int my_s = r * base + std::min(r, rem);
        return {my_s, my_s + my_n}; // [s,e)
    };

    const bool pm0 = (std::abs(polar_mu) < 1e-15);

    // ============================================================
    // (1) diagonalize H(k) on *my* i-block, and precompute occupations
    // ============================================================
    const int dim = 2 * st_.nLayer();

    Eigen::MatrixXd evals((int)NkTot_sz, dim);
    evals.setZero();

    // occ0 used when polar_mu==0
    Eigen::MatrixXd occ0((int)NkTot_sz, dim);
    occ0.setZero();

    // occ_up/occ_dn used when polar_mu!=0
    Eigen::MatrixXd occ_up((int)NkTot_sz, dim);
    Eigen::MatrixXd occ_dn((int)NkTot_sz, dim);
    occ_up.setZero();
    occ_dn.setZero();

    std::vector<unsigned char> has_diag(NkTot_sz, 0);

    const auto [i_s, i_e] = block_1d(N1, rank, nprocs);

    for (int i = i_s; i < i_e; ++i) {
        for (int j = 0; j < N2; ++j) {
            const int ik = idx_ij(i, j);
            if (!boundary_periodic && !inHex[(size_t)ik]) continue;

            Eigen::MatrixXcd Hk = build_Hk(kvec[(size_t)ik], enforce_hermitian);
            if (enforce_hermitian) Hk = 0.5 * (Hk + Hk.adjoint());

            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(Hk);
            if (es.info() != Eigen::Success)
                rgmpi::abort_all("cal_chi_grid_Ef: eigensolver failed");

            const auto& ev = es.eigenvalues();
            for (int n = 0; n < dim; ++n) {
                const double En = ev(n);
                evals(ik, n) = En;

                if (pm0) {
                    occ0(ik, n) = la::fermi(En, EF, T_K);
                } else {
                    occ_up(ik, n) = la::fermi(En + polar_mu, EF, T_K);
                    occ_dn(ik, n) = la::fermi(En - polar_mu, EF, T_K);
                }
            }
            has_diag[(size_t)ik] = 1;
        }
    }

#ifdef USE_MPI
    if (rgmpi::inited() && nprocs > 1) {
        MPI_Allreduce(MPI_IN_PLACE,
                      evals.data(),
                      (int)((long long)NkTot_sz * (long long)dim),
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        if (pm0) {
            MPI_Allreduce(MPI_IN_PLACE,
                          occ0.data(),
                          (int)((long long)NkTot_sz * (long long)dim),
                          MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        } else {
            MPI_Allreduce(MPI_IN_PLACE,
                          occ_up.data(),
                          (int)((long long)NkTot_sz * (long long)dim),
                          MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE,
                          occ_dn.data(),
                          (int)((long long)NkTot_sz * (long long)dim),
                          MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        }

        // assemble has_diag mask
        std::vector<int> has_i(NkTot_sz, 0), has_sum(NkTot_sz, 0);
        for (size_t t = 0; t < NkTot_sz; ++t) has_i[t] = (int)has_diag[t];

        MPI_Allreduce(has_i.data(), has_sum.data(), (int)NkTot_sz, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        for (size_t t = 0; t < NkTot_sz; ++t) has_diag[t] = (has_sum[t] > 0) ? 1 : 0;
    }
#endif

    // ============================================================
    // (2) allocate q grid around (q_center_iq, q_center_jq)
    // ============================================================
    const int qR = Nq;
    const int NqSide = 2 * qR + 1;
    const size_t Nout = (size_t)NqSide * (size_t)NqSide;

    core::GridData qgrid;
    qgrid.dx = kpatch.dx;
    qgrid.dy = kpatch.dy;
    qgrid.resize(Nout);

    auto& chiF   = qgrid.add<cd>("chi").v;
    auto& nPairF = qgrid.add<long long>("nKpair").v;
    auto& qxF    = qgrid.add<double>("qx").v;
    auto& qyF    = qgrid.add<double>("qy").v;

    {
        size_t out = 0;
        for (int djq = -qR; djq <= qR; ++djq) {
            for (int diq = -qR; diq <= qR; ++diq, ++out) {
                const int iq = q_center_iq + diq;
                const int jq = q_center_jq + djq;

                qgrid.iq[out] = iq;
                qgrid.jq[out] = jq;

                const Vec2 qvec = double(iq) * dq1 + double(jq) * dq2;
                qxF[out] = qvec.x();
                qyF[out] = qvec.y();

                chiF[out]   = cd(0.0, 0.0);
                nPairF[out] = 0;
            }
        }
    }

    // ============================================================
    // (3) compute chi(q): use precomputed occ arrays
    // ============================================================
    for (size_t out = 0; out < Nout; ++out) {
        const int iq = qgrid.iq[out];
        const int jq = qgrid.jq[out];

        cd chi_local(0.0, 0.0);
        long long nPair_local = 0;
        long long nNest_local = 0;

        if (boundary_periodic) {
            for (int i = i_s; i < i_e; ++i) {
                const int i2 = mod_pos(i + iq, N1);
                for (int j = 0; j < N2; ++j) {
                    const int j2 = mod_pos(j + jq, N2);

                    const int ik  = idx_ij(i,  j);
                    const int ik2 = idx_ij(i2, j2);
                    if (!has_diag[(size_t)ik] || !has_diag[(size_t)ik2]) continue;

                    for (int n = 0; n < dim; ++n) {
                        const double En = evals(ik, n);

                        if (pm0) {
                            const double fn = occ0(ik, n);

                            for (int m = 0; m < dim; ++m) {
                                const double Em = evals(ik2, m);
                                const double fm = occ0(ik2, m);

                                const double df = fm - fn;
                                const double dE = (En - Em);
                                const cd denom(dE, eta);

                                chi_local += df / denom;
                            }
                        } else {
                            const double fn_up = occ_up(ik, n);
                            const double fn_dn = occ_dn(ik, n);

                            for (int m = 0; m < dim; ++m) {
                                const double Em = evals(ik2, m);
                                const double fm_up = occ_up(ik2, m);
                                const double fm_dn = occ_dn(ik2, m);

                                const double df1 = fm_up - fn_up;
                                const double dE1 = (En - Em);
                                const cd denom1(dE1, eta);

                                const double df2 = fm_up - fn_dn;
                                const double dE2 = (En - Em - 2 * polar_mu);
                                const cd denom2(dE2, eta);

                                const double df3 = fm_dn - fn_up;
                                const double dE3 = (En - Em + 2 * polar_mu);
                                const cd denom3(dE3, eta);

                                const double df4 = fm_dn - fn_dn;
                                const double dE4 = (En - Em);
                                const cd denom4(dE4, eta);

                                chi_local += df1 / (4.0 * denom1);
                                chi_local += df2 / (4.0 * denom2);
                                chi_local += df3 / (4.0 * denom3);
                                chi_local += df4 / (4.0 * denom4);
                            }
                        }
                    }

                    nPair_local += (long long)dim * (long long)dim;
                }
            }
        } else {
            for (int i = i_s; i < i_e; ++i) {
                const int kiq = (i + iq_min);
                for (int j = 0; j < N2; ++j) {
                    const int kjq = (j + jq_min);

                    const int ik = idx_ij(i, j);
                    if (!has_diag[(size_t)ik]) continue;

                    const int kiq2 = kiq + iq;
                    const int kjq2 = kjq + jq;

                    const int i2 = kiq2 - iq_min;
                    const int j2 = kjq2 - jq_min;
                    if (i2 < 0 || i2 >= N1) continue;
                    if (j2 < 0 || j2 >= N2) continue;

                    const int ik2 = idx_ij(i2, j2);
                    if (!has_diag[(size_t)ik2]) continue;

                    for (int n = 0; n < dim; ++n) {
                        const double En = evals(ik, n);

                        if (pm0) {
                            const double fn = occ0(ik, n);

                            for (int m = 0; m < dim; ++m) {
                                const double Em = evals(ik2, m);
                                const double fm = occ0(ik2, m);

                                const double df = fm - fn;
                                const double dE = (En - Em);
                                const cd denom(dE, eta);

                                chi_local += df / denom;
                            }
                        } else {
                            const double fn_up = occ_up(ik, n);
                            const double fn_dn = occ_dn(ik, n);

                            for (int m = 0; m < dim; ++m) {
                                const double Em = evals(ik2, m);
                                const double fm_up = occ_up(ik2, m);
                                const double fm_dn = occ_dn(ik2, m);

                                const double df1 = fm_up - fn_up;
                                const double dE1 = (En - Em);
                                const cd denom1(dE1, eta);

                                const double df2 = fm_up - fn_dn;
                                const double dE2 = (En - Em - 2 * polar_mu);
                                const cd denom2(dE2, eta);

                                const double df3 = fm_dn - fn_up;
                                const double dE3 = (En - Em + 2 * polar_mu);
                                const cd denom3(dE3, eta);

                                const double df4 = fm_dn - fn_dn;
                                const double dE4 = (En - Em);
                                const cd denom4(dE4, eta);

                                chi_local += df1 / (4.0 * denom1);
                                chi_local += df2 / (4.0 * denom2);
                                chi_local += df3 / (4.0 * denom3);
                                chi_local += df4 / (4.0 * denom4);
                            }
                        }
                    }

                    nPair_local += (long long)dim * (long long)dim;
                }
            }
        }

        // reduce to rank0 (per-q point)
        double re = chi_local.real();
        double im = chi_local.imag();
        double re_sum = 0.0, im_sum = 0.0;
        long long np_sum = 0;

#ifdef USE_MPI
        if (rgmpi::inited() && nprocs > 1) {
            MPI_Reduce(&re, &re_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&im, &im_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&nPair_local, &np_sum, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        } else
#endif
        {
            re_sum = re;
            im_sum = im;
            np_sum = nPair_local;
        }

        if (rank == 0) {
            cd chi_global(re_sum, im_sum);
            chi_global /= double(NkTot_sz);

            chiF[out]   = chi_global;
            nPairF[out] = np_sum;
        }
    }

#ifdef USE_MPI
    // broadcast final chi to all ranks to keep behavior consistent
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
            for (size_t out = 0; out < Nout; ++out) {
                chiF[out] = cd(re_all[out], im_all[out]);
            }
        }
    }
#endif

    return qgrid;
}


inline core::GridData RG_ModelBase::cal_chi_grid_from_fermiPatch(
    const std::string& in_path,
    int Nq,
    int q_center_iq,
    int q_center_jq,
    double eta,
    bool enforce_hermitian,      // kept for signature compatibility; not used
    bool boundary_periodic) const
{
    int rank = 0, nprocs = 1;
    rgmpi::rank_size(rank, nprocs);

    using Vec2 = Eigen::Vector2d;
    using cd   = std::complex<double>;

    (void)enforce_hermitian;

    if (Nq < 0)        rgmpi::abort_all("RG_ModelBase::cal_chi_grid_from_fermiPatch_file_interband: Nq must be >= 0");
    if (!(eta > 0.0))  rgmpi::abort_all("RG_ModelBase::cal_chi_grid_from_fermiPatch_file_interband: eta must be > 0");

    // ----------------------------
    // helpers (lambda only)
    // ----------------------------
    auto trim = [](const std::string& s)->std::string {
        size_t a = 0, b = s.size();
        while (a < b && std::isspace((unsigned char)s[a])) ++a;
        while (b > a && std::isspace((unsigned char)s[b-1])) --b;
        return s.substr(a, b-a);
    };
    auto is_hash_line = [&](const std::string& s)->bool {
        std::string t = trim(s);
        return (!t.empty() && t[0] == '#');
    };
    auto mod_pos = [](int a, int m)->int {
        int r = a % m;
        return (r < 0) ? (r + m) : r;
    };
    auto idx_ij = [](int jstride, int i, int j)->int { return i * jstride + j; };
    auto block_1d = [](int n, int r, int p)->std::pair<int,int> {
        if (p <= 1) return {0, n};
        const int base = n / p;
        const int rem  = n % p;
        const int my_n = base + (r < rem ? 1 : 0);
        const int my_s = r * base + std::min(r, rem);
        return {my_s, my_s + my_n}; // [s,e)
    };

    // ----------------------------
    // (0) read numeric lines
    // Expected numeric columns:
    // idx iq jq kx ky occ E0..E{dim-1} f0..f{dim-1}
    // infer dim from ncol = 6 + 2*dim
    // ----------------------------
    std::ifstream fin(in_path);
    if (!fin) rgmpi::abort_all("...: cannot open input file");

    std::string line;
    int ncol = -1;
    int dim  = -1;

    std::vector<int>    iq_raw, jq_raw;
    std::vector<double> kx_raw, ky_raw;
    std::vector<double> E_raw; // flattened per row: [t*dim + b]
    std::vector<double> f_raw; // flattened per row: [t*dim + b]

    iq_raw.reserve(200000);
    jq_raw.reserve(200000);
    kx_raw.reserve(200000);
    ky_raw.reserve(200000);

    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        if (is_hash_line(line)) continue;

        std::istringstream iss(line);
        std::vector<double> v;
        v.reserve(64);
        double x;
        while (iss >> x) v.push_back(x);
        if (v.empty()) continue;

        if (ncol < 0) {
            ncol = (int)v.size();
            if (ncol < 8) rgmpi::abort_all("...: numeric line has <8 columns");
            if ((ncol - 6) % 2 != 0) rgmpi::abort_all("...: malformed, expect ncol = 6 + 2*dim");
            dim = (ncol - 6) / 2;
            if (dim <= 0) rgmpi::abort_all("...: inferred dim <= 0");
            E_raw.reserve((size_t)200000 * (size_t)dim);
            f_raw.reserve((size_t)200000 * (size_t)dim);
        } else {
            if ((int)v.size() != ncol) rgmpi::abort_all("...: inconsistent column count");
        }

        const int iqv = (int)std::llround(v[1]);
        const int jqv = (int)std::llround(v[2]);

        iq_raw.push_back(iqv);
        jq_raw.push_back(jqv);
        kx_raw.push_back(v[3]);
        ky_raw.push_back(v[4]);

        const int E0 = 6;
        const int F0 = 6 + dim;
        for (int b = 0; b < dim; ++b) E_raw.push_back(v[E0 + b]);
        for (int b = 0; b < dim; ++b) f_raw.push_back(v[F0 + b]);
    }
    fin.close();

    const size_t NkTot_sz = iq_raw.size();
    if (NkTot_sz == 0) rgmpi::abort_all("...: no numeric data lines");

    // ----------------------------
    // (1) infer iq/jq rectangle, build first-hit grid
    // ----------------------------
    int iq_min =  1000000000, iq_max = -1000000000;
    int jq_min =  1000000000, jq_max = -1000000000;
    for (size_t t = 0; t < NkTot_sz; ++t) {
        iq_min = std::min(iq_min, iq_raw[t]);
        iq_max = std::max(iq_max, iq_raw[t]);
        jq_min = std::min(jq_min, jq_raw[t]);
        jq_max = std::max(jq_max, jq_raw[t]);
    }
    const int N1 = iq_max - iq_min + 1;
    const int N2 = jq_max - jq_min + 1;
    if (N1 <= 0 || N2 <= 0) rgmpi::abort_all("...: infer N1/N2 failed");

    const size_t Nrect = (size_t)N1 * (size_t)N2;

    std::vector<unsigned char> has(Nrect, 0);
    std::vector<double> kxG(Nrect, std::numeric_limits<double>::quiet_NaN());
    std::vector<double> kyG(Nrect, std::numeric_limits<double>::quiet_NaN());
    std::vector<double> EG((size_t)dim * Nrect, std::numeric_limits<double>::quiet_NaN());
    std::vector<double> fG((size_t)dim * Nrect, std::numeric_limits<double>::quiet_NaN());

    auto atE = [&](int b, size_t lin)->double& { return EG[(size_t)b * Nrect + lin]; };
    auto atf = [&](int b, size_t lin)->double& { return fG[(size_t)b * Nrect + lin]; };

    for (size_t t = 0; t < NkTot_sz; ++t) {
        const int i = iq_raw[t] - iq_min;
        const int j = jq_raw[t] - jq_min;
        if (i < 0 || i >= N1 || j < 0 || j >= N2) continue;
        const size_t lin = (size_t)idx_ij(N2, i, j);
        if (has[lin]) continue;

        has[lin] = 1;
        kxG[lin] = kx_raw[t];
        kyG[lin] = ky_raw[t];

        const size_t base = t * (size_t)dim;
        for (int b = 0; b < dim; ++b) {
            atE(b, lin) = E_raw[base + (size_t)b];
            atf(b, lin) = f_raw[base + (size_t)b];
        }
    }

    std::vector<int> lin_list;
    lin_list.reserve(Nrect);
    for (size_t lin = 0; lin < Nrect; ++lin)
        if (has[lin]) lin_list.push_back((int)lin);
    if (lin_list.empty()) rgmpi::abort_all("...: no valid points after gridding");

    // ----------------------------
    // (2) infer dq1/dq2 from index neighbors around (0,0)
    // ----------------------------
    Vec2 dq1 = Vec2::Constant(std::numeric_limits<double>::quiet_NaN());
    Vec2 dq2 = Vec2::Constant(std::numeric_limits<double>::quiet_NaN());

    auto try_lin = [&](int iq, int jq)->std::pair<bool,size_t> {
        const int i = iq - iq_min;
        const int j = jq - jq_min;
        if (i < 0 || i >= N1 || j < 0 || j >= N2) return {false, 0};
        const size_t lin = (size_t)idx_ij(N2, i, j);
        if (!has[lin]) return {false, 0};
        return {true, lin};
    };

    {
        auto [ok00, lin00] = try_lin(0,0);
        auto [ok10, lin10] = try_lin(1,0);
        auto [ok01, lin01] = try_lin(0,1);
        if (ok00 && ok10 && ok01 &&
            std::isfinite(kxG[lin00]) && std::isfinite(kyG[lin00]) &&
            std::isfinite(kxG[lin10]) && std::isfinite(kyG[lin10]) &&
            std::isfinite(kxG[lin01]) && std::isfinite(kyG[lin01])) {

            const Vec2 k00(kxG[lin00], kyG[lin00]);
            dq1 = Vec2(kxG[lin10], kyG[lin10]) - k00;
            dq2 = Vec2(kxG[lin01], kyG[lin01]) - k00;
        }
    }

    // ----------------------------
    // (3) allocate q grid around center
    // ----------------------------
    const int qR = Nq;
    const int NqSide = 2 * qR + 1;
    const size_t Nout = (size_t)NqSide * (size_t)NqSide;

    core::GridData qgrid;
    qgrid.resize(Nout);

    auto& chiF   = qgrid.add<cd>("chi").v;
    auto& nPairF = qgrid.add<long long>("nKpair").v;
    auto& qxF    = qgrid.add<double>("qx").v;
    auto& qyF    = qgrid.add<double>("qy").v;

    {
        size_t out = 0;
        for (int djq = -qR; djq <= qR; ++djq) {
            for (int diq = -qR; diq <= qR; ++diq, ++out) {
                const int iq = q_center_iq + diq;
                const int jq = q_center_jq + djq;

                qgrid.iq[out] = iq;
                qgrid.jq[out] = jq;

                if (std::isfinite(dq1.x()) && std::isfinite(dq2.x())) {
                    const Vec2 qvec = double(iq) * dq1 + double(jq) * dq2;
                    qxF[out] = qvec.x();
                    qyF[out] = qvec.y();
                } else {
                    qxF[out] = std::numeric_limits<double>::quiet_NaN();
                    qyF[out] = std::numeric_limits<double>::quiet_NaN();
                }

                chiF[out]   = cd(0.0, 0.0);
                nPairF[out] = 0;
            }
        }
    }

    // ----------------------------
    // (4) MPI split i-block (same style)
    // ----------------------------
    const auto [i_s, i_e] = block_1d(N1, rank, nprocs);

    auto lin_to_ij = [&](int lin)->std::pair<int,int> {
        int i = lin / N2;
        int j = lin - i * N2;
        return {i,j};
    };

    // ----------------------------
    // (5) compute interband chi(q): only b != m
    // ----------------------------
    for (size_t out = 0; out < Nout; ++out) {
        const int iq = qgrid.iq[out];
        const int jq = qgrid.jq[out];

        cd chi_local(0.0, 0.0);
        long long nPair_local = 0;

        for (int lin : lin_list) {
            auto [i, j] = lin_to_ij(lin);
            if (i < i_s || i >= i_e) continue;

            int i2 = i + iq;
            int j2 = j + jq;

            if (boundary_periodic) {
                i2 = mod_pos(i2, N1);
                j2 = mod_pos(j2, N2);
            } else {
                if (i2 < 0 || i2 >= N1) continue;
                if (j2 < 0 || j2 >= N2) continue;
            }

            const size_t lin2 = (size_t)idx_ij(N2, i2, j2);
            if (!has[lin2]) continue;

            for (int b = 0; b < dim; ++b) {
                const double Eb = atE(b, (size_t)lin);
                const double fb = atf(b, (size_t)lin);
                if (!std::isfinite(Eb) || !std::isfinite(fb)) continue;

                for (int m = 0; m < dim; ++m) {
                    if (m == b) continue; // <<< interband only

                    const double Em = atE(m, lin2);
                    const double fm = atf(m, lin2);
                    if (!std::isfinite(Em) || !std::isfinite(fm)) continue;

                    const double df = fb - fm;
                    const double dE = Eb - Em;

                    chi_local += df / cd(dE, eta);
                    nPair_local += 1;
                }
            }
        }

        double re = chi_local.real(), im = chi_local.imag();
        double re_sum = 0.0, im_sum = 0.0;
        long long np_sum = 0;

#ifdef USE_MPI
        if (rgmpi::inited() && nprocs > 1) {
            MPI_Reduce(&re, &re_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&im, &im_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&nPair_local, &np_sum, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        } else
#endif
        {
            re_sum = re; im_sum = im; np_sum = nPair_local;
        }

        if (rank == 0) {
            chiF[out]   = cd(re_sum, im_sum);
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
            for (size_t out = 0; out < Nout; ++out) {
                chiF[out] = cd(re_all[out], im_all[out]);
            }
        }
    }
#endif

    return qgrid;
}




// ------------------------------------------------------------------
// MPI / non-MPI version: bisection on filling(EF)-filling_target
// ------------------------------------------------------------------
inline double RG_ModelBase::find_EF_from_filling(double filling_target,
                                                double T_K,
                                                const core::GridData& kpatch) const
{
    if (T_K == 0.0) {
        return find_EF_from_filling_T0(filling_target, kpatch);
    }

    constexpr double E0       = 0.2;
    constexpr double EMAX     = 10.0;
    constexpr double expand   = 2.0;
    constexpr int    max_iter = 80;
    constexpr double tol_f    = 1e-9;
    constexpr double tol_EF   = 1e-6;

    if (!(filling_target >= 0.0 && filling_target <= 1.0)) {
        rgmpi::abort_all("find_EF_from_filling: filling_target must be in [0,1]");
    }

    kpatch.assert_consistent();
    const size_t NkTot_sz = kpatch.size();
    if (NkTot_sz == 0) rgmpi::abort_all("find_EF_from_filling: empty kpatch");
    if (!(T_K >= 0.0)) rgmpi::abort_all("find_EF_from_filling: T_K must be >= 0");

    int rank = 0, nprocs = 1;
    rgmpi::rank_size(rank, nprocs);

    const bool enforce_hermitian = true;

    const int norb = static_cast<int>(st_.atoms().size());
    const auto& kvec = kpatch.get<rg::Vec2>("kvec").v;

    const unsigned char* inHexPtr = nullptr;
    try { inHexPtr = kpatch.get<unsigned char>("inside").v.data(); }
    catch (...) { inHexPtr = nullptr; }

    // -------- count valid k points (global, identical on all ranks) --------
    long long nK_global = 0;
    for (size_t i = 0; i < NkTot_sz; ++i) {
        if (!inHexPtr || inHexPtr[i]) ++nK_global;
    }
    if (nK_global <= 0) rgmpi::abort_all("find_EF_from_filling: no valid k points");

    // -------- distribute rr_global in [0, nK_global) across ranks --------
    const long long base = nK_global / nprocs;
    const long long rem  = nK_global % nprocs;
    const long long my_nK = base + (rank < rem ? 1 : 0);
    const long long my_s  = rank * base + std::min<long long>(rank, rem);
    const long long my_e  = my_s + my_nK;

    Eigen::MatrixXd evals_local((int)my_nK, norb);

    // -------- diagonalize only local share --------
    long long rr_global = 0;
    int rr_local = 0;

    for (size_t id = 0; id < NkTot_sz; ++id) {
        if (inHexPtr && !inHexPtr[id]) continue;

        if (rr_global >= my_s && rr_global < my_e) {
            Eigen::MatrixXcd Hk = build_Hk(kvec[id], enforce_hermitian);
            Hk = 0.5 * (Hk + Hk.adjoint());

            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(Hk);
            if (es.info() != Eigen::Success) {
                rgmpi::abort_all("find_EF_from_filling: eigen decomposition failed");
            }

            evals_local.row(rr_local++) = es.eigenvalues().transpose();
        }
        ++rr_global;
    }
    if (rr_local != (int)my_nK) {
        rgmpi::abort_all("find_EF_from_filling: MPI rr_local mismatch");
    }

    // global normalization
    const double invDen = 1.0 / (double(nK_global) * double(norb));

    auto filling_of_EF = [&](double EF) -> double {
        double acc_local = 0.0;
        for (int p = 0; p < (int)my_nK; ++p) {
            for (int b = 0; b < norb; ++b) {
                acc_local += la::fermi(evals_local(p, b), EF, T_K);
            }
        }
        double acc = 0.0;
        rgmpi::allreduce_sum(acc_local, acc);
        return acc * invDen;
    };

    auto g = [&](double EF) -> double { return filling_of_EF(EF) - filling_target; };

    // -------- bracket root --------
    double lo = -E0;
    double hi = +E0;
    double f_lo = g(lo);
    double f_hi = g(hi);

    while (f_lo * f_hi > 0.0 && (std::abs(lo) < EMAX) && (std::abs(hi) < EMAX)) {
        lo *= expand;
        hi *= expand;
        f_lo = g(lo);
        f_hi = g(hi);
    }

    if (f_lo * f_hi > 0.0) {
        std::ostringstream oss;
        oss << "find_EF_from_filling: failed to bracket root\n"
            << "  filling_target=" << filling_target << "\n"
            << "  last interval=[" << lo << ", " << hi << "]\n"
            << "  filling(lo)=" << (f_lo + filling_target) << "\n"
            << "  filling(hi)=" << (f_hi + filling_target) << "\n"
            << "  Consider increasing internal EMAX if bandwidth is larger.";
        rgmpi::abort_all(oss.str());
    }

    // -------- bisection --------
    for (int it = 0; it < max_iter; ++it) {
        if ((hi - lo) < tol_EF) return 0.5 * (hi + lo);

        const double mid   = 0.5 * (lo + hi);
        const double f_mid = g(mid);

        if (std::abs(f_mid) < tol_f) return mid;

        if (f_lo * f_mid <= 0.0) {
            hi = mid; f_hi = f_mid;
        } else {
            lo = mid; f_lo = f_mid;
        }
    }

    return 0.5 * (lo + hi);
}

// ------------------------------------------------------------------
// MPI T=0 version (no sorting): bisection on count(E<EF)
// ------------------------------------------------------------------
inline double RG_ModelBase::find_EF_from_filling_T0(double filling_target,
                                                   const core::GridData& kpatch) const
{
    constexpr int    max_iter = 80;
    constexpr double tol_EF   = 1e-9;
    constexpr double pad_E    = 1e-6;

    if (!(filling_target >= 0.0 && filling_target <= 1.0)) {
        rgmpi::abort_all("find_EF_from_filling_T0: filling_target must be in [0,1]");
    }

    kpatch.assert_consistent();
    const size_t NkTot_sz = kpatch.size();
    if (NkTot_sz == 0) rgmpi::abort_all("find_EF_from_filling_T0: empty kpatch");

    int rank = 0, nprocs = 1;
    rgmpi::rank_size(rank, nprocs);

    const bool enforce_hermitian = true;

    const int norb = static_cast<int>(st_.atoms().size());
    const auto& kvec = kpatch.get<rg::Vec2>("kvec").v;

    const unsigned char* inHexPtr = nullptr;
    try { inHexPtr = kpatch.get<unsigned char>("inside").v.data(); }
    catch (...) { inHexPtr = nullptr; }

    long long nK_global = 0;
    for (size_t i = 0; i < NkTot_sz; ++i) {
        if (!inHexPtr || inHexPtr[i]) ++nK_global;
    }
    if (nK_global <= 0) rgmpi::abort_all("find_EF_from_filling_T0: no valid k points");

    const long long Ntot = nK_global * (long long)norb;
    if (Ntot <= 0) rgmpi::abort_all("find_EF_from_filling_T0: invalid total states");

    // -------- distribute valid-k rr_global across ranks --------
    const long long base = nK_global / nprocs;
    const long long rem  = nK_global % nprocs;
    const long long my_nK = base + (rank < rem ? 1 : 0);
    const long long my_s  = rank * base + std::min<long long>(rank, rem);
    const long long my_e  = my_s + my_nK;

    // Precompute local eigenvalues, and local Emin/Emax
    Eigen::MatrixXd evals_local((int)my_nK, norb);

    double Emin_local = +std::numeric_limits<double>::infinity();
    double Emax_local = -std::numeric_limits<double>::infinity();

    long long rr_global = 0;
    int rr_local = 0;

    for (size_t id = 0; id < NkTot_sz; ++id) {
        if (inHexPtr && !inHexPtr[id]) continue;

        if (rr_global >= my_s && rr_global < my_e) {
            Eigen::MatrixXcd Hk = build_Hk(kvec[id], enforce_hermitian);
            Hk = 0.5 * (Hk + Hk.adjoint());

            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(Hk);
            if (es.info() != Eigen::Success) {
                rgmpi::abort_all("find_EF_from_filling_T0: eigen decomposition failed");
            }

            const auto ev = es.eigenvalues();
            evals_local.row(rr_local++) = ev.transpose();

            Emin_local = std::min(Emin_local, ev.minCoeff());
            Emax_local = std::max(Emax_local, ev.maxCoeff());
        }
        ++rr_global;
    }

    if (rr_local != (int)my_nK) {
        rgmpi::abort_all("find_EF_from_filling_T0: MPI rr_local mismatch");
    }

    double Emin = 0.0, Emax = 0.0;
    rgmpi::allreduce_min(Emin_local, Emin);
    rgmpi::allreduce_max(Emax_local, Emax);

    // trivial ends (use global Emin/Emax)
    if (filling_target <= 0.0) return Emin - pad_E;
    if (filling_target >= 1.0) return Emax + pad_E;

    const double invDen = 1.0 / (double)Ntot;

    auto count_below = [&](double EF) -> double {
        long long cnt_local = 0;
        for (int p = 0; p < (int)my_nK; ++p) {
            for (int b = 0; b < norb; ++b) {
                if (evals_local(p, b) < EF) ++cnt_local;
            }
        }
        long long cnt = 0;
        rgmpi::allreduce_sum(cnt_local, cnt);
        return double(cnt) * invDen;
    };

    auto g = [&](double EF) -> double { return count_below(EF) - filling_target; };

    double lo = Emin - pad_E;
    double hi = Emax + pad_E;

    double f_lo = g(lo);
    double f_hi = g(hi);

    // (rare) safety expand
    int safe_iter = 0;
    while (f_lo * f_hi > 0.0 && safe_iter++ < 10) {
        const double w = (hi - lo);
        lo -= w;
        hi += w;
        f_lo = g(lo);
        f_hi = g(hi);
    }
    if (f_lo * f_hi > 0.0) {
        rgmpi::abort_all("find_EF_from_filling_T0: failed to bracket root");
    }

    for (int it = 0; it < max_iter; ++it) {
        if ((hi - lo) < tol_EF) return 0.5 * (hi + lo);

        const double mid   = 0.5 * (lo + hi);
        const double f_mid = g(mid);

        if (f_lo * f_mid <= 0.0) {
            hi = mid; f_hi = f_mid;
        } else {
            lo = mid; f_lo = f_mid;
        }
    }

    return 0.5 * (lo + hi);
}

// ------------------------------------------------------------------
// find_filling_from_EF (finite T): compute filling(EF) = < f(E,EF,T) >
// MPI-aware: diagonalize only local share; reduce sum across ranks
// ------------------------------------------------------------------
inline double RG_ModelBase::find_filling_from_EF(double EF,
                                           double T_K,
                                           const core::GridData& kpatch) const
{
    if (T_K == 0.0) {
        return find_filling_from_EF_T0(EF, kpatch);
    }

    kpatch.assert_consistent();
    const size_t NkTot_sz = kpatch.size();
    if (NkTot_sz == 0) rgmpi::abort_all("filling_from_EF: empty kpatch");
    if (!(T_K > 0.0)) rgmpi::abort_all("filling_from_EF: T_K must be > 0 for finite-T version");

    int rank = 0, nprocs = 1;
    rgmpi::rank_size(rank, nprocs);

    const bool enforce_hermitian = true;

    const int norb = static_cast<int>(st_.atoms().size());
    const auto& kvec = kpatch.get<rg::Vec2>("kvec").v;

    const unsigned char* inHexPtr = nullptr;
    try { inHexPtr = kpatch.get<unsigned char>("inside").v.data(); }
    catch (...) { inHexPtr = nullptr; }

    // count valid k (global, same on all ranks)
    long long nK_global = 0;
    for (size_t i = 0; i < NkTot_sz; ++i) {
        if (!inHexPtr || inHexPtr[i]) ++nK_global;
    }
    if (nK_global <= 0) rgmpi::abort_all("filling_from_EF: no valid k points");

    // distribute rr_global in [0, nK_global) across ranks
    const long long base  = nK_global / nprocs;
    const long long rem   = nK_global % nprocs;
    const long long my_nK = base + (rank < rem ? 1 : 0);
    const long long my_s  = rank * base + std::min<long long>(rank, rem);
    const long long my_e  = my_s + my_nK;

    // local accumulation of sum_b f(E_b(k),EF,T)
    double acc_local = 0.0;

    long long rr_global = 0;
    for (size_t id = 0; id < NkTot_sz; ++id) {
        if (inHexPtr && !inHexPtr[id]) continue;

        if (rr_global >= my_s && rr_global < my_e) {
            Eigen::MatrixXcd Hk = build_Hk(kvec[id], enforce_hermitian);
            Hk = 0.5 * (Hk + Hk.adjoint());

            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(Hk);
            if (es.info() != Eigen::Success) {
                rgmpi::abort_all("filling_from_EF: eigen decomposition failed");
            }

            const auto ev = es.eigenvalues();
            for (int b = 0; b < norb; ++b) {
                acc_local += la::fermi(ev(b), EF, T_K);
            }
        }
        ++rr_global;
    }

    // Allreduce sum across ranks
    double acc = 0.0;
    rgmpi::allreduce_sum(acc_local, acc);

    // normalize by (nK_global * norb)
    const double invDen = 1.0 / (double(nK_global) * double(norb));
    return acc * invDen;
}

// ------------------------------------------------------------------
// find_filling_from_EF_T0: filling(EF) = count(E < EF) / (nK_global*norb)
// MPI-aware: diagonalize local share and reduce count
// ------------------------------------------------------------------
inline double RG_ModelBase::find_filling_from_EF_T0(double EF,
                                              const core::GridData& kpatch) const
{
    kpatch.assert_consistent();
    const size_t NkTot_sz = kpatch.size();
    if (NkTot_sz == 0) rgmpi::abort_all("filling_from_EF_T0: empty kpatch");

    int rank = 0, nprocs = 1;
    rgmpi::rank_size(rank, nprocs);

    const bool enforce_hermitian = true;

    const int norb = static_cast<int>(st_.atoms().size());
    const auto& kvec = kpatch.get<rg::Vec2>("kvec").v;

    const unsigned char* inHexPtr = nullptr;
    try { inHexPtr = kpatch.get<unsigned char>("inside").v.data(); }
    catch (...) { inHexPtr = nullptr; }

    long long nK_global = 0;
    for (size_t i = 0; i < NkTot_sz; ++i) {
        if (!inHexPtr || inHexPtr[i]) ++nK_global;
    }
    if (nK_global <= 0) rgmpi::abort_all("filling_from_EF_T0: no valid k points");

    const long long Ntot = nK_global * (long long)norb;
    if (Ntot <= 0) rgmpi::abort_all("filling_from_EF_T0: invalid total states");

    // distribute valid-k rr_global across ranks
    const long long base  = nK_global / nprocs;
    const long long rem   = nK_global % nprocs;
    const long long my_nK = base + (rank < rem ? 1 : 0);
    const long long my_s  = rank * base + std::min<long long>(rank, rem);
    const long long my_e  = my_s + my_nK;

    long long cnt_local = 0;

    long long rr_global = 0;
    for (size_t id = 0; id < NkTot_sz; ++id) {
        if (inHexPtr && !inHexPtr[id]) continue;

        if (rr_global >= my_s && rr_global < my_e) {
            Eigen::MatrixXcd Hk = build_Hk(kvec[id], enforce_hermitian);
            Hk = 0.5 * (Hk + Hk.adjoint());

            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(Hk);
            if (es.info() != Eigen::Success) {
                rgmpi::abort_all("filling_from_EF_T0: eigen decomposition failed");
            }

            const auto ev = es.eigenvalues();
            for (int b = 0; b < norb; ++b) {
                if (ev(b) < EF) ++cnt_local;
            }
        }
        ++rr_global;
    }

    long long cnt = 0;
    rgmpi::allreduce_sum(cnt_local, cnt);

    return double(cnt) / double(Ntot);
}



// RG_ModelBase::cal_fermi_patch_from_mu
inline core::GridData RG_ModelBase::cal_fermi_patch_from_mu(double EF, double T_K, const core::GridData& kpatch) const
{
    constexpr bool debug_print_summary = false;

    int rank = 0, nprocs = 1;
    rgmpi::rank_size(rank, nprocs);

    kpatch.assert_consistent();
    const size_t NkTot_sz = kpatch.size();
    if (NkTot_sz == 0) rgmpi::abort_all("cal_fermi_patch_from_mu: empty patch");
    if (!(T_K >= 0.0)) rgmpi::abort_all("cal_fermi_patch_from_mu: T_K must be >= 0");
    if (!std::isfinite(EF)) rgmpi::abort_all("cal_fermi_patch_from_mu: EF must be finite");

    const auto& kvecP = kpatch.get<rg::Vec2>("kvec").v;

    // inside mask optional
    const unsigned char* inHexPtr = nullptr;
    try { inHexPtr = kpatch.get<unsigned char>("inside").v.data(); }
    catch (...) { inHexPtr = nullptr; }

    // infer band dimension from H(k)
    Eigen::MatrixXcd H0 = build_Hk(kvecP[0], true);
    const int dim = int(H0.rows());
    if (H0.cols() != dim) rgmpi::abort_all("cal_fermi_patch_from_mu: H is not square");

    // output grid
    core::GridData grid;
    grid.resize(NkTot_sz);

    auto& kvecF   = grid.add<rg::Vec2>("kvec").v;
    auto& inHexF  = grid.add<unsigned char>("inside").v;
    auto& EFfield = grid.add<double>("EF_used").v;
    auto& occAvgF = grid.add<double>("occ_k_avg").v;

    // debug per k
    auto& evalsF  = grid.add<std::vector<double>>("evals").v;      // N, each dim
    auto& occbF   = grid.add<std::vector<double>>("occ_band").v;   // N, each dim

    const double NaN = std::numeric_limits<double>::quiet_NaN();

    // (A) deterministic fields: do on all ranks
    for (size_t i = 0; i < NkTot_sz; ++i) {
        grid.iq[i] = kpatch.iq[i];
        grid.jq[i] = kpatch.jq[i];
        kvecF[i]   = kvecP[i];

        const unsigned char inside = (inHexPtr ? inHexPtr[i] : (unsigned char)1);
        inHexF[i]  = inside;

        EFfield[i] = EF;
        occAvgF[i] = NaN;

        evalsF[i].assign((size_t)dim, NaN);
        occbF[i].assign((size_t)dim, NaN);
    }

    // work partition
    const long long N = (long long)NkTot_sz;
    const long long base = N / (long long)nprocs;
    const long long rem  = N % (long long)nprocs;
    const long long my_n = base + (rank < rem ? 1 : 0);
    const long long my_s = (long long)rank * base + std::min<long long>(rank, rem);
    const long long my_e = my_s + my_n;

    std::vector<double> occ_local((size_t)my_n, NaN);
    std::vector<double> evals_local((size_t)my_n * (size_t)dim, NaN);
    std::vector<double> occb_local ((size_t)my_n * (size_t)dim, NaN);

    // local stats (in-hex only)
    long long n_inhex_local = 0;
    double occ_min_local = +std::numeric_limits<double>::infinity();
    double occ_max_local = -std::numeric_limits<double>::infinity();
    double occ_acc_local = 0.0;

    for (long long ii = my_s; ii < my_e; ++ii) {
        const size_t i   = (size_t)ii;
        const size_t loc = (size_t)(ii - my_s);

        if (!inHexF[i]) {
            occ_local[loc] = NaN;
            continue;
        }

        Eigen::MatrixXcd Hk = build_Hk(kvecP[i], true);
        Hk = 0.5 * (Hk + Hk.adjoint());

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(Hk);
        if (es.info() != Eigen::Success)
            rgmpi::abort_all("cal_fermi_patch_from_mu: eigensolver failed");

        const auto& ev = es.eigenvalues();

        double occ_sum = 0.0;
        for (int b = 0; b < dim; ++b) {
            const double fb = la::fermi(ev(b), EF, T_K);
            occ_sum += fb;

            evals_local[loc * (size_t)dim + (size_t)b] = ev(b);
            occb_local [loc * (size_t)dim + (size_t)b] = fb;
        }

        const double occ_avg = occ_sum / double(dim);
        occ_local[loc] = occ_avg;

        ++n_inhex_local;
        occ_acc_local += occ_avg;
        occ_min_local = std::min(occ_min_local, occ_avg);
        occ_max_local = std::max(occ_max_local, occ_avg);
    }

#ifdef USE_MPI
    if (rgmpi::inited() && nprocs > 1) {

        // counts/displs for Gatherv
        std::vector<int> counts(nprocs, 0), displs(nprocs, 0);
        for (int r = 0; r < nprocs; ++r) {
            const long long rn = base + (r < rem ? 1 : 0);
            const long long rs = (long long)r * base + std::min<long long>(r, rem);
            counts[r] = (int)rn;
            displs[r] = (int)rs;
        }

        // (C1) gather occ_k_avg to rank0 then Bcast back
        std::vector<double> occ_global;
        if (rank == 0) occ_global.assign((size_t)N, NaN);

        MPI_Gatherv(occ_local.data(), (int)my_n, MPI_DOUBLE,
                    (rank == 0 ? occ_global.data() : nullptr),
                    counts.data(), displs.data(), MPI_DOUBLE,
                    0, MPI_COMM_WORLD);

        if (rank == 0) {
            for (long long ii = 0; ii < N; ++ii) occAvgF[(size_t)ii] = occ_global[(size_t)ii];
        }
        MPI_Bcast(occAvgF.data(), (int)N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // (C2) gather per-band arrays (flat) then Bcast back
        std::vector<int> counts_dim(nprocs, 0), displs_dim(nprocs, 0);
        for (int r = 0; r < nprocs; ++r) {
            counts_dim[r] = counts[r] * dim;
            displs_dim[r] = displs[r] * dim;
        }

        std::vector<double> evals_global, occb_global;
        if (rank == 0) {
            evals_global.assign((size_t)N * (size_t)dim, NaN);
            occb_global .assign((size_t)N * (size_t)dim, NaN);
        }

        MPI_Gatherv(evals_local.data(), (int)(my_n * dim), MPI_DOUBLE,
                    (rank == 0 ? evals_global.data() : nullptr),
                    counts_dim.data(), displs_dim.data(), MPI_DOUBLE,
                    0, MPI_COMM_WORLD);

        MPI_Gatherv(occb_local.data(), (int)(my_n * dim), MPI_DOUBLE,
                    (rank == 0 ? occb_global.data() : nullptr),
                    counts_dim.data(), displs_dim.data(), MPI_DOUBLE,
                    0, MPI_COMM_WORLD);

        if (rank != 0) {
            evals_global.assign((size_t)N * (size_t)dim, NaN);
            occb_global .assign((size_t)N * (size_t)dim, NaN);
        }

        MPI_Bcast(evals_global.data(), (int)(N * dim), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(occb_global.data(),  (int)(N * dim), MPI_DOUBLE, 0, MPI_COMM_WORLD);

        for (long long ii = 0; ii < N; ++ii) {
            const size_t i = (size_t)ii;
            for (int b = 0; b < dim; ++b) {
                evalsF[i][(size_t)b] = evals_global[i*(size_t)dim + (size_t)b];
                occbF[i][(size_t)b]  = occb_global [i*(size_t)dim + (size_t)b];
            }
        }

        // summary stats
        long long n_inhex = 0;
        double occ_acc = 0.0;
        rgmpi::allreduce_sum(n_inhex_local, n_inhex);
        rgmpi::allreduce_sum(occ_acc_local, occ_acc);

        double min_in = (n_inhex_local > 0) ? occ_min_local : +std::numeric_limits<double>::infinity();
        double max_in = (n_inhex_local > 0) ? occ_max_local : -std::numeric_limits<double>::infinity();

        double occ_min = 0.0, occ_max = 0.0;
        rgmpi::allreduce_min(min_in, occ_min);
        rgmpi::allreduce_max(max_in, occ_max);

        if (debug_print_summary && rank == 0) {
            const double occ_mean = (n_inhex > 0 ? occ_acc / double(n_inhex) : NaN);
            std::cout << std::setprecision(16);
            std::cout << "\n=== cal_fermi_patch_from_mu: occ_k_avg summary (MPI) ===\n";
            std::cout << "  EF(input)=" << EF << " eV  T_K=" << T_K << " K\n";
            std::cout << "  n_inhex=" << n_inhex << "  dim(nband)=" << dim << "  ranks=" << nprocs << "\n";
            std::cout << "  occ_k_avg: min=" << occ_min
                      << "  max=" << occ_max
                      << "  mean=" << occ_mean << "\n";
            std::cout << "  (occ_k_avg(k) = (1/nband)*sum_b f(E_b(k),EF,T))\n";
            std::cout << "==================================================\n";
        }

        return grid;
    }
#endif

    // ---- serial path ----
    for (long long ii = my_s; ii < my_e; ++ii) {
        const size_t i   = (size_t)ii;
        const size_t loc = (size_t)(ii - my_s);

        occAvgF[i] = occ_local[loc];
        for (int b = 0; b < dim; ++b) {
            evalsF[i][(size_t)b] = evals_local[loc*(size_t)dim + (size_t)b];
            occbF[i][(size_t)b]  = occb_local [loc*(size_t)dim + (size_t)b];
        }
    }

    if (debug_print_summary) {
        long long n_inhex = 0;
        double occ_min = +std::numeric_limits<double>::infinity();
        double occ_max = -std::numeric_limits<double>::infinity();
        double occ_acc = 0.0;

        for (size_t i = 0; i < NkTot_sz; ++i) {
            if (!inHexF[i]) continue;
            const double x = occAvgF[i];
            ++n_inhex;
            occ_acc += x;
            occ_min = std::min(occ_min, x);
            occ_max = std::max(occ_max, x);
        }

        const double occ_mean = (n_inhex > 0 ? occ_acc / double(n_inhex) : NaN);
        std::cout << "\n=== cal_fermi_patch_from_mu finished ===\n";
        std::cout << "  EF(input)=" << EF << " eV  T_K=" << T_K << " K\n";
        std::cout << "  n_inhex=" << n_inhex << "  dim(nband)=" << dim << "\n";
        std::cout << "  occ_k_avg: min=" << occ_min
                  << "  max=" << occ_max
                  << "  mean=" << occ_mean << "\n";
    }

    return grid;
}

inline core::GridData RG_ModelBase::cal_fermi_patch_from_filling(double filling, double T_K, const core::GridData& kpatch) const
{
    if (!(T_K >= 0.0))
        rgmpi::abort_all("cal_fermi_patch_from_filling: T_K must be >= 0");

    const double EF = find_EF_from_filling(filling, T_K, kpatch);

    return cal_fermi_patch_from_mu(EF, T_K, kpatch);
}









} // namespace rg
