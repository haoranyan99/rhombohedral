#pragma once

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <vector>

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

    void   set_Dfield(double D) { Dfield_ = D; }
    double get_Dfield() const   { return Dfield_; }

    const RG_Structure& get_structure() const { return st_; }

    virtual Eigen::MatrixXcd build_Hk(
        const Vec2& k,
        bool enforce_hermitian = true
    ) const = 0;

    virtual Eigen::VectorXd bands_at_k(
        const Vec2& k,
        bool enforce_hermitian = true
    ) const;

    virtual Eigen::MatrixXd bands_along_path(
        const std::vector<Vec2>& klist,
        bool enforce_hermitian = true
    ) const;

    virtual Eigen::MatrixXd bands_along_path(
        const core::PathData& path,
        bool enforce_hermitian = true
    ) const;

    virtual core::SeriesData cal_dos_gaussian(
        const core::GridData& kpatch,
        double e_low,
        double e_high,
        int num_e,
        double eta,
        double T_K,
        bool enforce_hermitian = true
    ) const;

    virtual double find_EF_from_filling(
        double filling_target,
        double T_K,
        const core::GridData& kpatch
    ) const;

    virtual double find_EF_from_filling_T0(
        double filling_target,
        const core::GridData& kpatch
    ) const;

    virtual double find_filling_from_EF(
        double EF,
        double T_K,
        const core::GridData& kpatch
    ) const;

    virtual double find_filling_from_EF_T0(
        double EF,
        const core::GridData& kpatch
    ) const;

    virtual core::GridData cal_fermi_patch_from_mu(
        double EF,
        double T_K,
        const core::GridData& kpatch
    ) const;

    virtual core::GridData cal_chi_grid_Ef(
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
        double nest_thr = 50
    ) const;

    virtual void write_info(
        std::ostream& os,
        const std::string& prefix = "# "
    ) const
    {
        (void)os;
        (void)prefix;
    }

protected:
    static Eigen::MatrixXcd add_Dfield_(
        int norb,
        int Nl,
        double D
    );

    void write_info_common_(
        std::ostream& os,
        const std::string& prefix
    ) const;

protected:
    RG_Structure st_;
    double Dfield_ = 0.0;
};

// ============================================================
// small helpers
// ============================================================

inline Eigen::MatrixXcd RG_ModelBase::add_Dfield_(
    int norb,
    int Nl,
    double D
) {
    Eigen::MatrixXcd M =
        Eigen::MatrixXcd::Zero(norb, norb);

    std::vector<int> diag;

    if (Nl == 4) {
        diag = {-1,-1,0,0,1,1,2,2};
    } else if (Nl == 5) {
        diag = {-2,-2,-1,-1,0,0,1,1,2,2};
    } else if (Nl == 6) {
        diag = {-2,-2,-1,-1,0,0,1,1,2,2,3,3};
    } else {
        throw std::runtime_error(
            "RG_ModelBase::add_Dfield_: unsupported layer number"
        );
    }

    if (static_cast<int>(diag.size()) != norb) {
        throw std::runtime_error(
            "RG_ModelBase::add_Dfield_: diag size != norb"
        );
    }

    for (int i = 0; i < norb; ++i) {
        M(i, i) = static_cast<double>(diag[i]) * D;
    }

    return M;
}

inline Eigen::MatrixXcd hermitianize_(
    const Eigen::MatrixXcd& H
) {
    return 0.5 * (H + H.adjoint());
}

inline Eigen::VectorXd RG_ModelBase::bands_at_k(
    const Vec2& k,
    bool enforce_hermitian
) const {
    Eigen::MatrixXcd H = build_Hk(k, enforce_hermitian);

    if (enforce_hermitian) {
        H = hermitianize_(H);
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(H);

    if (es.info() != Eigen::Success) {
        throw std::runtime_error(
            "RG_ModelBase::bands_at_k: eigensolver failed"
        );
    }

    return es.eigenvalues();
}

inline Eigen::MatrixXd RG_ModelBase::bands_along_path(
    const std::vector<Vec2>& klist,
    bool enforce_hermitian
) const {
    if (klist.empty()) {
        return Eigen::MatrixXd();
    }

    const Eigen::VectorXd e0 =
        bands_at_k(klist.front(), enforce_hermitian);

    const int dim = static_cast<int>(e0.size());

    Eigen::MatrixXd E(
        static_cast<int>(klist.size()),
        dim
    );

    E.row(0) = e0.transpose();

    for (int i = 1; i < static_cast<int>(klist.size()); ++i) {
        E.row(i) =
            bands_at_k(klist[static_cast<size_t>(i)], enforce_hermitian)
            .transpose();
    }

    return E;
}

inline Eigen::MatrixXd RG_ModelBase::bands_along_path(
    const core::PathData& path,
    bool enforce_hermitian
) const {
    return bands_along_path(
        path.k_list,
        enforce_hermitian
    );
}

// ============================================================
// DOS
// ============================================================

inline core::SeriesData RG_ModelBase::cal_dos_gaussian(
    const core::GridData& kpatch,
    double e_low,
    double e_high,
    int num_e,
    double eta,
    double T_K,
    bool enforce_hermitian
) const {
    if (num_e < 2) {
        rgmpi::abort_all("cal_dos_gaussian: num_e must be >= 2");
    }

    if (!(e_high > e_low)) {
        rgmpi::abort_all("cal_dos_gaussian: e_high must be > e_low");
    }

    kpatch.assert_consistent();

    const size_t NkTot = kpatch.size();

    if (NkTot == 0) {
        rgmpi::abort_all("cal_dos_gaussian: empty kpatch");
    }

    int rank = 0;
    int nprocs = 1;
    rgmpi::rank_size(rank, nprocs);

    const auto& kvec = kpatch.get<Vec2>("kvec").v;

    const int dim =
        static_cast<int>(st_.atoms().size());

    if (dim <= 0) {
        rgmpi::abort_all("cal_dos_gaussian: dim <= 0");
    }

    const double dE =
        (e_high - e_low)
        / static_cast<double>(num_e - 1);

    if (!(eta > 0.0)) {
        eta = 0.6 * dE;
    }

    if (!(eta > 0.0)) {
        rgmpi::abort_all("cal_dos_gaussian: eta must be > 0");
    }

    core::SeriesData out;
    out.resize(static_cast<size_t>(num_e));
    out.dx = dE;

    auto& dosF =
        out.add<double>("dos").v;

    auto& fillingF =
        out.add<double>("filling").v;

    auto& dopingF =
        out.add<double>("doping").v;

    for (int ie = 0; ie < num_e; ++ie) {
        out.x[static_cast<size_t>(ie)] =
            e_low + dE * static_cast<double>(ie);

        dosF[static_cast<size_t>(ie)] = 0.0;
        fillingF[static_cast<size_t>(ie)] = 0.0;
        dopingF[static_cast<size_t>(ie)] = 0.0;
    }

    std::vector<double> dos_local(
        static_cast<size_t>(num_e),
        0.0
    );

    std::vector<double> filling_local(
        static_cast<size_t>(num_e),
        0.0
    );

    const auto [i0, i1] =
        rgmpi::block_1d_int(
            static_cast<int>(NkTot),
            rank,
            nprocs
        );

    const double inv_sqrt2pi =
        1.0 / std::sqrt(2.0 * la::pi);

    const double pref =
        inv_sqrt2pi / eta;

    const double inv2eta2 =
        1.0 / (2.0 * eta * eta);

    for (int ii = i0; ii < i1; ++ii) {
        const size_t ik =
            static_cast<size_t>(ii);

        Eigen::MatrixXcd Hk =
            build_Hk(kvec[ik], enforce_hermitian);

        if (enforce_hermitian) {
            Hk = hermitianize_(Hk);
        }

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(Hk);

        if (es.info() != Eigen::Success) {
            rgmpi::abort_all(
                "cal_dos_gaussian: eigensolver failed"
            );
        }

        const auto& ev = es.eigenvalues();

        for (int ie = 0; ie < num_e; ++ie) {
            const double E =
                out.x[static_cast<size_t>(ie)];

            double acc_dos = 0.0;
            double acc_fill = 0.0;

            for (int b = 0; b < dim; ++b) {
                const double eb = ev(b);
                const double de = E - eb;

                acc_dos +=
                    std::exp(-(de * de) * inv2eta2)
                    * pref;

                acc_fill +=
                    la::fermi(eb, E, T_K);
            }

            dos_local[static_cast<size_t>(ie)] += acc_dos;
            filling_local[static_cast<size_t>(ie)] += acc_fill;
        }
    }

    rgmpi::allreduce_sum_vector(dos_local);
    rgmpi::allreduce_sum_vector(filling_local);

    const double invNk =
        1.0 / static_cast<double>(NkTot);

    const double invNkDim =
        1.0
        / (
            static_cast<double>(NkTot)
          * static_cast<double>(dim)
        );

    const double a_Ang = st_.a();

    for (int ie = 0; ie < num_e; ++ie) {
        const size_t i =
            static_cast<size_t>(ie);

        dosF[i] =
            dos_local[i] * invNk;

        fillingF[i] =
            filling_local[i] * invNkDim;

        dopingF[i] =
            la::filling_to_doping(
                fillingF[i],
                a_Ang
            );
    }

    out.assert_consistent();
    return out;
}

// ============================================================
// filling from EF
// ============================================================

inline double RG_ModelBase::find_filling_from_EF(
    double EF,
    double T_K,
    const core::GridData& kpatch
) const {
    if (T_K == 0.0) {
        return find_filling_from_EF_T0(EF, kpatch);
    }

    kpatch.assert_consistent();

    const size_t NkTot = kpatch.size();

    if (NkTot == 0) {
        rgmpi::abort_all("find_filling_from_EF: empty kpatch");
    }

    if (!(T_K > 0.0)) {
        rgmpi::abort_all("find_filling_from_EF: T_K must be > 0");
    }

    int rank = 0;
    int nprocs = 1;
    rgmpi::rank_size(rank, nprocs);

    const auto& kvec =
        kpatch.get<Vec2>("kvec").v;

    const int dim =
        static_cast<int>(st_.atoms().size());

    const auto [i0, i1] =
        rgmpi::block_1d_int(
            static_cast<int>(NkTot),
            rank,
            nprocs
        );

    double occ_local = 0.0;

    for (int ii = i0; ii < i1; ++ii) {
        const size_t ik =
            static_cast<size_t>(ii);

        Eigen::MatrixXcd Hk =
            hermitianize_(build_Hk(kvec[ik], true));

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(Hk);

        if (es.info() != Eigen::Success) {
            rgmpi::abort_all(
                "find_filling_from_EF: eigensolver failed"
            );
        }

        const auto& ev = es.eigenvalues();

        for (int b = 0; b < dim; ++b) {
            occ_local += la::fermi(ev(b), EF, T_K);
        }
    }

    double occ_sum = 0.0;
    rgmpi::allreduce_sum(occ_local, occ_sum);

    return
        occ_sum
        / (
            static_cast<double>(NkTot)
          * static_cast<double>(dim)
        );
}

inline double RG_ModelBase::find_filling_from_EF_T0(
    double EF,
    const core::GridData& kpatch
) const {
    kpatch.assert_consistent();

    const size_t NkTot = kpatch.size();

    if (NkTot == 0) {
        rgmpi::abort_all("find_filling_from_EF_T0: empty kpatch");
    }

    int rank = 0;
    int nprocs = 1;
    rgmpi::rank_size(rank, nprocs);

    const auto& kvec =
        kpatch.get<Vec2>("kvec").v;

    const int dim =
        static_cast<int>(st_.atoms().size());

    const auto [i0, i1] =
        rgmpi::block_1d_int(
            static_cast<int>(NkTot),
            rank,
            nprocs
        );

    long long count_local = 0;

    for (int ii = i0; ii < i1; ++ii) {
        const size_t ik =
            static_cast<size_t>(ii);

        Eigen::MatrixXcd Hk =
            hermitianize_(build_Hk(kvec[ik], true));

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(Hk);

        if (es.info() != Eigen::Success) {
            rgmpi::abort_all(
                "find_filling_from_EF_T0: eigensolver failed"
            );
        }

        const auto& ev = es.eigenvalues();

        for (int b = 0; b < dim; ++b) {
            if (ev(b) < EF) {
                ++count_local;
            }
        }
    }

    long long count_sum = 0;
    rgmpi::allreduce_sum(count_local, count_sum);

    return
        static_cast<double>(count_sum)
        / (
            static_cast<double>(NkTot)
          * static_cast<double>(dim)
        );
}

// ============================================================
// EF from filling
// ============================================================

inline double RG_ModelBase::find_EF_from_filling(
    double filling_target,
    double T_K,
    const core::GridData& kpatch
) const {
    if (T_K == 0.0) {
        return find_EF_from_filling_T0(
            filling_target,
            kpatch
        );
    }

    if (!(filling_target >= 0.0 && filling_target <= 1.0)) {
        rgmpi::abort_all(
            "find_EF_from_filling: filling_target must be in [0,1]"
        );
    }

    kpatch.assert_consistent();

    const size_t NkTot = kpatch.size();

    if (NkTot == 0) {
        rgmpi::abort_all("find_EF_from_filling: empty kpatch");
    }

    int rank = 0;
    int nprocs = 1;
    rgmpi::rank_size(rank, nprocs);

    const auto& kvec =
        kpatch.get<Vec2>("kvec").v;

    const int dim =
        static_cast<int>(st_.atoms().size());

    const auto [i0, i1] =
        rgmpi::block_1d_int(
            static_cast<int>(NkTot),
            rank,
            nprocs
        );

    std::vector<double> evals_local;
    evals_local.reserve(
        static_cast<size_t>(i1 - i0)
      * static_cast<size_t>(dim)
    );

    double Emin_local =
        +std::numeric_limits<double>::infinity();

    double Emax_local =
        -std::numeric_limits<double>::infinity();

    for (int ii = i0; ii < i1; ++ii) {
        const size_t ik =
            static_cast<size_t>(ii);

        Eigen::MatrixXcd Hk =
            hermitianize_(build_Hk(kvec[ik], true));

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(Hk);

        if (es.info() != Eigen::Success) {
            rgmpi::abort_all(
                "find_EF_from_filling: eigensolver failed"
            );
        }

        const auto& ev = es.eigenvalues();

        for (int b = 0; b < dim; ++b) {
            const double e = ev(b);

            evals_local.push_back(e);

            Emin_local = std::min(Emin_local, e);
            Emax_local = std::max(Emax_local, e);
        }
    }

    double Emin = 0.0;
    double Emax = 0.0;

    rgmpi::allreduce_min(Emin_local, Emin);
    rgmpi::allreduce_max(Emax_local, Emax);

    const double norm =
        1.0
        / (
            static_cast<double>(NkTot)
          * static_cast<double>(dim)
        );

    auto filling_of_EF =
        [&](double EF) -> double {
            double acc_local = 0.0;

            for (double e : evals_local) {
                acc_local += la::fermi(e, EF, T_K);
            }

            double acc = 0.0;
            rgmpi::allreduce_sum(acc_local, acc);

            return acc * norm;
        };

    double lo = Emin - 1e-6;
    double hi = Emax + 1e-6;

    double f_lo =
        filling_of_EF(lo) - filling_target;

    double f_hi =
        filling_of_EF(hi) - filling_target;

    if (f_lo * f_hi > 0.0) {
        rgmpi::abort_all(
            "find_EF_from_filling: failed to bracket root"
        );
    }

    constexpr int max_iter = 100;
    constexpr double tol_EF = 1e-10;
    constexpr double tol_f = 1e-10;

    for (int it = 0; it < max_iter; ++it) {
        const double mid =
            0.5 * (lo + hi);

        const double f_mid =
            filling_of_EF(mid) - filling_target;

        if (std::abs(f_mid) < tol_f) {
            return mid;
        }

        if ((hi - lo) < tol_EF) {
            return mid;
        }

        if (f_lo * f_mid <= 0.0) {
            hi = mid;
            f_hi = f_mid;
        } else {
            lo = mid;
            f_lo = f_mid;
        }
    }

    return 0.5 * (lo + hi);
}

inline double RG_ModelBase::find_EF_from_filling_T0(
    double filling_target,
    const core::GridData& kpatch
) const {
    if (!(filling_target >= 0.0 && filling_target <= 1.0)) {
        rgmpi::abort_all(
            "find_EF_from_filling_T0: filling_target must be in [0,1]"
        );
    }

    kpatch.assert_consistent();

    const size_t NkTot = kpatch.size();

    if (NkTot == 0) {
        rgmpi::abort_all("find_EF_from_filling_T0: empty kpatch");
    }

    int rank = 0;
    int nprocs = 1;
    rgmpi::rank_size(rank, nprocs);

    const auto& kvec =
        kpatch.get<Vec2>("kvec").v;

    const int dim =
        static_cast<int>(st_.atoms().size());

    const auto [i0, i1] =
        rgmpi::block_1d_int(
            static_cast<int>(NkTot),
            rank,
            nprocs
        );

    std::vector<double> evals_local;
    evals_local.reserve(
        static_cast<size_t>(i1 - i0)
      * static_cast<size_t>(dim)
    );

    double Emin_local =
        +std::numeric_limits<double>::infinity();

    double Emax_local =
        -std::numeric_limits<double>::infinity();

    for (int ii = i0; ii < i1; ++ii) {
        const size_t ik =
            static_cast<size_t>(ii);

        Eigen::MatrixXcd Hk =
            hermitianize_(build_Hk(kvec[ik], true));

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(Hk);

        if (es.info() != Eigen::Success) {
            rgmpi::abort_all(
                "find_EF_from_filling_T0: eigensolver failed"
            );
        }

        const auto& ev = es.eigenvalues();

        for (int b = 0; b < dim; ++b) {
            const double e = ev(b);

            evals_local.push_back(e);

            Emin_local = std::min(Emin_local, e);
            Emax_local = std::max(Emax_local, e);
        }
    }

    double Emin = 0.0;
    double Emax = 0.0;

    rgmpi::allreduce_min(Emin_local, Emin);
    rgmpi::allreduce_max(Emax_local, Emax);

    if (filling_target <= 0.0) {
        return Emin - 1e-6;
    }

    if (filling_target >= 1.0) {
        return Emax + 1e-6;
    }

    const double norm =
        1.0
        / (
            static_cast<double>(NkTot)
          * static_cast<double>(dim)
        );

    auto filling_of_EF =
        [&](double EF) -> double {
            long long cnt_local = 0;

            for (double e : evals_local) {
                if (e < EF) {
                    ++cnt_local;
                }
            }

            long long cnt = 0;
            rgmpi::allreduce_sum(cnt_local, cnt);

            return static_cast<double>(cnt) * norm;
        };

    double lo = Emin - 1e-6;
    double hi = Emax + 1e-6;

    constexpr int max_iter = 100;
    constexpr double tol_EF = 1e-10;

    for (int it = 0; it < max_iter; ++it) {
        const double mid =
            0.5 * (lo + hi);

        const double f_mid =
            filling_of_EF(mid) - filling_target;

        if ((hi - lo) < tol_EF) {
            return mid;
        }

        if (f_mid >= 0.0) {
            hi = mid;
        } else {
            lo = mid;
        }
    }

    return 0.5 * (lo + hi);
}

// ============================================================
// fermi patch with eigenvectors
// ============================================================

inline core::GridData RG_ModelBase::cal_fermi_patch_from_mu(
    double EF,
    double T_K,
    const core::GridData& kpatch
) const {
    kpatch.assert_consistent();

    const size_t NkTot = kpatch.size();

    if (NkTot == 0) {
        rgmpi::abort_all("cal_fermi_patch_from_mu: empty kpatch");
    }

    if (!(T_K >= 0.0)) {
        rgmpi::abort_all("cal_fermi_patch_from_mu: T_K must be >= 0");
    }

    if (!std::isfinite(EF)) {
        rgmpi::abort_all("cal_fermi_patch_from_mu: EF must be finite");
    }

    int rank = 0;
    int nprocs = 1;
    rgmpi::rank_size(rank, nprocs);

    const auto& kvecP =
        kpatch.get<Vec2>("kvec").v;

    Eigen::MatrixXcd H0 =
        hermitianize_(build_Hk(kvecP[0], true));

    const int dim =
        static_cast<int>(H0.rows());

    if (H0.cols() != dim) {
        rgmpi::abort_all(
            "cal_fermi_patch_from_mu: H(k) is not square"
        );
    }

    core::GridData out;
    out.resize(NkTot);

    out.dx = kpatch.dx;
    out.dy = kpatch.dy;
    out.mesh_type = kpatch.mesh_type;

    auto& kvecF =
        out.add<Vec2>("kvec").v;

    auto& EFfield =
        out.add<double>("EF_used").v;

    auto& occAvgF =
        out.add<double>("occ_k_avg").v;

    auto& evalsF =
        out.add<std::vector<double>>("evals").v;

    auto& occbF =
        out.add<std::vector<double>>("occ_band").v;

    auto& evecReF =
        out.add<std::vector<double>>("evec_re").v;

    auto& evecImF =
        out.add<std::vector<double>>("evec_im").v;

    const double NaN =
        std::numeric_limits<double>::quiet_NaN();

    for (size_t i = 0; i < NkTot; ++i) {
        out.iq[i] = kpatch.iq[i];
        out.jq[i] = kpatch.jq[i];

        kvecF[i] = kvecP[i];

        EFfield[i] = EF;
        occAvgF[i] = NaN;

        evalsF[i].assign(static_cast<size_t>(dim), NaN);
        occbF[i].assign(static_cast<size_t>(dim), NaN);

        evecReF[i].assign(
            static_cast<size_t>(dim) * static_cast<size_t>(dim),
            NaN
        );

        evecImF[i].assign(
            static_cast<size_t>(dim) * static_cast<size_t>(dim),
            NaN
        );
    }

    const auto [i0, i1] =
        rgmpi::block_1d_int(
            static_cast<int>(NkTot),
            rank,
            nprocs
        );

    std::vector<double> occ_all(NkTot, 0.0);

    std::vector<double> evals_all(
        NkTot * static_cast<size_t>(dim),
        0.0
    );

    std::vector<double> occb_all(
        NkTot * static_cast<size_t>(dim),
        0.0
    );

    std::vector<double> evecRe_all(
        NkTot
      * static_cast<size_t>(dim)
      * static_cast<size_t>(dim),
        0.0
    );

    std::vector<double> evecIm_all(
        NkTot
      * static_cast<size_t>(dim)
      * static_cast<size_t>(dim),
        0.0
    );

    for (int ii = i0; ii < i1; ++ii) {
        const size_t ik =
            static_cast<size_t>(ii);

        Eigen::MatrixXcd Hk =
            hermitianize_(build_Hk(kvecP[ik], true));

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(Hk);

        if (es.info() != Eigen::Success) {
            rgmpi::abort_all(
                "cal_fermi_patch_from_mu: eigensolver failed"
            );
        }

        const auto& ev =
            es.eigenvalues();

        const auto& U =
            es.eigenvectors();

        double occ_sum = 0.0;

        for (int b = 0; b < dim; ++b) {
            const double fb =
                la::fermi(ev(b), EF, T_K);

            occ_sum += fb;

            evals_all[
                ik * static_cast<size_t>(dim)
              + static_cast<size_t>(b)
            ] = ev(b);

            occb_all[
                ik * static_cast<size_t>(dim)
              + static_cast<size_t>(b)
            ] = fb;
        }

        occ_all[ik] =
            occ_sum / static_cast<double>(dim);

        for (int a = 0; a < dim; ++a) {
            for (int b = 0; b < dim; ++b) {
                const size_t pos =
                    ik * static_cast<size_t>(dim) * static_cast<size_t>(dim)
                  + static_cast<size_t>(a) * static_cast<size_t>(dim)
                  + static_cast<size_t>(b);

                evecRe_all[pos] = U(a, b).real();
                evecIm_all[pos] = U(a, b).imag();
            }
        }
    }

    rgmpi::allreduce_sum_vector(occ_all);
    rgmpi::allreduce_sum_vector(evals_all);
    rgmpi::allreduce_sum_vector(occb_all);
    rgmpi::allreduce_sum_vector(evecRe_all);
    rgmpi::allreduce_sum_vector(evecIm_all);

    for (size_t ik = 0; ik < NkTot; ++ik) {
        occAvgF[ik] = occ_all[ik];

        for (int b = 0; b < dim; ++b) {
            evalsF[ik][static_cast<size_t>(b)] =
                evals_all[
                    ik * static_cast<size_t>(dim)
                  + static_cast<size_t>(b)
                ];

            occbF[ik][static_cast<size_t>(b)] =
                occb_all[
                    ik * static_cast<size_t>(dim)
                  + static_cast<size_t>(b)
                ];
        }

        for (int a = 0; a < dim; ++a) {
            for (int b = 0; b < dim; ++b) {
                const size_t global_pos =
                    ik * static_cast<size_t>(dim) * static_cast<size_t>(dim)
                  + static_cast<size_t>(a) * static_cast<size_t>(dim)
                  + static_cast<size_t>(b);

                const size_t local_pos =
                    static_cast<size_t>(a) * static_cast<size_t>(dim)
                  + static_cast<size_t>(b);

                evecReF[ik][local_pos] =
                    evecRe_all[global_pos];

                evecImF[ik][local_pos] =
                    evecIm_all[global_pos];
            }
        }
    }

    out.assert_consistent();
    return out;
}

// ============================================================
// direct chi from model
// NOTE:
// This is kept only for quick direct calculations.
// For production bin fermiPatch chi, use cal_susceptibility.h.
// ============================================================

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
    double nest_thr
) const {
    (void)nest_thr;

    if (Nq < 0) {
        rgmpi::abort_all("cal_chi_grid_Ef: Nq must be >= 0");
    }

    if (!(eta > 0.0)) {
        rgmpi::abort_all("cal_chi_grid_Ef: eta must be > 0");
    }

    if (!(T_K >= 0.0)) {
        rgmpi::abort_all("cal_chi_grid_Ef: T_K must be >= 0");
    }

    kpatch.assert_consistent();

    const size_t NkTot = kpatch.size();

    if (NkTot == 0) {
        rgmpi::abort_all("cal_chi_grid_Ef: empty kpatch");
    }

    int rank = 0;
    int nprocs = 1;
    rgmpi::rank_size(rank, nprocs);

    const auto& kvec =
        kpatch.get<Vec2>("kvec").v;

    const int dim =
        static_cast<int>(st_.atoms().size());

    const bool pm0 =
        std::abs(polar_mu) < 1e-15;

    const auto [i0, i1] =
        rgmpi::block_1d_int(
            static_cast<int>(NkTot),
            rank,
            nprocs
        );

    std::vector<double> evals_all(
        NkTot * static_cast<size_t>(dim),
        0.0
    );

    std::vector<double> occ_up_all(
        NkTot * static_cast<size_t>(dim),
        0.0
    );

    std::vector<double> occ_dn_all(
        NkTot * static_cast<size_t>(dim),
        0.0
    );

    for (int ii = i0; ii < i1; ++ii) {
        const size_t ik =
            static_cast<size_t>(ii);

        Eigen::MatrixXcd Hk =
            build_Hk(kvec[ik], enforce_hermitian);

        if (enforce_hermitian) {
            Hk = hermitianize_(Hk);
        }

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(Hk);

        if (es.info() != Eigen::Success) {
            rgmpi::abort_all(
                "cal_chi_grid_Ef: eigensolver failed"
            );
        }

        const auto& ev =
            es.eigenvalues();

        for (int b = 0; b < dim; ++b) {
            const size_t pos =
                ik * static_cast<size_t>(dim)
              + static_cast<size_t>(b);

            const double E = ev(b);

            evals_all[pos] = E;

            if (pm0) {
                const double f =
                    la::fermi(E, EF, T_K);

                occ_up_all[pos] = f;
                occ_dn_all[pos] = f;
            } else {
                occ_up_all[pos] =
                    la::fermi(E - 0.5 * polar_mu, EF, T_K);

                occ_dn_all[pos] =
                    la::fermi(E + 0.5 * polar_mu, EF, T_K);
            }
        }
    }

    rgmpi::allreduce_sum_vector(evals_all);
    rgmpi::allreduce_sum_vector(occ_up_all);
    rgmpi::allreduce_sum_vector(occ_dn_all);

    const int qSide =
        2 * Nq + 1;

    const size_t Nout =
        static_cast<size_t>(qSide)
      * static_cast<size_t>(qSide);

    core::GridData qgrid;
    qgrid.resize(Nout);
    qgrid.dx = kpatch.dx;
    qgrid.dy = kpatch.dy;
    qgrid.mesh_type = "q_square";

    auto& chiF =
        qgrid.add<cd>("chi").v;

    auto& nPairF =
        qgrid.add<long long>("nKpair").v;

    auto& qxF =
        qgrid.add<double>("qx").v;

    auto& qyF =
        qgrid.add<double>("qy").v;

    const Vec2 dq1 =
        kpatch.dx * st_.b1();

    const Vec2 dq2 =
        kpatch.dy * st_.b2();

    size_t out_idx = 0;

    for (int djq = -Nq; djq <= Nq; ++djq) {
        for (int diq = -Nq; diq <= Nq; ++diq, ++out_idx) {
            const int iq =
                q_center_iq + diq;

            const int jq =
                q_center_jq + djq;

            qgrid.iq[out_idx] = iq;
            qgrid.jq[out_idx] = jq;

            const Vec2 qvec =
                static_cast<double>(iq) * dq1
              + static_cast<double>(jq) * dq2;

            qxF[out_idx] = qvec.x();
            qyF[out_idx] = qvec.y();

            chiF[out_idx] = cd(0.0, 0.0);
            nPairF[out_idx] = 0;
        }
    }

    std::vector<double> chi_re(Nout, 0.0);
    std::vector<double> chi_im(Nout, 0.0);
    std::vector<long long> nPair(Nout, 0);

    auto wrap_existing =
        [&](int iq, int jq, size_t& idx) -> bool {
            if (!boundary_periodic) {
                return kpatch.ij_to_idx(iq, jq, idx);
            }

            if (
                kpatch.mesh_type == "hex"
             || kpatch.mesh_type == "unknown"
            ) {
                return kpatch.ij_to_idx(iq, jq, idx);
            }

            int iq_min = 1000000000;
            int iq_max = -1000000000;
            int jq_min = 1000000000;
            int jq_max = -1000000000;

            for (size_t t = 0; t < NkTot; ++t) {
                iq_min = std::min(iq_min, kpatch.iq[t]);
                iq_max = std::max(iq_max, kpatch.iq[t]);
                jq_min = std::min(jq_min, kpatch.jq[t]);
                jq_max = std::max(jq_max, kpatch.jq[t]);
            }

            const int N1 = iq_max - iq_min + 1;
            const int N2 = jq_max - jq_min + 1;

            const int wiq =
                iq_min + la::mod_pos_int(iq - iq_min, N1);

            const int wjq =
                jq_min + la::mod_pos_int(jq - jq_min, N2);

            return kpatch.ij_to_idx(wiq, wjq, idx);
        };

    for (size_t out = 0; out < Nout; ++out) {
        const int dq_iq =
            qgrid.iq[out];

        const int dq_jq =
            qgrid.jq[out];

        cd chi_local(0.0, 0.0);
        long long nPair_local = 0;

        for (int ii = i0; ii < i1; ++ii) {
            const size_t ik =
                static_cast<size_t>(ii);

            const int iq1 =
                kpatch.iq[ik];

            const int jq1 =
                kpatch.jq[ik];

            size_t ik2 = 0;

            if (!wrap_existing(
                    iq1 + dq_iq,
                    jq1 + dq_jq,
                    ik2
                ))
            {
                continue;
            }

            for (int b = 0; b < dim; ++b) {
                const size_t p1 =
                    ik * static_cast<size_t>(dim)
                  + static_cast<size_t>(b);

                const double E1 =
                    evals_all[p1];

                if (pm0) {
                    const double f1 =
                        occ_up_all[p1];

                    for (int m = 0; m < dim; ++m) {
                        const size_t p2 =
                            ik2 * static_cast<size_t>(dim)
                          + static_cast<size_t>(m);

                        const double E2 =
                            evals_all[p2];

                        const double f2 =
                            occ_up_all[p2];

                        chi_local +=
                            (f2 - f1)
                            / cd(E1 - E2, eta);

                        ++nPair_local;
                    }
                } else {
                    const double f1u =
                        occ_up_all[p1];

                    const double f1d =
                        occ_dn_all[p1];

                    const double E1u =
                        E1 - 0.5 * polar_mu;

                    const double E1d =
                        E1 + 0.5 * polar_mu;

                    for (int m = 0; m < dim; ++m) {
                        const size_t p2 =
                            ik2 * static_cast<size_t>(dim)
                          + static_cast<size_t>(m);

                        const double E2 =
                            evals_all[p2];

                        const double f2u =
                            occ_up_all[p2];

                        const double f2d =
                            occ_dn_all[p2];

                        const double E2u =
                            E2 - 0.5 * polar_mu;

                        const double E2d =
                            E2 + 0.5 * polar_mu;

                        chi_local +=
                            0.25 * (f2u - f1u)
                            / cd(E1u - E2u, eta);

                        chi_local +=
                            0.25 * (f2d - f1d)
                            / cd(E1d - E2d, eta);

                        chi_local +=
                            0.25 * (f2u - f1d)
                            / cd(E1d - E2u, eta);

                        chi_local +=
                            0.25 * (f2d - f1u)
                            / cd(E1u - E2d, eta);

                        nPair_local += 4;
                    }
                }
            }
        }

        chi_re[out] = chi_local.real();
        chi_im[out] = chi_local.imag();
        nPair[out] = nPair_local;
    }

    rgmpi::allreduce_sum_vector(chi_re);
    rgmpi::allreduce_sum_vector(chi_im);
    rgmpi::allreduce_sum_vector(nPair);

    const double norm =
        1.0 / static_cast<double>(NkTot);

    for (size_t out = 0; out < Nout; ++out) {
        chiF[out] =
            norm * cd(chi_re[out], chi_im[out]);

        nPairF[out] = nPair[out];
    }

    qgrid.assert_consistent();
    return qgrid;
}

} // namespace rg