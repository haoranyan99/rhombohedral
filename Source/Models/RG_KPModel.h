// File: Source/Models/RG_KPModel.h
#pragma once

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cmath>
#include <complex>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <utility>

#include "Common/DataContainers.h"
#include "PhysStruct/RG_Structure.h"
#include "LinearAlgebra/Constants.h"
#include "LinearAlgebra/MathFunctions.h"
#include "Source/Models/RG_ModelBase.h"

namespace rg {

class RG_KPModel : public RG_ModelBase {
public:
    using Vec2 = Eigen::Vector2d;
    using Vec3 = Eigen::Vector3d;
    using cd   = std::complex<double>;

    using RG_ModelBase::get_structure;
    using RG_ModelBase::get_Dfield;
    using RG_ModelBase::set_Dfield;
    using RG_ModelBase::add_Dfield_;
    using RG_ModelBase::bands_at_k;
    using RG_ModelBase::bands_along_path;
    using RG_ModelBase::cal_dos_gaussian;
    using RG_ModelBase::cal_chi_grid_Ef;
    using RG_ModelBase::find_EF_from_filling;
    using RG_ModelBase::cal_fermi_patch_from_filling;

    // SWMcC parameters (eV)
    struct RGKPParams {
        double gamma0 = 3.16;
        double gamma1 = 0.502;
        double gamma3 = -0.377;
        double gamma4 = -0.099;
        double gamma2 = -0.0171;
        double gamma5 = 0.0;
        double gamma6 = 0.0;

        double epsilon = -0.0014;

        bool include_gamma2 = true;
        bool include_gamma3 = true;
        bool include_gamma4 = true;

        // (optional but small) search radius around rounded (m0,n0)
        // minimal safe default:
        int G_neighbor = 1;

        double t_min = 1e-12;
    };


public:
    RG_KPModel() = default;
    explicit RG_KPModel(const RG_Structure& st) : RG_ModelBase(st) {}

    RGKPParams&       get_params()       { return p_; }
    const RGKPParams& get_params() const { return p_; }

    Eigen::MatrixXcd build_Hk(const Vec2& k, bool enforce_hermitian = true) const;

    void write_info(std::ostream& os, const std::string& prefix = "# ") const override;


private:
    static int idxA_(int layer) { return 2 * layer + 0; }
    static int idxB_(int layer) { return 2 * layer + 1; }

    // Return (q, xi) where q = k - (K_or_Kp + G) is the *closest* valley-relative momentum.
    // xi = +1 for K, xi = -1 for Kp (chirality label used in π = xi*q_x + i q_y).
    std::pair<Vec2,int> nearest_valley_q_(const Vec2& k) const;

    // π = xi*qx + i qy ;  π† = xi*qx - i qy
    static cd pi_(const Vec2& q, int xi)     { return cd(double(xi) * q.x(),  q.y()); }
    static cd pi_dag_(const Vec2& q, int xi) { return cd(double(xi) * q.x(), -q.y()); }

    // Convert gamma_i to v_i coefficient in linearized k·p (keeping your coordinate units).
    // v ~ (sqrt(3)/2)*a*gamma  (this matches standard linearization scale).
    double v_from_gamma_(double gamma) const {
        const double a = st_.a();  // your in-plane lattice constant
        return (std::sqrt(3.0) * 0.5) * a * gamma;
    }


private:
    RGKPParams p_;
};

inline std::pair<Vec2,int> RG_KPModel::nearest_valley_q_(const Vec2& k) const
{
    const Vec2 b1 = st_.b1();
    const Vec2 b2 = st_.b2();
    const Vec2 K  = st_.K();
    const Vec2 Kp = st_.Kp();

    Eigen::Matrix2d B;
    B.col(0) = b1;
    B.col(1) = b2;

    const double det = B.determinant();
    if (std::abs(det) < 1e-14)
        throw std::runtime_error("RG_KPModel::nearest_valley_q_: b1,b2 nearly singular");
    const Eigen::Matrix2d Binv = B.inverse();

    auto best_q_for_center = [&](const Vec2& center)->std::pair<Vec2,double> {
        const Eigen::Vector2d uv = Binv * (k - center);
        const int m0 = (int)std::llround(uv.x());
        const int n0 = (int)std::llround(uv.y());

        const int R = std::max(1, p_.G_neighbor); // or hardcode R=1 for "ultra-minimal"
        Vec2 best_q(0,0);
        double best_n2 = std::numeric_limits<double>::infinity();

        for (int dm = -R; dm <= R; ++dm) {
            for (int dn = -R; dn <= R; ++dn) {
                const int m = m0 + dm;
                const int n = n0 + dn;
                const Vec2 G = double(m)*b1 + double(n)*b2;
                const Vec2 q = k - (center + G);
                const double n2 = q.squaredNorm();
                if (n2 < best_n2) { best_n2 = n2; best_q = q; }
            }
        }
        return {best_q, best_n2};
    };

    auto [qK,  dK2 ] = best_q_for_center(K);
    auto [qKp, dKp2] = best_q_for_center(Kp);

    // xi: +1 for K, -1 for Kp (chirality label)
    if (dK2 <= dKp2) return {qK,  +1};
    else             return {qKp, -1};
}


inline Eigen::MatrixXcd RG_KPModel::build_Hk(const Vec2& k, bool enforce_hermitian) const
{
    const int Nl   = st_.nLayer();
    const int norb = 2 * Nl;
    if ((int)st_.atoms().size() != norb)
        throw std::runtime_error("RG_KPModel::build_Hk: orbit count != 2*Nl");

    const auto [q, xi] = nearest_valley_q_(k);

    const cd  pi    = cd(double(xi)*q.x(),  q.y());
    const cd  pidag = cd(double(xi)*q.x(), -q.y());

    const double v0 = (std::sqrt(3.0)*0.5) * st_.a() * p_.gamma0;
    const double v3 = (std::sqrt(3.0)*0.5) * st_.a() * p_.gamma3;
    const double v4 = (std::sqrt(3.0)*0.5) * st_.a() * p_.gamma4;

    Eigen::MatrixXcd H = Eigen::MatrixXcd::Zero(norb, norb);

    // intralayer Dirac
    for (int l = 0; l < Nl; ++l) {
        const int ia = idxA_(l), ib = idxB_(l);
        H(ia, ib) += v0 * pidag;
        H(ib, ia) += v0 * pi;
    }

    // interlayer ABC couplings
    for (int l = 0; l < Nl - 1; ++l) {
        const int A_l  = idxA_(l);
        const int B_l  = idxB_(l);
        const int A_lp = idxA_(l + 1);
        const int B_lp = idxB_(l + 1);

        H(B_l,  A_lp) += p_.gamma1;
        H(A_lp, B_l ) += p_.gamma1;

        if (p_.include_gamma3 && std::abs(p_.gamma3) > p_.t_min) {
            H(A_l,  B_lp) += v3 * pi;
            H(B_lp, A_l ) += v3 * pidag;
        }

        if (p_.include_gamma4 && std::abs(p_.gamma4) > p_.t_min) {
            H(A_l,  A_lp) += v4 * pidag;
            H(A_lp, A_l ) += v4 * pi;

            H(B_l,  B_lp) += v4 * pidag;
            H(B_lp, B_l ) += v4 * pi;
        }
    }

    // gamma2: A1 <-> B_N
    if (p_.include_gamma2 && std::abs(p_.gamma2) > p_.t_min) {
        const int A1 = idxA_(0);
        const int BN = idxB_(Nl - 1);
        H(A1, BN) += p_.gamma2;
        H(BN, A1) += p_.gamma2;
    }

    if (std::abs(p_.epsilon) > 0.0) {
        H(idxA_(0),      idxA_(0))     += p_.epsilon;
        H(idxB_(Nl - 1), idxB_(Nl - 1))+= p_.epsilon;
    }

    if (Dfield_ != 0.0) H += add_Dfield_(norb, Nl, Dfield_);
    if (enforce_hermitian) H = 0.5 * (H + H.adjoint());
    return H;
}

inline void RG_KPModel::write_info(std::ostream& os, const std::string& prefix) const
{
    os << prefix << "[Model] RG_KPModel\n";
    os << prefix << "Dfield_eV = " << Dfield_ << "\n";

    os << prefix << "kp_params (SWMcC):\n";
    os << prefix << "  gamma0   = " << p_.gamma0 << "\n";
    os << prefix << "  gamma1   = " << p_.gamma1 << "\n";
    os << prefix << "  gamma3   = " << p_.gamma3 << "\n";
    os << prefix << "  gamma4   = " << p_.gamma4 << "\n";
    os << prefix << "  gamma2   = " << p_.gamma2 << "\n";
    os << prefix << "  gamma5   = " << p_.gamma5 << "\n";
    os << prefix << "  gamma6   = " << p_.gamma6 << "\n";
    os << prefix << "  epsilon  = " << p_.epsilon << "\n";
    os << prefix << "  inc_g2   = " << (p_.include_gamma2 ? 1 : 0) << "\n";
    os << prefix << "  inc_g3   = " << (p_.include_gamma3 ? 1 : 0) << "\n";
    os << prefix << "  inc_g4   = " << (p_.include_gamma4 ? 1 : 0) << "\n";
    os << prefix << "  t_min    = " << p_.t_min << "\n";
}















} // namespace rg
