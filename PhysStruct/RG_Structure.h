// File: PhysStruct/RG_Structure.h
#pragma once

#include <Eigen/Dense>
#include <vector>
#include <string>
#include <stdexcept>
#include <cmath>

#include "LinearAlgebra/Constants.h"
#include "Common/DataContainers.h"


namespace rg {

// -------------------- Types --------------------
using Vec2 = Eigen::Vector2d;
using Vec3 = Eigen::Vector3d;
using Mat3 = Eigen::Matrix3d;
using PathData = core::PathData;


enum class Sublattice : int { A = 0, B = 1 };

struct AtomSite {
    int layer;          // 0 ... nLayer-1
    Sublattice sub;     // A/B within that layer
    Vec3 r;             // Cartesian coordinate
};

// -------------------- RG_Structure --------------------
class RG_Structure {
public:
    RG_Structure() = default;

    RG_Structure(int nLayer,
                 double pressure,
                 double a0,
                 double d0,
                 double vacuum = 10.0)
    {
        build(nLayer, pressure, a0, d0, vacuum);
    }

    void build(int nLayer,
               double pressure,
               double a0,
               double d0,
               double vacuum = 10.0)
    {
        if (nLayer <= 0) throw std::runtime_error("RG_Structure: nLayer must be > 0");
        if (!(a0 > 0.0)) throw std::runtime_error("RG_Structure: a0 must be > 0");
        if (!(d0 > 0.0)) throw std::runtime_error("RG_Structure: d0 must be > 0");
        if (!(vacuum >= 0.0)) throw std::runtime_error("RG_Structure: vacuum must be >= 0");

        nLayer_   = nLayer;
        pressure_ = pressure;
        a0_       = a0;
        d0_       = d0;
        vacuum_   = vacuum;

        // --- 1) pressure -> a,d (match your python linear model) ---
        // a(P) = a0 - 0.00197*P
        // d(P) = d0 - 0.05710*P
        a_ = a0_ - 0.00197 * pressure_;
        d_ = d0_ - 0.05710 * pressure_;
        if (!(a_ > 0.0)) throw std::runtime_error("RG_Structure: a(P) <= 0, pressure too large?");
        if (!(d_ > 0.0)) throw std::runtime_error("RG_Structure: d(P) <= 0, pressure too large?");

        // --- 2) graphene primitive vectors (your convention) ---
        // a1 = (0, -a)
        // a2 = (sqrt(3)/2*a, 1/2*a)
        const double sqrt3 = std::sqrt(3.0);
        a1_ = Vec2(0.0, -a_);
        a2_ = Vec2(0.5 * sqrt3 * a_, 0.5 * a_);

        // --- 2.5) reciprocal vectors + high-sym points ---
        init_reciprocal_and_hsymm_();

        // --- 3) fractional points for ABC stacking ---
        // A=(0,0), B=(1/3,2/3), C=(2/3,1/3)
        fracA_ = Vec2(0.0, 0.0);
        fracB_ = Vec2(1.0/3.0, 2.0/3.0);
        fracC_ = Vec2(2.0/3.0, 1.0/3.0);

        // --- 4) build atoms (2 per layer) using layer%3 rule ---
        atoms_.clear();
        atoms_.reserve(static_cast<size_t>(2 * nLayer_));

        for (int l = 0; l < nLayer_; ++l) {
            Vec2 f0, f1;
            fractionalPairForLayer_(l, f0, f1); // (A,B) or (B,C) or (C,A)

            const double z = static_cast<double>(l) * d_;

            AtomSite s0;
            s0.layer = l;
            s0.sub   = Sublattice::A;
            s0.r     = Vec3(fracToXY_(f0).x(), fracToXY_(f0).y(), z);
            atoms_.push_back(s0);

            AtomSite s1;
            s1.layer = l;
            s1.sub   = Sublattice::B;
            s1.r     = Vec3(fracToXY_(f1).x(), fracToXY_(f1).y(), z);
            atoms_.push_back(s1);
        }

        // --- 5) 3x3 lattice matrix with vacuum along z ---
        lattice_.setZero();
        lattice_.row(0) = Vec3(a1_.x(), a1_.y(), 0.0);
        lattice_.row(1) = Vec3(a2_.x(), a2_.y(), 0.0);
        lattice_.row(2) = Vec3(0.0, 0.0, static_cast<double>(nLayer_) * d_ + vacuum_);
    }

    // -------------------- Accessors --------------------
    int nLayer() const { return nLayer_; }
    double pressure() const { return pressure_; }

    double a0() const { return a0_; }
    double d0() const { return d0_; }

    double a() const { return a_; }   // pressure-corrected
    double d() const { return d_; }   // pressure-corrected
    double vacuum() const { return vacuum_; }

    const Vec2& a1() const { return a1_; }
    const Vec2& a2() const { return a2_; }
    const Vec2& b1() const { return b1_; }
    const Vec2& b2() const { return b2_; }

    const Mat3& lattice() const { return lattice_; }
    const std::vector<AtomSite>& atoms() const { return atoms_; }

    // High-symmetry points (your convention)
    const Vec2& Gamma() const { return Gamma_; }
    const Vec2& M()     const { return M_; }
    const Vec2& K()     const { return K_; }
    const Vec2& Kp()    const { return Kp_; }

    void write_info(std::ostream& os, const std::string& prefix = "# ") const;

    // Convenient: get Cartesian position by (layer, sub)
    Vec3 sitePosition(int layer, Sublattice sub) const {
        if (layer < 0 || layer >= nLayer_)
            throw std::out_of_range("RG_Structure::sitePosition layer out of range");
        const int idx = 2 * layer + (sub == Sublattice::A ? 0 : 1);
        return atoms_.at(static_cast<size_t>(idx)).r;
    }

    core::GridData generate_localK_kmesh_square(int n_k, double dk_frac) const
    {
        using Vec2 = Eigen::Vector2d;

        if (n_k < 1)
            throw std::runtime_error("RG_Structure::generate_localK_kmesh_square: n_k must be >= 1");
        if (!(dk_frac > 0.0))
            throw std::runtime_error("RG_Structure::generate_localK_kmesh_square: dk_frac must be > 0");

        core::GridData grid;

        const int n_side = 2 * n_k + 1;
        const size_t N = size_t(n_side) * size_t(n_side);
        grid.resize(N);

        // store step as dimensionless fraction
        grid.dx = dk_frac;
        grid.dy = dk_frac;

        auto& kvec = grid.add<Vec2>("kvec").v;

        // optional convenience
        auto& kxF = grid.add<double>("kx").v;
        auto& kyF = grid.add<double>("ky").v;

        const Vec2 center = K_;
        const Vec2 g1 = dk_frac * b1_;
        const Vec2 g2 = dk_frac * b2_;

        size_t idx = 0;
        for (int ix = 0; ix < n_side; ++ix) {
            const int iq = ix - n_k;       // -n_k .. +n_k

            for (int iy = 0; iy < n_side; ++iy, ++idx) {
                const int jq = iy - n_k;   // -n_k .. +n_k

                grid.iq[idx] = iq;
                grid.jq[idx] = jq;

                const Vec2 k = center + double(iq) * g1 + double(jq) * g2;
                kvec[idx] = k;

                kxF[idx] = k.x();
                kyF[idx] = k.y();
            }
        }

        return grid;
    }

    core::GridData generate_localK_kmesh_b1b2(int n_k, double dk_frac) const
    {
        using Vec2 = Eigen::Vector2d;

        if (n_k < 1)
            throw std::runtime_error("RG_Structure::generate_localK_kmesh_b1b2: n_k must be >= 1");
        if (!(dk_frac > 0.0))
            throw std::runtime_error("RG_Structure::generate_localK_kmesh_b1b2: dk_frac must be > 0");

        core::GridData grid;

        const int n_side = 2 * n_k + 1;                 // points per axis
        const size_t N = size_t(n_side) * size_t(n_side);
        grid.resize(N);

        // step size: one integer step corresponds to dk_frac * b1/b2
        grid.dx = dk_frac;
        grid.dy = dk_frac;

        auto& kvec = grid.add<Vec2>("kvec").v;
        const Vec2 center = K_;
        const Vec2 g1 = dk_frac * b1_;
        const Vec2 g2 = dk_frac * b2_;

        size_t idx = 0;
        for (int ix = 0; ix < n_side; ++ix) {
            const int iq = ix - n_k;   // -n_k .. +n_k

            for (int iy = 0; iy < n_side; ++iy, ++idx) {
                const int jq = iy - n_k; // -n_k .. +n_k

                grid.iq[idx] = iq;
                grid.jq[idx] = jq;

                const Vec2 k = center
                            + double(iq) * g1
                            + double(jq) * g2;

                kvec[idx] = k;
            }
        }

        return grid;
    }

    core::GridData generate_localK_kmesh_hex_b1b2(int n_k, double dk_frac) const
    {
        using Vec2 = Eigen::Vector2d;

        if (n_k < 1)
            throw std::runtime_error("RG_Structure::generate_localK_kmesh_hex_b1b2: n_k must be >= 1");
        if (!(dk_frac > 0.0))
            throw std::runtime_error("RG_Structure::generate_localK_kmesh_hex_b1b2: dk_frac must be > 0");

        core::GridData grid;

        const int n_side = 2 * n_k + 1;                 // points per axis (square container)
        const size_t N = size_t(n_side) * size_t(n_side);
        grid.resize(N);

        // one integer step in iq/jq corresponds to dk_frac*b1 or dk_frac*b2
        grid.dx = dk_frac;
        grid.dy = dk_frac;

        auto& kvec  = grid.add<Vec2>("kvec").v;
        auto& inHex = grid.add<unsigned char>("inside").v;

        const Vec2 center = K_;
        const Vec2 b1 = b1_;
        const Vec2 b2 = b2_;

        const Vec2 g1 = dk_frac * b1;
        const Vec2 g2 = dk_frac * b2;

        // axial-hex metric: max(|i|,|j|,|i+j|) <= n_k
        auto inside_hex = [n_k](int iq, int jq) -> bool {
            const int a = std::abs(iq);
            const int b = std::abs(jq);
            const int c = std::abs(iq + jq);
            const int m = std::max(a, std::max(b, c));
            return (m <= n_k);
        };

        size_t idx = 0;
        for (int ix = 0; ix < n_side; ++ix) {
            const int iq = ix - n_k;     // -n_k .. +n_k

            for (int iy = 0; iy < n_side; ++iy, ++idx) {
                const int jq = iy - n_k; // -n_k .. +n_k

                grid.iq[idx] = iq;
                grid.jq[idx] = jq;

                const Vec2 k = center + double(iq) * g1 + double(jq) * g2;
                kvec[idx] = k;

                inHex[idx] = inside_hex(iq, jq) ? (unsigned char)1 : (unsigned char)0;
            }
        }

        return grid;
    }

    core::GridData generate_BZ_kmesh(int Nk1, int Nk2) const
    {
        using Vec2 = Eigen::Vector2d;

        if (Nk1 <= 0 || Nk2 <= 0)
            throw std::runtime_error("generate_BZ_kmesh: Nk1,Nk2 must be > 0");

        const int N1 = 2 * Nk1 + 1;
        const int N2 = 2 * Nk2 + 1;

        core::GridData grid;
        const size_t N = size_t(N1) * size_t(N2);
        grid.resize(N);
        grid.dx = 1.0 / double(N1);
        grid.dy = 1.0 / double(N2);

        auto& kvec = grid.add<Vec2>("kvec").v;
        auto& kxF  = grid.add<double>("kx").v;
        auto& kyF  = grid.add<double>("ky").v;

        const Vec2 b1 = b1_;
        const Vec2 b2 = b2_;

        size_t idx = 0;
        for (int i = 0; i < N1; ++i) {
            for (int j = 0; j < N2; ++j, ++idx) {

                const int di = i - Nk1;   // in [-Nk1, Nk1]
                const int dj = j - Nk2;   // in [-Nk2, Nk2]
                grid.iq[idx] = di;
                grid.jq[idx] = dj;

                // fractional coords: u,v in [-0.5, 0.5]
                const double u = double(di) / double(N1);
                const double v = double(dj) / double(N2);
                const Vec2 k = u * b1 + v * b2;

                kvec[idx] = k;
                kxF[idx]  = k.x();
                kyF[idx]  = k.y();
            }
        }

        return grid;
    }



    PathData generate_GMKG(int Nk) const
    {
        if (Nk < 2) throw std::runtime_error("RG_Structure::generate_GMKG: Nk >= 2 required");

        const Vec2 sec[4] = { Gamma_, M_, K_, Gamma_ };
        const int num_seg = 3;
        const int total   = num_seg * Nk + 1;

        PathData path;
        path.resize(total);

        int idx = 0;
        double s_acc = 0.0;

        auto append_segment = [&](const Vec2& k0, const Vec2& k1) {
            const Vec2 dk = k1 - k0;
            const double seglen = dk.norm();

            for (int ik = 0; ik < Nk; ++ik) {
                const double t = double(ik) / double(Nk); // [0,1)
                path.k_list[idx] = k0 + t * dk;
                path.kline[idx]  = s_acc + t * seglen;
                ++idx;
            }
            s_acc += seglen;
        };

        append_segment(sec[0], sec[1]); // Γ->M
        append_segment(sec[1], sec[2]); // M->K
        append_segment(sec[2], sec[3]); // K->Γ

        // final endpoint
        path.k_list[idx] = sec[3];
        path.kline[idx]  = s_acc;

        // xticks at segment boundaries
        path.xtick_pos = { 0.0, path.kline[Nk], path.kline[2*Nk], path.kline[3*Nk] };
        path.xtick_lab = { "Γ", "M", "K", "Γ" };

        return path;
    }

    PathData generate_localK_MKKp(int Nk, double frac) const
    {
        if (Nk < 2)  throw std::runtime_error("RG_Structure::generate_localK_MKKp: Nk >= 2");
        if (!(frac > 0.0)) throw std::runtime_error("RG_Structure::generate_localK_MKKp: frac > 0 required");

        const Vec2 k_start = K_ + frac * (M_  - K_);
        const Vec2 k_mid   = K_;
        const Vec2 k_end   = K_ + frac * (Kp_ - K_);

        const int num_seg = 2;
        const int total   = num_seg * Nk + 1;

        PathData path;
        path.resize(total);

        int idx = 0;
        double s_acc = 0.0;

        auto append_segment = [&](const Vec2& k0, const Vec2& k1) {
            const Vec2 dk = k1 - k0;
            const double seglen = dk.norm();

            for (int ik = 0; ik < Nk; ++ik) {
                const double t = double(ik) / double(Nk); // [0,1)
                path.k_list[idx] = k0 + t * dk;
                path.kline[idx]  = s_acc + t * seglen;
                ++idx;
            }
            s_acc += seglen;
        };

        append_segment(k_start, k_mid);
        append_segment(k_mid,   k_end);

        path.k_list[idx] = k_end;
        path.kline[idx]  = s_acc;

        path.xtick_pos = { 0.0, path.kline[Nk], path.kline[2*Nk] };
        path.xtick_lab = { "M", "K", "K'" };

        return path;
    }


private:
    Vec2 fracToXY_(const Vec2& f) const { return f.x() * a1_ + f.y() * a2_; }

    void fractionalPairForLayer_(int layer, Vec2& out0, Vec2& out1) const {
        const int m = layer % 3;
        if (m == 0)      { out0 = fracA_; out1 = fracB_; }
        else if (m == 1) { out0 = fracB_; out1 = fracC_; }
        else             { out0 = fracC_; out1 = fracA_; }
    }

    void init_reciprocal_and_hsymm_() {
        const double area = a1_.x() * a2_.y() - a1_.y() * a2_.x();
        if (std::abs(area) < 1e-15)
            throw std::runtime_error("RG_Structure: a1,a2 area too small (cannot build b1,b2)");

        const double twoPi = 2.0 * la::pi;

        b1_ = twoPi * Vec2( a2_.y(), -a2_.x()) / area;
        b2_ = twoPi * Vec2(-a1_.y(),  a1_.x()) / area;

        Gamma_ = Vec2(0.0, 0.0);
        M_     = 0.5 * b1_;
        K_     = (b1_ + b2_) / 3.0;
        Kp_    = 2.0 * (b1_ + b2_) / 3.0;
    }


private:
    int nLayer_ = 0;
    double pressure_ = 0.0;
    double a0_ = 0.0, d0_ = 0.0;
    double a_  = 0.0, d_  = 0.0;
    double vacuum_ = 10.0;

    Vec2 a1_{0,0}, a2_{0,0};
    Vec2 b1_{0,0}, b2_{0,0};

    Vec2 Gamma_{0,0}, M_{0,0}, K_{0,0}, Kp_{0,0};

    Vec2 fracA_{0,0}, fracB_{0,0}, fracC_{0,0};

    Mat3 lattice_ = Mat3::Zero();
    std::vector<AtomSite> atoms_;
};

inline void RG_Structure::write_info(std::ostream& os, const std::string& prefix) const
{
    os << prefix << "=== RG_Structure ===\n";
    os << prefix << "nLayer   = " << nLayer_ << "\n";
    os << prefix << "pressure = " << pressure_ << " (GPa)\n";
    os << prefix << "a0       = " << a0_ << " (A)\n";
    os << prefix << "d0       = " << d0_ << " (A)\n";
    os << prefix << "vacuum   = " << vacuum_ << " (A)\n";
    os << prefix << "a(P)     = " << a_ << " (A)\n";
    os << prefix << "d(P)     = " << d_ << " (A)\n";

    os << prefix << "a1 = (" << a1_.x() << ", " << a1_.y() << ")\n";
    os << prefix << "a2 = (" << a2_.x() << ", " << a2_.y() << ")\n";
    os << prefix << "b1 = (" << b1_.x() << ", " << b1_.y() << ")\n";
    os << prefix << "b2 = (" << b2_.x() << ", " << b2_.y() << ")\n";

    os << prefix << "Gamma = (" << Gamma_.x() << ", " << Gamma_.y() << ")\n";
    os << prefix << "M     = (" << M_.x() << ", " << M_.y() << ")\n";
    os << prefix << "K     = (" << K_.x() << ", " << K_.y() << ")\n";
    os << prefix << "K'    = (" << Kp_.x() << ", " << Kp_.y() << ")\n";
}


} // namespace rg
