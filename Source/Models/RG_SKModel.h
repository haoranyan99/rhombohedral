// File: Source/Models/RG_SKModel.h
#pragma once

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>
#include <sstream>
#include <limits>
#include <fstream>
#include <filesystem>
#include <nlohmann/json.hpp>
#include <iostream>

#include "Common/DataContainers.h"
#include "PhysStruct/RG_Structure.h"
#include "LinearAlgebra/Constants.h"
#include "LinearAlgebra/MathFunctions.h"
#include "Source/Models/RG_ModelBase.h"

namespace rg {

class RG_SKModel : public RG_ModelBase {
public:
    using RG_ModelBase::Vec2;
    using RG_ModelBase::Vec3;
    using RG_ModelBase::cd;
    using RG_ModelBase::get_structure;
    using RG_ModelBase::get_Dfield;
    using RG_ModelBase::set_Dfield;
    using RG_ModelBase::add_Dfield_;
    using RG_ModelBase::bands_at_k;
    using RG_ModelBase::bands_along_path;
    using RG_ModelBase::cal_dos_gaussian;
    using RG_ModelBase::cal_chi_grid_Ef;
    using RG_ModelBase::find_EF_from_filling;

    // Integer-format hopping record: m n x y z t
    struct HmnRInt {
        int m = 0;
        int n = 0;
        int x = 0;
        int y = 0;
        int z = 0;
        double t = 0.0;
    };

    // enlarged-atom mapping: enlarged index -> (home orbital n, translation x,y,z)
    struct NRInfo {
        int n = 0;
        int x = 0;
        int y = 0;
        int z = 0;
    };

    // SK parameters (pz-only) -- MATCH Python rgsk_tb.py
    struct RGSKParams {
        // V_{pi}, V_{sigma} in eV
        double Vpp_pi0    = -2.81; // eV
        double Vpp_sigma0 =  0.48; // eV

        // Wu-paper radial params (dimensionless q, length a in Angstrom)
        double q_pi    = 3.1451;
        double a_pi    = 1.418;
        double q_sigma = 7.428;
        double a_sigma = 3.349;

        // soft cutoff: 1/(1+exp((r-rc)/lc))
        double r_c = 6.14;  // Angstrom
        double l_c = 0.265; // Angstrom

        // neighbour enumeration cutoff (3D)
        double r_cut = 8.5; // Angstrom

        int    nR_max = 6;  // translations x,y in [-nR_max, nR_max]

        bool include_intralayer  = true;
        bool include_interlayer  = true;

        double r_min = 1e-6;
        double t_min = 1e-12;
    };

public:
    RG_SKModel() = default;

    explicit RG_SKModel(const RG_Structure& st);
    explicit RG_SKModel(const RG_Structure& st, const std::string& para_path);

    void set_para_file(const std::string& para_path) {
        if (para_path.empty()) {
            load_params_from_default_json_();
        } else {
            load_params_from_json_(std::filesystem::path(para_path));
        }
        hoppings_cached_.clear();
        hoppings_built_ = false;
    }

    RGSKParams&       get_params()       { return p_; }
    const RGSKParams& get_params() const { return p_; }

    const std::vector<HmnRInt>& ensure_hoppings_cached() const {
        if (!hoppings_built_ || hoppings_cached_.empty()) {
            build_hoppings_cached_();
        }
        return hoppings_cached_;
    };

    const std::vector<HmnRInt>& get_hoppings_cached() const { return hoppings_cached_; }

    // model-specific H(k)
    Eigen::MatrixXcd build_Hk(const Vec2& k, bool enforce_hermitian = true) const override;

    std::vector<HmnRInt> generate_HmnRInt_list_with_shear(
        const Vec2& first_layer_shear
    ) const;

    Eigen::MatrixXcd build_Hk_from_hoppings(
        const Vec2& k,
        const std::vector<HmnRInt>& hops,
        bool enforce_hermitian = true,
        bool include_Dfield = true
    ) const;

    void write_info(std::ostream& os, const std::string& prefix = "# ") const override;

private:
    static std::filesystem::path resolve_default_para_path_();

    void load_params_from_json_(const std::filesystem::path& path);

    void load_params_from_default_json_();

    std::vector<Vec3> atoms_cart_() const;

    void build_enlarged_atoms_(std::vector<Vec3>& enlarge_atoms,
                              std::vector<NRInfo>& nR_list) const;

    void build_enlarged_atoms_from_base_(const std::vector<Vec3>& base,
                                         std::vector<Vec3>& enlarge_atoms,
                                         std::vector<NRInfo>& nR_list) const;

    void build_pair_lists_cutoff_(const std::vector<Vec3>& atoms_cart,
                                  const std::vector<Vec3>& enlarge_atoms,
                                  const std::vector<NRInfo>& nR_list,
                                  std::vector<int>& ind1_list,
                                  std::vector<int>& ind2_list) const;

    void build_hoppings_cached_() const {
        if (st_.nLayer() <= 0 || st_.atoms().empty())
            throw std::runtime_error("RG_SKModel: structure not set or empty.");
        hoppings_cached_ = generate_HmnRInt_list_();
        hoppings_built_  = true;
    };                 // NOTE: const (mutable cache)

    std::vector<HmnRInt> cal_HmnR_list_sk_int_(const std::vector<Vec3>& atoms_cart,
                                              const std::vector<Vec3>& enlarge_atoms,
                                              const std::vector<int>& ind1_list,
                                              const std::vector<int>& ind2_list,
                                              const std::vector<NRInfo>& nR_list) const;

    std::vector<HmnRInt> generate_HmnRInt_list_() const;

    double soft_cut_(double r) const;          // 1/(1+exp((r-rc)/lc))
    double Vpp_pi_(double r) const;            // VPI * exp(qpi*(1-r/api)) * soft_cut
    double Vpp_sigma_(double r) const;         // VSIGMA * exp(qsigma*(1-r/asigma)) * soft_cut
    double sk_integral_pz_pz_(const Vec3& dr, double r) const;

private:
    RGSKParams p_;
    mutable std::vector<HmnRInt> hoppings_cached_;
    mutable bool hoppings_built_ = false;
};


// ============================
// inline/definitions
// ============================

inline std::vector<RG_SKModel::Vec3> RG_SKModel::atoms_cart_() const {
    const auto& a = st_.atoms();
    std::vector<Vec3> ret;
    ret.reserve(a.size());
    for (const auto& s : a) ret.push_back(s.r);
    return ret;
}

inline Eigen::MatrixXcd RG_SKModel::build_Hk(const Vec2& k, bool enforce_hermitian) const {
    return build_Hk_from_hoppings(
        k,
        ensure_hoppings_cached(),
        enforce_hermitian,
        true
    );
}

inline Eigen::MatrixXcd RG_SKModel::build_Hk_from_hoppings(
    const Vec2& k,
    const std::vector<HmnRInt>& hops,
    bool enforce_hermitian,
    bool include_Dfield
) const {
    const int norb = static_cast<int>(st_.atoms().size());
    Eigen::MatrixXcd H = Eigen::MatrixXcd::Zero(norb, norb);

    const Vec2 a1 = st_.a1();
    const Vec2 a2 = st_.a2();

    const std::complex<double> I(0.0, 1.0);

    // ---- hopping part ----
    for (const auto& h : hops) {
        const Vec2 R2 = double(h.x) * a1 + double(h.y) * a2;
        const double kdotR = k.dot(R2);
        const std::complex<double> phase = std::exp(I * kdotR);
        H(h.m, h.n) += h.t * phase;
    }

    // ---- static D-field part (use base's Dfield_) ----
    if (include_Dfield && Dfield_ != 0.0) {
        const int Nl = st_.nLayer();
        H += RG_ModelBase::add_Dfield_(norb, Nl, Dfield_);
    }

    if (enforce_hermitian) H = 0.5 * (H + H.adjoint());
    return H;
}

inline void RG_SKModel::build_enlarged_atoms_(std::vector<Vec3>& enlarge_atoms,
                                             std::vector<NRInfo>& nR_list) const
{
    if (st_.nLayer() <= 0) throw std::runtime_error("RG_SKModel: structure not set");
    const auto base = atoms_cart_();
    build_enlarged_atoms_from_base_(base, enlarge_atoms, nR_list);
}

inline void RG_SKModel::build_enlarged_atoms_from_base_(
    const std::vector<Vec3>& base,
    std::vector<Vec3>& enlarge_atoms,
    std::vector<NRInfo>& nR_list
) const
{
    const int norb = static_cast<int>(base.size());

    enlarge_atoms.clear();
    nR_list.clear();

    const Vec2 a1 = st_.a1();
    const Vec2 a2 = st_.a2();

    const int nR = p_.nR_max;
    enlarge_atoms.reserve(static_cast<size_t>((2*nR + 1) * (2*nR + 1) * norb));
    nR_list.reserve(enlarge_atoms.capacity());

    for (int x = -nR; x <= nR; ++x) {
        for (int y = -nR; y <= nR; ++y) {
            const Vec2 R2 = static_cast<double>(x) * a1 + static_cast<double>(y) * a2;
            const Vec3 R(R2.x(), R2.y(), 0.0);

            for (int n = 0; n < norb; ++n) {
                enlarge_atoms.push_back(base[n] + R);
                nR_list.push_back(NRInfo{n, x, y, 0});
            }
        }
    }
}

inline std::vector<RG_SKModel::HmnRInt>
RG_SKModel::generate_HmnRInt_list_with_shear(
    const Vec2& first_layer_shear
) const {
    std::vector<Vec3> base = atoms_cart_();
    const auto& meta = st_.atoms();

    if (base.size() != meta.size()) {
        throw std::runtime_error(
            "generate_HmnRInt_list_with_shear: atom metadata size mismatch"
        );
    }

    for (size_t i = 0; i < base.size(); ++i) {
        if (meta[i].layer == 0) {
            base[i].x() += first_layer_shear.x();
            base[i].y() += first_layer_shear.y();
        }
    }

    std::vector<Vec3> enlarge_atoms;
    std::vector<NRInfo> nR_list;
    build_enlarged_atoms_from_base_(base, enlarge_atoms, nR_list);

    std::vector<int> ind1_list, ind2_list;
    build_pair_lists_cutoff_(base, enlarge_atoms, nR_list, ind1_list, ind2_list);

    return cal_HmnR_list_sk_int_(base, enlarge_atoms, ind1_list, ind2_list, nR_list);
}

inline void RG_SKModel::build_pair_lists_cutoff_(const std::vector<Vec3>& atoms_cart,
                                                const std::vector<Vec3>& enlarge_atoms,
                                                const std::vector<NRInfo>& nR_list,
                                                std::vector<int>& ind1_list,
                                                std::vector<int>& ind2_list) const
{
    if (nR_list.size() != enlarge_atoms.size())
        throw std::runtime_error("build_pair_lists_cutoff_: nR_list size mismatch");

    const auto& meta = st_.atoms();
    if (meta.size() != atoms_cart.size())
        throw std::runtime_error("build_pair_lists_cutoff_: atoms size mismatch");

    ind1_list.clear();
    ind2_list.clear();

    const int norb = static_cast<int>(atoms_cart.size());
    const int nE   = static_cast<int>(enlarge_atoms.size());

    for (int m = 0; m < norb; ++m) {
        const Vec3 r1 = atoms_cart[m];

        for (int e = 0; e < nE; ++e) {
            const Vec3 dr = enlarge_atoms[e] - r1;
            const double r = dr.norm();

            if (r < p_.r_min) continue;
            if (r > p_.r_cut) continue;

            const int n = nR_list[e].n;
            const bool same_layer = (meta[m].layer == meta[n].layer);

            if (same_layer && !p_.include_intralayer) continue;
            if (!same_layer && !p_.include_interlayer) continue;

            ind1_list.push_back(m);
            ind2_list.push_back(e);
        }
    }
}

inline std::vector<RG_SKModel::HmnRInt>
RG_SKModel::cal_HmnR_list_sk_int_(const std::vector<Vec3>& atoms_cart,
                                 const std::vector<Vec3>& enlarge_atoms,
                                 const std::vector<int>& ind1_list,
                                 const std::vector<int>& ind2_list,
                                 const std::vector<NRInfo>& nR_list) const
{
    if (ind1_list.size() != ind2_list.size())
        throw std::runtime_error("cal_HmnR_list_sk_int_: pair list size mismatch");
    if (nR_list.size() != enlarge_atoms.size())
        throw std::runtime_error("cal_HmnR_list_sk_int_: nR_list size mismatch");

    const int num_pairs = static_cast<int>(ind1_list.size());
    std::vector<HmnRInt> ret;
    ret.reserve(static_cast<size_t>(num_pairs));

    for (int i = 0; i < num_pairs; ++i) {
        const int m = ind1_list[i];
        const int e = ind2_list[i];

        const Vec3 dr = enlarge_atoms[e] - atoms_cart[m];
        const double r = dr.norm();

        if (r < p_.r_min || r > p_.r_cut) continue;

        const double hopp = sk_integral_pz_pz_(dr, r);
        if (std::abs(hopp) < p_.t_min) continue;

        const auto info = nR_list[e];
        ret.push_back(HmnRInt{m, info.n, info.x, info.y, info.z, hopp});
    }
    return ret;
}

inline std::vector<RG_SKModel::HmnRInt> RG_SKModel::generate_HmnRInt_list_() const {
    const auto base = atoms_cart_();

    std::vector<Vec3> enlarge_atoms;
    std::vector<NRInfo> nR_list;
    build_enlarged_atoms_(enlarge_atoms, nR_list);

    std::vector<int> ind1_list, ind2_list;
    build_pair_lists_cutoff_(base, enlarge_atoms, nR_list, ind1_list, ind2_list);

    return cal_HmnR_list_sk_int_(base, enlarge_atoms, ind1_list, ind2_list, nR_list);
}

inline double RG_SKModel::soft_cut_(double r) const
{
    const double x = (r - p_.r_c) / p_.l_c;
    if (x >  50.0) return 0.0;
    if (x < -50.0) return 1.0;
    return 1.0 / (1.0 + std::exp(x));
}

inline double RG_SKModel::Vpp_pi_(double r) const
{
    const double expo = p_.q_pi * (1.0 - r / p_.a_pi);
    return p_.Vpp_pi0 * std::exp(expo) * soft_cut_(r);
}

inline double RG_SKModel::Vpp_sigma_(double r) const
{
    const double expo = p_.q_sigma * (1.0 - r / p_.a_sigma);
    return p_.Vpp_sigma0 * std::exp(expo) * soft_cut_(r);
}

inline double RG_SKModel::sk_integral_pz_pz_(const Vec3& dr, double r) const // r = |dr|
{
    const double z = dr.z();
    const double z2_over_r2 = (z*z) / (r*r);
    const double vpi = Vpp_pi_(r);
    const double vs  = Vpp_sigma_(r);
    return vpi * (1.0 - z2_over_r2) + vs * z2_over_r2;
}


// ============================
// constructors + json loader
// ============================

inline RG_SKModel::RG_SKModel(const RG_Structure& st)
    : RG_ModelBase(st)
{
    load_params_from_default_json_();
}

inline RG_SKModel::RG_SKModel(const RG_Structure& st, const std::string& para_path)
    : RG_ModelBase(st)
{
    if (para_path.empty()) load_params_from_default_json_();
    else                   load_params_from_json_(std::filesystem::path(para_path));
}


// Try a few common relative locations (minimal robust).
inline std::filesystem::path RG_SKModel::resolve_default_para_path_()
{
    const std::vector<std::filesystem::path> candidates = {
        "Source/Parameters/RG_para.json",
        "../Source/Parameters/RG_para.json",
        "../../Source/Parameters/RG_para.json",
        "../../../Source/Parameters/RG_para.json"
    };

    for (const auto& p : candidates) {
        std::error_code ec;
        if (std::filesystem::exists(p, ec) && !ec) return p;
    }
    return candidates.front(); // for error message
}

inline void RG_SKModel::load_params_from_default_json_()
{
    const auto path = resolve_default_para_path_();
    load_params_from_json_(path);
}

inline void RG_SKModel::load_params_from_json_(const std::filesystem::path& path)
{
    std::ifstream ifs(path);
    if (!ifs) {
        std::cerr
            << "[RG_SKModel] Warning: cannot open para json: "
            << path.string() << "\n"
            << "  -> fallback to default SK parameters\n";
        return;
    }

    nlohmann::json j;
    try {
        ifs >> j;
    } catch (...) {
        std::cerr
            << "[RG_SKModel] Warning: failed to parse json: "
            << path.string() << "\n"
            << "  -> fallback to default SK parameters\n";
        return;
    }

    if (!j.contains("sk_params") || !j["sk_params"].is_object()) {
        std::cerr
            << "[RG_SKModel] Warning: json missing object \"sk_params\": "
            << path.string() << "\n"
            << "  -> fallback to default SK parameters\n";
        return;
    }

    const auto& s = j["sk_params"];

    auto getd = [&](const char* key, double& ref) {
        if (s.contains(key) && s[key].is_number()) ref = s[key].get<double>();
    };
    auto geti = [&](const char* key, int& ref) {
        if (s.contains(key) && s[key].is_number_integer()) ref = s[key].get<int>();
    };
    auto getb = [&](const char* key, bool& ref) {
        if (s.contains(key) && s[key].is_boolean()) ref = s[key].get<bool>();
    };

    // ---- override defaults only when keys exist ----
    getd("Vpp_pi0",    p_.Vpp_pi0);
    getd("Vpp_sigma0", p_.Vpp_sigma0);

    getd("q_pi",    p_.q_pi);
    getd("a_pi",    p_.a_pi);
    getd("q_sigma", p_.q_sigma);
    getd("a_sigma", p_.a_sigma);

    getd("r_c", p_.r_c);
    getd("l_c", p_.l_c);

    getd("r_cut", p_.r_cut);
    geti("nR_max", p_.nR_max);

    getb("include_intralayer", p_.include_intralayer);
    getb("include_interlayer", p_.include_interlayer);

    getd("r_min", p_.r_min);
    getd("t_min", p_.t_min);
}

inline void RG_SKModel::write_info(std::ostream& os,
                                   const std::string& prefix) const
{
    os << prefix << "[Model] RG_SKModel\n";
    os << prefix << "Dfield_eV = " << Dfield_ << "\n";

    os << prefix << "sk_params:\n";
    os << prefix << "  Vpp_pi0        = " << p_.Vpp_pi0 << "\n";
    os << prefix << "  Vpp_sigma0     = " << p_.Vpp_sigma0 << "\n";
    os << prefix << "  q_pi           = " << p_.q_pi << "\n";
    os << prefix << "  a_pi           = " << p_.a_pi << "\n";
    os << prefix << "  q_sigma        = " << p_.q_sigma << "\n";
    os << prefix << "  a_sigma        = " << p_.a_sigma << "\n";
    os << prefix << "  r_c            = " << p_.r_c << "\n";
    os << prefix << "  l_c            = " << p_.l_c << "\n";
    os << prefix << "  r_cut          = " << p_.r_cut << "\n";
    os << prefix << "  nR_max         = " << p_.nR_max << "\n";
    os << prefix << "  include_intra  = " << (p_.include_intralayer ? 1 : 0) << "\n";
    os << prefix << "  include_inter  = " << (p_.include_interlayer ? 1 : 0) << "\n";
    os << prefix << "  r_min          = " << p_.r_min << "\n";
    os << prefix << "  t_min          = " << p_.t_min << "\n";
}

} // namespace rg
