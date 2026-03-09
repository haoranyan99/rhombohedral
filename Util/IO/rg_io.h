// File: Util/IO/rg_io.h
#pragma once

#include <Eigen/Dense>

#include <filesystem>
#include <algorithm>
#include <cctype>
#include <limits>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <complex> 
#include <chrono>
#include <ctime>


#include "Common/DataContainers.h"
#include "PhysStruct/RG_Structure.h"
#include "Util/Config/config_reader.h"
#include "LinearAlgebra/MathFunctions.h"
#include "Source/Models/RG_ModelBase.h"

namespace rgio {

#pragma region ===================== Small TEXT utilities =====================

inline void cout_vec2(const std::string& name, const rg::Vec2& v) {
    std::cout << std::setw(10) << name << " = ("
              << std::setw(14) << v.x() << ", "
              << std::setw(14) << v.y() << ")\n";
}

inline void cout_vec3(const std::string& name, const rg::Vec3& v) {
    std::cout << std::setw(10) << name << " = ("
              << std::setw(14) << v.x() << ", "
              << std::setw(14) << v.y() << ", "
              << std::setw(14) << v.z() << ")\n";
}

inline std::string make_time_stamp(const std::string& prefix)
{
    using namespace std::chrono;

    const auto now = system_clock::now();
    const std::time_t tt = system_clock::to_time_t(now);

    std::tm tm{};
#ifdef _WIN32
    localtime_s(&tm, &tt);
#else
    localtime_r(&tt, &tm);
#endif

    std::ostringstream oss;
    oss << prefix << "_"
        << (tm.tm_year + 1900) << "_"
        << (tm.tm_mon + 1) << "_"
        << tm.tm_mday << "_"
        << tm.tm_hour << "_"
        << tm.tm_min << "_"
        << tm.tm_sec;

    return oss.str();
}

inline std::string trim(const std::string& s) {
    size_t a = 0, b = s.size();
    while (a < b && std::isspace((unsigned char)s[a])) ++a;
    while (b > a && std::isspace((unsigned char)s[b - 1])) --b;
    return s.substr(a, b - a);
}

inline bool is_hash_line(const std::string& s) {
    std::string t = trim(s);
    return (!t.empty() && t[0] == '#');
}

inline std::string tag_double_fixed(double x, int prec) {
    std::ostringstream oss;
    oss.setf(std::ios::fixed);
    oss << std::setprecision(prec) << x;
    return oss.str();
}
inline std::string tag3(double m)      { return tag_double_fixed(m, 3); }
inline std::string tag4(double m)      { return tag_double_fixed(m, 4); }
inline std::string tag5(double m)      { return tag_double_fixed(m, 5); }
inline std::string tag6(double m)      { return tag_double_fixed(m, 6); }
inline std::string tag7(double m)      { return tag_double_fixed(m, 7); }
inline std::string tag8(double m)      { return tag_double_fixed(m, 8); }
inline std::string tag9(double m)      { return tag_double_fixed(m, 9); }

inline std::string make_key(double D, double T, double x) {
    // folder tags are consistent => to_string is enough for keying
    return std::to_string(D) + "|" + std::to_string(T) + "|" + std::to_string(x);
}

#pragma endregion






#pragma region ===================== parse utilities =====================


inline bool parse_header_int(const std::string& nohash,
                                    const std::string& key,
                                    int& out)
{
    const auto pos = nohash.find(key);
    if (pos == std::string::npos) return false;
    const auto eq = nohash.find('=', pos);
    if (eq == std::string::npos) return false;

    std::string rhs = trim(nohash.substr(eq + 1));
    std::istringstream iss(rhs);
    int v;
    if (!(iss >> v)) return false;
    out = v;
    return true;
}

inline bool parse_header_double(const std::string& nohash,
                                     const std::string& key,
                                     double& out)
{
    const auto pos = nohash.find(key);
    if (pos == std::string::npos) return false;
    const auto eq = nohash.find('=', pos);
    if (eq == std::string::npos) return false;

    std::string rhs = trim(nohash.substr(eq + 1));
    std::istringstream iss(rhs);
    double v;
    if (!(iss >> v)) return false;
    out = v;
    return true;
}

// parse "# b1 = (x, y)" tolerant
inline bool parse_header_vec2(const std::string& nohash,
                                   const std::string& key,
                                   rg::Vec2& out)
{
    const auto pos = nohash.find(key);
    if (pos == std::string::npos) return false;
    const auto eq = nohash.find('=', pos);
    if (eq == std::string::npos) return false;

    std::string rhs = trim(nohash.substr(eq + 1));
    auto lp = rhs.find('(');
    auto rp = rhs.find(')');
    if (lp == std::string::npos || rp == std::string::npos || rp <= lp) return false;

    std::string inside = rhs.substr(lp + 1, rp - lp - 1);
    for (char& c : inside) if (c == ',') c = ' ';

    std::istringstream iss(inside);
    double x, y;
    if (!(iss >> x >> y)) return false;
    if (!std::isfinite(x) || !std::isfinite(y)) return false;

    out = rg::Vec2(x, y);
    return true;
}

// parse header like: "kmesh = BZ Nk=200 dk_frac=0.0003"
inline bool parse_kmesh_line(const std::string& nohash,
                                  std::string& kmesh,
                                  int& Nk,
                                  double& dk_frac)
{
    if (nohash.find("kmesh") == std::string::npos) return false;

    std::istringstream iss(nohash);
    std::string tok;
    std::vector<std::string> toks;
    while (iss >> tok) toks.push_back(tok);

    bool got = false;
    for (size_t i = 0; i < toks.size(); ++i) {
        if (toks[i] == "kmesh") {
            if (i + 2 < toks.size() && toks[i + 1] == "=") {
                kmesh = toks[i + 2];
                got = true;
            }
        } else if (toks[i].rfind("Nk=", 0) == 0) {
            Nk = std::stoi(toks[i].substr(3));
            got = true;
        } else if (toks[i].rfind("dk_frac=", 0) == 0) {
            dk_frac = std::stod(toks[i].substr(8));
            got = true;
        }
    }
    return got;
}

// parse folder tag like "D0.067", "T0.000", "mu0.945000", "doping0.0500"
inline bool parse_tag_value(const std::string& tag, const char* prefix, double& out) {
    const std::string p(prefix);
    if (tag.rfind(p, 0) != 0) return false; // not starts with prefix
    try {
        out = std::stod(tag.substr(p.size()));
        return true;
    } catch (...) {
        return false;
    }
}



#pragma endregion






#pragma region ===================== write output header =====================


inline std::string write_cal_fermi_header(
    const config::CalFermiConfig& cfg,
    const rg::RG_Structure& st,
    const rg::RG_ModelBase& model,
    double T_K,
    double EF,
    double doping,
    double filling
) {
    std::ostringstream hdr;

    hdr << std::setprecision(8);

    // ---- basic info ----
    hdr << "# --- calculate_fermi ---\n";
    hdr << "# model = " << cfg.model << "\n";
    hdr << "# mode = "  << cfg.mode  << "\n";

    // ---- structure & model ----
    st.write_info(hdr, "# ");
    model.write_info(hdr, "# ");

    // ---- physical parameters ----
    hdr << "# Dfield = " << cfg.Dfield_eV << "\n";
    hdr << "# T = " << T_K << "\n";
    hdr << "# EF = " << EF << "\n";
    hdr << "# doping = " << doping << "\n";
    hdr << "# filling = " << filling << "\n";

    // ---- kmesh ----
    hdr << "# kmesh = " << cfg.kmesh.type
        << " Nk=" << cfg.kmesh.Nk;

    if (cfg.kmesh.type != "BZ")
        hdr << " dk_frac=" << cfg.kmesh.dk_frac;

    hdr << "\n";

    return hdr.str();
}



inline std::string write_cal_chi_header(
    const config::CalChiConfig& cfg,
    const rg::RG_Structure& st,
    const rg::RG_ModelBase& model,
    double T_K,
    double EF,
    double doping,
    double filling,
    double polar_mu   // NEW
) {
    std::ostringstream hdr;
    hdr << std::setprecision(8);

    hdr << "# --- calculate_chi ---\n";
    hdr << "# model = "    << cfg.model << "\n";
    hdr << "# mode = "     << cfg.mode      << "\n";
    hdr << "# boundary = " << cfg.boundary << "\n";

    st.write_info(hdr, "# ");
    model.write_info(hdr, "# ");

    hdr << "# Dfield = "   << cfg.Dfield_eV << "\n";
    hdr << "# polar_mu = " << polar_mu      << "\n"; // CHANGED
    hdr << "# eta = "      << cfg.eta       << "\n";
    hdr << "# iq_range = [" << cfg.iq_min << ", " << cfg.iq_max << "]\n";
    hdr << "# jq_range = [" << cfg.jq_min << ", " << cfg.jq_max << "]\n";
    hdr << "# T = "      << T_K    << "\n";
    hdr << "# doping = " << doping << "\n";
    hdr << "# EF = "     << EF     << "\n";
    hdr << "# filling = " << filling << "\n";
    hdr << "# fermiPatch_path = " << cfg.fermi_patch_path << "\n";

    return hdr.str();
}


#pragma endregion




#pragma region ===================== write output datatxt =====================

inline void write_bands_txt(const std::string& out_path,
                            const core::PathData& path,
                            const Eigen::MatrixXd& E)
{
    const int Nk = static_cast<int>(path.k_list.size());
    if (Nk <= 0)
        throw std::runtime_error("rgio::write_bands_txt: empty path.k_list");

    if (E.rows() != Nk)
        throw std::runtime_error("rgio::write_bands_txt: E.rows != path.k_list.size()");

    const int nb = static_cast<int>(E.cols());
    if (nb <= 0)
        throw std::runtime_error("rgio::write_bands_txt: E.cols <= 0");

    // IMPORTANT: append (header is written by caller)
    std::ofstream fout(out_path, std::ios::app);
    if (!fout)
        throw std::runtime_error("rgio::write_bands_txt: cannot append " + out_path);

    // ---- columns (band-format info) ----
    if (!path.kline.empty())
        fout << "# columns: ik s kx ky";
    else
        fout << "# columns: ik kx ky";

    for (int b = 0; b < nb; ++b) fout << " E" << b;
    fout << "\n";

    // ---- data ----
    fout << std::setprecision(15);
    for (int ik = 0; ik < Nk; ++ik) {
        fout << ik << " ";

        if (!path.kline.empty())
            fout << path.kline[ik] << " ";

        fout << path.k_list[ik].x() << " "
             << path.k_list[ik].y();

        for (int b = 0; b < nb; ++b)
            fout << " " << E(ik, b);

        fout << "\n";
    }
}

inline void write_dos_txt(const std::string& out_path,
                          const core::SeriesData& dosS,
                          const std::string& header = "")
{
    dosS.assert_consistent();

    const auto& E   = dosS.x;
    const auto& dos = dosS.get<double>("dos").v;

    if (E.size() != dosS.size() || dos.size() != dosS.size())
        throw std::runtime_error("write_dos_txt: field size mismatch");

    // optional filling column
    const double* fill_ptr = nullptr;
    try {
        const auto& fill = dosS.get<double>("filling").v;
        if (fill.size() != dosS.size())
            throw std::runtime_error("write_dos_txt: filling field size mismatch");
        fill_ptr = fill.data();
    } catch (...) {
        fill_ptr = nullptr;
    }

    // optional doping column
    const double* dop_ptr = nullptr;
    try {
        const auto& dop = dosS.get<double>("doping").v;
        if (dop.size() != dosS.size())
            throw std::runtime_error("write_dos_txt: doping field size mismatch");
        dop_ptr = dop.data();
    } catch (...) {
        dop_ptr = nullptr;
    }

    std::ofstream fout(out_path);
    if (!fout)
        throw std::runtime_error("write_dos_txt: cannot write: " + out_path);

    fout << std::setprecision(15);

    // -------- header (multi-line, ensure '#') --------
    if (!header.empty()) {
        std::istringstream iss(header);
        std::string line;
        while (std::getline(iss, line)) {
            if (line.empty()) continue;
            if (line[0] == '#') fout << line << "\n";
            else                fout << "# " << line << "\n";
        }
    }

    // -------- fixed writer header --------
    fout << "# RG DOS data\n";

    // columns (auto-adapt)
    if (fill_ptr && dop_ptr) {
        fout << "# columns: i E(eV) filling doping_1e12cm^-2 DOS\n";
    } else if (fill_ptr && !dop_ptr) {
        fout << "# columns: i E(eV) filling DOS\n";
    } else if (!fill_ptr && dop_ptr) {
        fout << "# columns: i E(eV) doping_1e12cm^-2 DOS\n";
    } else {
        fout << "# columns: i E(eV) DOS\n";
    }

    fout << "# dE = " << dosS.dx << "\n";

    // -------- data --------
    for (size_t i = 0; i < dosS.size(); ++i) {
        fout << i << " " << E[i] << " ";
        if (fill_ptr) fout << fill_ptr[i] << " ";
        if (dop_ptr)  fout << dop_ptr[i]  << " ";
        fout << dos[i] << "\n";
    }
}

inline void write_chi_txt(const std::string& out_path,
                          const core::GridData& qgrid,
                          const std::string& header = "")
{
    using cd = std::complex<double>;

    qgrid.assert_consistent();

    const auto& chi   = qgrid.get<cd>("chi").v;
    const auto& qx    = qgrid.get<double>("qx").v;
    const auto& qy    = qgrid.get<double>("qy").v;
    const auto& nPair = qgrid.get<long long>("nKpair").v;

    const size_t N = qgrid.size();
    if (chi.size() != N || qx.size() != N || qy.size() != N || nPair.size() != N)
        throw std::runtime_error("write_chi_txt: field size mismatch");

    std::ofstream fout(out_path);
    if (!fout)
        throw std::runtime_error("write_chi_txt: cannot write: " + out_path);

    fout << std::setprecision(15);

    // -------- header (multi-line, ensure '#') --------
    if (!header.empty()) {
        std::istringstream iss(header);
        std::string line;
        while (std::getline(iss, line)) {
            if (line.empty()) continue;
            if (line[0] == '#') fout << line << "\n";
            else                fout << "# " << line << "\n";
        }
    }

    // -------- fixed writer header --------
    fout << "# RG chi(q) data\n";
    fout << "# qx,qy in 1/Angstrom; chi complex (arb.)\n";
    fout << "# dq_repr = " << qgrid.dx << "\n";

    // -------- data --------
    for (size_t i = 0; i < N; ++i) {
        fout << i << " "
             << qgrid.iq[i] << " " << qgrid.jq[i] << " "
             << qx[i] << " " << qy[i] << " "
             << chi[i].real() << " " << chi[i].imag() << " "
             << nPair[i];

        fout << "\n";
    }
}

inline void write_fermi_patch_txt(const std::string& out_path,
                                 const core::GridData& kgrid,
                                 const std::string& header = "",
                                 bool debug_bands = true)   // default: output E_b and f_b if present
{
    kgrid.assert_consistent();

    const auto& kvec    = kgrid.get<rg::Vec2>("kvec").v;
    const auto& occ_avg = kgrid.get<double>("occ_k_avg").v;

    const size_t N = kgrid.size();
    if (kvec.size() != N || occ_avg.size() != N)
        throw std::runtime_error("write_fermi_patch_txt: field size mismatch");

    // -------- inside field: optional --------
    // If missing, treat all points as inside.
    std::vector<unsigned char> inside_fallback;
    const std::vector<unsigned char>* inPtr = nullptr;

    try {
        inPtr = &kgrid.get<unsigned char>("inside").v;
        if (inPtr->size() != N)
            throw std::runtime_error("write_fermi_patch_txt: inside size mismatch");
    } catch (...) {
        inside_fallback.assign(N, (unsigned char)1);
        inPtr = &inside_fallback;
    }
    const auto& inList = *inPtr;

    // -------- optional band debug fields --------
    const std::vector<std::vector<double>>* evalsPtr = nullptr; // "evals": per-k vector<double>(dim)
    const std::vector<std::vector<double>>* occbPtr  = nullptr; // "occ_band": per-k vector<double>(dim)
    int dim = -1;

    if (debug_bands) {
        try {
            evalsPtr = &kgrid.get<std::vector<double>>("evals").v;
            occbPtr  = &kgrid.get<std::vector<double>>("occ_band").v;

            if (evalsPtr->size() != N || occbPtr->size() != N)
                throw std::runtime_error("write_fermi_patch_txt: evals/occ_band size mismatch");

            // infer dim from first inside point that has data
            for (size_t i = 0; i < N; ++i) {
                if (!inList[i]) continue;
                dim = (int)(*evalsPtr)[i].size();
                if (dim <= 0) break;
                if ((int)(*occbPtr)[i].size() != dim)
                    throw std::runtime_error("write_fermi_patch_txt: evals vs occ_band dim mismatch");
                break;
            }

            if (dim <= 0) { evalsPtr = nullptr; occbPtr = nullptr; }
        } catch (...) {
            evalsPtr = nullptr;
            occbPtr  = nullptr;
            dim = -1;
        }
    }

    std::ofstream fout(out_path);
    if (!fout)
        throw std::runtime_error("write_fermi_patch_txt: cannot write: " + out_path);

    // default numeric style (high precision for k, occ, f)
    fout << std::setprecision(15);

    // -------- header (multi-line) --------
    // header can be multiple lines; we ensure each line begins with '#'
    if (!header.empty()) {
        std::istringstream iss(header);
        std::string line;
        while (std::getline(iss, line)) {
            if (line.empty()) continue;
            if (line[0] == '#') fout << line << "\n";
            else                fout << "# " << line << "\n";
        }
    }

    // -------- fixed writer header --------
    fout << "# RG fermi patch (per-k averaged occupation)\n";
    fout << "# occ_k_avg(k) = (1 / N_band) * sum_b f(E_b(k), EF, T)\n";

    if (evalsPtr && occbPtr) {
        fout << "# debug: exporting per-band eigenvalues and occupations\n";
        fout << "# eigenvalues E_b are printed with fixed 4 decimals\n";
        fout << "# columns:\n";
        fout << "# idx  iq  jq  kx  ky  occ_k_avg";
        for (int b = 0; b < dim; ++b) fout << "  E" << b;
        for (int b = 0; b < dim; ++b) fout << "  f" << b;
        fout << "\n";
    } else {
        fout << "# columns:\n";
        fout << "# idx  iq  jq  kx  ky  occ_k_avg\n";
    }

    // -------- data --------
    size_t idx = 0;
    for (size_t i = 0; i < N; ++i) {
        if (!inList[i]) continue;   // only output inside
        ++idx;                      // 1-based index in output

        // keep high precision for k and occupations
        fout << idx << " "
             << kgrid.iq[i] << " "
             << kgrid.jq[i] << " "
             << kvec[i].x() << " "
             << kvec[i].y() << " "
             << occ_avg[i];

        if (evalsPtr && occbPtr) {
            const auto& E = (*evalsPtr)[i];
            const auto& f = (*occbPtr)[i];
            if ((int)E.size() != dim || (int)f.size() != dim)
                throw std::runtime_error("write_fermi_patch_txt: per-point band vector size mismatch");

            // eigenvalues: fixed 4 decimals
            fout << std::fixed << std::setprecision(4);
            for (int b = 0; b < dim; ++b) fout << " " << E[b];

            // restore to non-fixed + high precision for f_b
            fout.unsetf(std::ios::floatfield);
            fout << std::setprecision(15);
            for (int b = 0; b < dim; ++b) fout << " " << f[b];
        }

        fout << "\n";
    }

    if (idx == 0) {
        std::cerr << "[write_fermi_patch_txt] WARNING: no points written\n";
    } else {
        std::cout << "[write_fermi_patch_txt] wrote " << idx
                  << " points to " << out_path << "\n";
    }
}

#pragma endregion










} // namespace rgio
