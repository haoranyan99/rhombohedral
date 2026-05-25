#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef USE_MPI
  #include <mpi.h>
#endif

#include "Util/MPI/mpi.h"
#include "Util/Config/config_reader.h"
#include "Util/IO/rg_io.h"

#include "PhysStruct/RG_Structure.h"
#include "MeasureEngines/cal_susceptibility.h"

namespace fs = std::filesystem;

namespace {

inline std::string to_lower(std::string s)
{
    for (char& c : s) {
        c = static_cast<char>(
            std::tolower(static_cast<unsigned char>(c))
        );
    }
    return s;
}

inline std::vector<std::string> list_fermiPatch_bins_in_dir(const fs::path& dir)
{
    if (!fs::exists(dir) || !fs::is_directory(dir)) {
        return {};
    }

    std::vector<std::string> out;

    for (const auto& ent : fs::directory_iterator(dir)) {
        if (!ent.is_regular_file()) continue;

        const fs::path p = ent.path();

        if (p.extension() != ".bin") continue;

        const std::string fname = p.filename().string();

        const bool is_patch =
            fname.rfind("fermiPatch", 0) == 0
         || fname.rfind("fermi_spin", 0) == 0;

        if (!is_patch) continue;

        out.push_back(p.string());
    }

    std::sort(out.begin(), out.end());
    return out;
}

inline std::string make_mu_tag(double mu)
{
    if (std::abs(mu) < 5e-13) mu = 0.0;
    return std::string("mu") + rgio::tag6(mu);
}

inline std::string compact_double_tag(double x)
{
    if (std::abs(x) < 5e-13) {
        x = 0.0;
    }

    std::ostringstream oss;
    oss << std::fixed << std::setprecision(10) << x;

    std::string s = oss.str();

    while (!s.empty() && s.back() == '0') {
        s.pop_back();
    }

    if (!s.empty() && s.back() == '.') {
        s.pop_back();
    }

    if (s == "-0") {
        s = "0";
    }

    return s;
}

inline bool looks_like_timestamp_token(const std::string& token)
{
    int underscores = 0;

    if (token.empty()) {
        return false;
    }

    for (char c : token) {
        if (c == '_') {
            ++underscores;
            continue;
        }

        if (!std::isdigit(static_cast<unsigned char>(c))) {
            return false;
        }
    }

    return underscores >= 5;
}

inline std::string strip_trailing_timestamp(const std::string& s)
{
    std::string out = s;

    for (size_t pos = out.find('_');
         pos != std::string::npos;
         pos = out.find('_', pos + 1))
    {
        const std::string tail = out.substr(pos + 1);
        if (looks_like_timestamp_token(tail)) {
            out.erase(pos);
            break;
        }
    }

    return out;
}

inline std::string infer_fermi_suffix_from_path(
    const std::string& fermi_patch_path,
    const std::string& model
) {
    fs::path p(fermi_patch_path);

    if (p.filename().string().rfind("D", 0) == 0 && p.has_parent_path()) {
        p = p.parent_path();
    }

    std::string name = p.filename().string();

    const std::string prefix =
        std::string("fermi_") + model + "_mu_";

    if (name.rfind(prefix, 0) == 0) {
        name = name.substr(prefix.size());
    }

    name = strip_trailing_timestamp(name);

    return name.empty() ? "unknownFermi" : name;
}

inline std::string auto_chi_output_suffix(const config::CalChiConfig& cfg)
{
    std::ostringstream oss;

    oss << infer_fermi_suffix_from_path(
            cfg.fermi_patch_path,
            cfg.model
        )
        << "_q"
        << cfg.iq
        << "_"
        << cfg.jq
        << "_broadening"
        << compact_double_tag(1000.0 * cfg.eta)
        << "meV";

    if (cfg.diagonal_band_only) {
        oss << "_diagband";
    }

    return oss.str();
}

inline std::string chi_output_suffix(const config::CalChiConfig& cfg)
{
    std::string suffix =
        auto_chi_output_suffix(cfg);

    if (
        !cfg.output_suffix.empty()
     && cfg.output_suffix != "auto"
    ) {
        suffix += "_" + cfg.output_suffix;
    }

    return suffix;
}

} // namespace


int main(int argc, char** argv)
{
#ifdef USE_MPI
    MPI_Init(&argc, &argv);
#endif

    int rank = 0;
    int nprocs = 1;
    rgmpi::rank_size(rank, nprocs);

    try {
        const std::string config_file =
            (argc > 1) ? std::string(argv[1]) : "config.json";

        config::CalChiConfig cfg =
            config::read_cal_chi_config(config_file);

        const auto para =
            config::read_rg_para(cfg.para_file);

        cfg.model    = to_lower(cfg.model);
        cfg.boundary = to_lower(cfg.boundary);

        const bool periodic =
            (cfg.boundary == "periodic");

        rg::RG_Structure st(
            para.layer_num,
            para.pressure,
            para.a0,
            para.d0,
            para.vacuum
        );

        const std::string run_stamp =
            rgio::make_time_stamp("chi");

        const std::string output_suffix =
            chi_output_suffix(cfg);

        for (double polar_mu : cfg.polar_list) {

            const fs::path run_dir =
                fs::path(cfg.data_dir)
                / (
                    std::string("chi_")
                  + cfg.model
                  + "_mu_"
                  + output_suffix
                )
                / (std::string("D") + rgio::tag3(cfg.Dfield_eV))
                / (
                    std::string("polar_meV")
                  + rgio::tag3(1000.0 * polar_mu)
                );

            if (rank == 0) {
                fs::create_directories(run_dir);

                std::cout << "=== calculate_chi from fermiPatch bin ===\n";
                std::cout << "config          = " << config_file << "\n";
                std::cout << "model           = " << cfg.model << "\n";
                std::cout << "input_format    = bin\n";
                std::cout << "boundary        = " << cfg.boundary << "\n";
                std::cout << "Dfield_eV       = " << cfg.Dfield_eV << "\n";
                std::cout << "polar_meV       = " << 1000.0 * polar_mu << "\n";
                std::cout << "q integer       = (" << cfg.iq << ", "
                          << cfg.jq << ")\n";
                std::cout << "eta             = " << cfg.eta << "\n";
                std::cout << "use_form_factor = "
                          << (cfg.use_form_factor ? "true" : "false")
                          << "\n";
                std::cout << "diagonal_band_only = "
                          << (cfg.diagonal_band_only ? "true" : "false")
                          << "\n";
                std::cout << "suffix          = " << output_suffix << "\n";
                std::cout << "fermi root      = " << cfg.fermi_patch_path << "\n";
                std::cout << "out root        = " << run_dir.string() << "\n";
                std::cout << "MPI ranks       = " << nprocs << "\n";
            }

            rgmpi::barrier();

            for (double T_K : cfg.temperature_list) {

                const fs::path T_src =
                    fs::path(cfg.fermi_patch_path)
                    / (std::string("T") + rgio::tag3(T_K));

                for (double mu_raw : cfg.mu_list) {

                    double mu = mu_raw;
                    if (std::abs(mu) < 5e-13) {
                        mu = 0.0;
                    }

                    const fs::path mu_src =
                        T_src / make_mu_tag(mu);

                    if (rank == 0) {
                        std::cout << "------------------------------------------------------------\n";
                        std::cout << "[pt] T=" << rgio::tag3(T_K)
                                  << " mu=" << rgio::tag6(mu)
                                  << "\n";
                    }

                    const std::vector<std::string> fermi_bins =
                        list_fermiPatch_bins_in_dir(mu_src);

                    if (fermi_bins.empty()) {
                        if (rank == 0) {
                            std::cout << "[base] MISS dir="
                                      << mu_src.string()
                                      << "\n";
                        }

                        rgmpi::barrier();
                        continue;
                    }

                    if (rank == 0) {
                        std::cout << "[base] FOUND "
                                  << fermi_bins.size()
                                  << " fermiPatch file(s)\n";

                        for (const auto& fermi_bin : fermi_bins) {
                            std::cout << "       " << fermi_bin << "\n";
                        }
                    }

                    rgmpi::barrier();

                    const fs::path out_T =
                        run_dir / (std::string("T") + rgio::tag3(T_K));

                    const fs::path out_mu =
                        out_T / make_mu_tag(mu);

                    if (rank == 0) {
                        fs::create_directories(out_mu);
                    }

                    rgmpi::barrier();

                    rgio::cal_chi_param base_param;

                    base_param.iq = cfg.iq;
                    base_param.jq = cfg.jq;

                    base_param.T_K = T_K;
                    base_param.Ef  = mu;

                    base_param.polar_mu = polar_mu;
                    base_param.eta = cfg.eta;

                    base_param.boundary_periodic = periodic;
                    base_param.use_form_factor = cfg.use_form_factor;
                    base_param.diagonal_band_only = cfg.diagonal_band_only;

                    base_param.lattice_a = st.a();

                    base_param.area_density = 0.0;

                    rgio::cal_chi_param param = base_param;
                    rgio::cal_chi_result res;
                    bool have_res = false;
                    double doping_acc = 0.0;
                    double filling_acc = 0.0;
                    double occ_up_acc = 0.0;
                    double occ_dn_acc = 0.0;
                    std::vector<rgio::cal_chi_param> patch_params;

                    for (const auto& fermi_bin : fermi_bins) {
                        rgio::cal_chi_param patch_param = base_param;

                        rgio::cal_chi_result patch_res =
                            rg::cal_chi_grid_from_fermiPatch(
                                fermi_bin,
                                patch_param
                            );

                        if (!have_res) {
                            res = patch_res;
                            param = patch_param;
                            have_res = true;
                        } else {
                            core::GridData& qgrid_avg = res.qgrid;
                            const core::GridData& qgrid_patch =
                                patch_res.qgrid;

                            if (qgrid_avg.size() != qgrid_patch.size()) {
                                throw std::runtime_error(
                                    "calculate_chi: qgrid size mismatch while averaging patches"
                                );
                            }

                            auto& chi_avg =
                                qgrid_avg.get<std::complex<double>>("chi").v;

                            const auto& chi_patch =
                                qgrid_patch.get<std::complex<double>>("chi").v;

                            auto& nPair_avg =
                                qgrid_avg.get<long long>("nKpair").v;

                            const auto& nPair_patch =
                                qgrid_patch.get<long long>("nKpair").v;

                            for (size_t i = 0; i < qgrid_avg.size(); ++i) {
                                chi_avg[i] += chi_patch[i];
                                nPair_avg[i] += nPair_patch[i];
                            }
                        }

                        doping_acc += patch_res.doping;
                        filling_acc += patch_res.filling;
                        occ_up_acc += patch_res.occTot_up;
                        occ_dn_acc += patch_res.occTot_dn;
                        patch_params.push_back(patch_param);
                    }

                    if (!have_res) {
                        throw std::runtime_error(
                            "calculate_chi: no fermiPatch results to average"
                        );
                    }

                    const double inv_patch_count =
                        1.0 / static_cast<double>(fermi_bins.size());

                    res.doping = doping_acc * inv_patch_count;
                    res.filling = filling_acc * inv_patch_count;
                    res.occTot_up = occ_up_acc * inv_patch_count;
                    res.occTot_dn = occ_dn_acc * inv_patch_count;

                    auto& chi_avg =
                        res.qgrid.get<std::complex<double>>("chi").v;

                    auto& nPair_avg =
                        res.qgrid.get<long long>("nKpair").v;

                    for (size_t i = 0; i < res.qgrid.size(); ++i) {
                        chi_avg[i] *= inv_patch_count;
                        nPair_avg[i] =
                            static_cast<long long>(
                                std::llround(
                                    static_cast<double>(nPair_avg[i])
                                  * inv_patch_count
                                )
                            );
                    }

                    core::GridData& qgrid =
                        res.qgrid;

                    if (rank == 0) {
                        const fs::path out_path =
                            out_mu
                            / (std::string("chi_") + run_stamp + ".txt");

                        std::ofstream ofs(out_path.string());

                        if (!ofs) {
                            throw std::runtime_error(
                                "Cannot open output chi file: "
                                + out_path.string()
                            );
                        }

                        ofs << std::setprecision(15);

                        ofs << "# --- calculate_chi ---\n";
                        ofs << "# model = " << cfg.model << "\n";
                        ofs << "# input_format = bin\n";
                        ofs << "# boundary = " << cfg.boundary << "\n";
                        ofs << "# Dfield_eV = " << cfg.Dfield_eV << "\n";
                        ofs << "# polar_mu_eV = " << polar_mu << "\n";
                        ofs << "# polar_meV = " << 1000.0 * polar_mu << "\n";
                        ofs << "# T_K = " << param.T_K << "\n";
                        ofs << "# EF = " << param.Ef << "\n";
                        ofs << "# doping = " << res.doping << "\n";
                        ofs << "# filling = " << res.filling << "\n";
                        ofs << "# eta = " << cfg.eta << "\n";
                        ofs << "# use_form_factor = "
                            << (cfg.use_form_factor ? "true" : "false")
                            << "\n";
                        ofs << "# diagonal_band_only = "
                            << (cfg.diagonal_band_only ? "true" : "false")
                            << "\n";
                        ofs << "# band_scattering = "
                            << (cfg.diagonal_band_only ? "diagonal" : "all")
                            << "\n";
                        ofs << "# form_factor_kind = "
                            << (
                                cfg.use_form_factor
                                ? "density_overlap"
                                : "none"
                               )
                            << "\n";
                        ofs << "# iq = " << cfg.iq << "\n";
                        ofs << "# jq = " << cfg.jq << "\n";
                        ofs << "# q_unit = iq * dx * b1 + jq * dy * b2\n";
                        ofs << "# normalization = sum_k * area_density\n";
                        ofs << "# fermiPatch_count = "
                            << fermi_bins.size() << "\n";
                        ofs << "# fermiPatch_average = arithmetic_mean\n";

                        for (size_t ibin = 0; ibin < fermi_bins.size(); ++ibin) {
                            ofs << "# fermiPatch_bin_" << ibin
                                << " = " << fermi_bins[ibin] << "\n";

                            if (ibin < patch_params.size()
                             && patch_params[ibin].has_spin_metadata)
                            {
                                ofs << "# fermiPatch_spin_sign_" << ibin
                                    << " = " << patch_params[ibin].spin_sign
                                    << "\n";
                                ofs << "# fermiPatch_doping_spin_" << ibin
                                    << " = " << patch_params[ibin].doping_spin
                                    << "\n";
                                ofs << "# fermiPatch_filling_spin_" << ibin
                                    << " = " << patch_params[ibin].filling_spin
                                    << "\n";
                            }
                        }

                        if (qgrid.dx > 0.0 && qgrid.dy > 0.0) {
                            ofs << "# dx = " << qgrid.dx << "\n";
                            ofs << "# dy = " << qgrid.dy << "\n";
                            ofs << "# area_density = "
                                << qgrid.dx * qgrid.dy << "\n";
                        }

                        ofs << "\n";
                        ofs << "# iq jq qx qy chi_real chi_imag nKpair nK\n";

                        const auto& chiF =
                            qgrid.get<std::complex<double>>("chi").v;

                        const auto& qxF =
                            qgrid.get<double>("qx").v;

                        const auto& qyF =
                            qgrid.get<double>("qy").v;

                        const auto& nPairF =
                            qgrid.get<long long>("nKpair").v;

                        for (size_t i = 0; i < qgrid.size(); ++i) {
                            long long nK = 0;

                            if (param.dim > 0) {
                                const long long band_pair_factor =
                                    cfg.diagonal_band_only
                                    ? static_cast<long long>(param.dim)
                                    : (
                                        static_cast<long long>(param.dim)
                                      * static_cast<long long>(param.dim)
                                      );

                                nK =
                                    nPairF[i]
                                    / band_pair_factor;
                            }

                            ofs << qgrid.iq[i] << " "
                                << qgrid.jq[i] << " "
                                << qxF[i] << " "
                                << qyF[i] << " "
                                << chiF[i].real() << " "
                                << chiF[i].imag() << " "
                                << nPairF[i] << " "
                                << nK << "\n";
                        }

                        std::cout << "[out] Wrote: "
                                  << out_path.string()
                                  << "\n";
                    }

                    rgmpi::barrier();
                }
            }
        }

        if (rank == 0) {
            std::cout << "DONE.\n";
        }

        rgmpi::finalize();
        return 0;

    } catch (const std::exception& ex) {
        if (rank == 0) {
            std::cerr << "[ERROR] " << ex.what() << "\n";
        }

        try {
            rgmpi::abort_all(
                std::string("calculate_chi FAILED: ") + ex.what(),
                1
            );
        } catch (...) {}

        rgmpi::finalize();
        return 1;
    }
}
