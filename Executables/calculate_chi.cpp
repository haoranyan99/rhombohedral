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

inline std::string newest_fermiPatch_bin_in_dir(const fs::path& dir)
{
    if (!fs::exists(dir) || !fs::is_directory(dir)) {
        return {};
    }

    fs::file_time_type best_time{};
    fs::path best_path;
    bool found = false;

    for (const auto& ent : fs::directory_iterator(dir)) {
        if (!ent.is_regular_file()) continue;

        const fs::path p = ent.path();

        if (p.extension() != ".bin") continue;

        const std::string fname = p.filename().string();

        if (fname.rfind("fermiPatch_", 0) != 0) continue;

        const auto t = fs::last_write_time(p);

        if (!found || t > best_time) {
            best_time = t;
            best_path = p;
            found = true;
        }
    }

    return found ? best_path.string() : std::string{};
}

inline std::string make_mu_tag(double mu)
{
    if (std::abs(mu) < 5e-13) mu = 0.0;
    return std::string("mu") + rgio::tag6(mu);
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

        for (double polar_mu : cfg.polar_list) {

            const fs::path run_dir =
                fs::path(cfg.data_dir)
                / (
                    std::string("chi_")
                  + cfg.model
                  + "_mu_"
                  + cfg.output_suffix
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

                    const std::string fermi_bin =
                        newest_fermiPatch_bin_in_dir(mu_src);

                    if (fermi_bin.empty()) {
                        if (rank == 0) {
                            std::cout << "[base] MISS dir="
                                      << mu_src.string()
                                      << "\n";
                        }

                        rgmpi::barrier();
                        continue;
                    }

                    if (rank == 0) {
                        std::cout << "[base] FOUND file="
                                  << fermi_bin
                                  << "\n";
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

                    rgio::cal_chi_param param;

                    param.iq = cfg.iq;
                    param.jq = cfg.jq;

                    param.T_K = T_K;
                    param.Ef  = mu;

                    param.polar_mu = polar_mu;
                    param.eta = cfg.eta;

                    param.boundary_periodic = periodic;
                    param.use_form_factor = cfg.use_form_factor;

                    param.lattice_a = st.a();

                    param.area_density = 0.0;

                    rgio::cal_chi_result res =
                        rg::cal_chi_grid_from_fermiPatch(
                            fermi_bin,
                            param
                        );

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
                        ofs << "# iq = " << cfg.iq << "\n";
                        ofs << "# jq = " << cfg.jq << "\n";
                        ofs << "# q_unit = iq * dx * b1 + jq * dy * b2\n";
                        ofs << "# normalization = sum_k * area_density\n";
                        ofs << "# fermiPatch_bin = " << fermi_bin << "\n";

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
                                nK =
                                    nPairF[i]
                                    / (
                                        static_cast<long long>(param.dim)
                                      * static_cast<long long>(param.dim)
                                    );
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