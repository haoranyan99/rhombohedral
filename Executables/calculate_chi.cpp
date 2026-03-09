// File: Executables/calculate_chi.cpp
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <memory>
#include <filesystem>
#include <cmath>

#ifdef USE_MPI
  #include <mpi.h>
#endif

#include "Util/MPI/mpi.h"
#include "Util/Config/config_reader.h"
#include "Util/IO/rg_io.h"
#include "Util/IO/rg_fermiPatch.h"
#include "PhysStruct/RG_Structure.h"
#include "Source/Models/RG_ModelBase.h"
#include "Source/Models/RG_SKModel.h"
#include "Source/Models/RG_KPModel.h"
#include "MeasureEngines/cal_susceptibility.h"

namespace fs = std::filesystem;

namespace {

inline std::string to_lower(std::string s)
{
    for (char& c : s) c = (char)std::tolower((unsigned char)c);
    return s;
}

// pick newest fermiPatch_*.txt in a folder
inline std::string newest_fermiPatch_in_dir(const fs::path& dir)
{
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
}

inline void print_sep(int rank)
{
    if (rank == 0) std::cout << "------------------------------------------------------------\n";
}

} // namespace

int main(int argc, char** argv)
{
#ifdef USE_MPI
    MPI_Init(&argc, &argv);
#endif

    int rank = 0, nprocs = 1;
    rgmpi::rank_size(rank, nprocs);

    try {
        const std::string config_file = (argc > 1) ? std::string(argv[1]) : "config.json";

        config::CalChiConfig cfg = config::read_cal_chi_config(config_file);
        const auto para = config::read_rg_para(cfg.para_file);

        cfg.mode     = cfg.mode.empty()     ? std::string("doping") : to_lower(cfg.mode);
        cfg.boundary = cfg.boundary.empty() ? std::string("open")   : to_lower(cfg.boundary);

        if (cfg.mode != "doping" && cfg.mode != "mu")
            rgmpi::abort_all("calculate_chi: cfg.mode must be \"doping\" or \"mu\"");

        const bool periodic = (cfg.boundary == "periodic");

        if (cfg.temperature_list.empty()) rgmpi::abort_all("calculate_chi: empty temperature_list");
        if (cfg.polar_list.empty())       rgmpi::abort_all("calculate_chi: empty polar_list");
        if (cfg.fermi_patch_path.empty()) rgmpi::abort_all("calculate_chi: empty fermi_patch_path");

        if (cfg.mode == "doping") {
            if (cfg.doping_list.empty()) rgmpi::abort_all("calculate_chi: empty doping_list (mode=doping)");
        } else {
            if (cfg.mu_list.empty()) rgmpi::abort_all("calculate_chi: empty mu_list (mode=mu)");
        }

        // structure + model (header info only)
        const double pressure = 0.0;
        rg::RG_Structure st(para.layer_num, pressure, para.a0, para.d0, para.vacuum);

        std::unique_ptr<rg::RG_ModelBase> model_ptr;
        if (cfg.model == "sk")      model_ptr = std::make_unique<rg::RG_SKModel>(st, cfg.para_file);
        else if (cfg.model == "kp") model_ptr = std::make_unique<rg::RG_KPModel>(st);
        else rgmpi::abort_all("calculate_chi: model must be \"sk\" or \"kp\"");

        rg::RG_ModelBase& model = *model_ptr;
        if (cfg.model == "sk") {
            if (auto* sk = dynamic_cast<rg::RG_SKModel*>(model_ptr.get()))
                sk->ensure_hoppings_cached();
        }
        model.set_Dfield(cfg.Dfield_eV);

        const std::vector<double>& X_list = (cfg.mode == "mu") ? cfg.mu_list : cfg.doping_list;
        const std::string run_stamp = rgio::make_time_stamp("chi");

        for (double polar_mu : cfg.polar_list) {
            if (!std::isfinite(polar_mu)) continue;

            const fs::path run_dir =
                (cfg.data_dir.empty() ? fs::path("data") : fs::path(cfg.data_dir))
                / (std::string("chi_") + cfg.model + "_" + cfg.mode + "_" + cfg.output_suffix)
                / (std::string("D") + rgio::tag3(cfg.Dfield_eV))
                / (std::string("polar_meV") + rgio::tag3(1000.0 * polar_mu));

            if (rank == 0) {
                std::error_code ec;
                fs::create_directories(run_dir, ec);
                if (ec) throw std::runtime_error("Failed to create output directory: " + run_dir.string());
                std::cout << "[chi] polar_meV=" << rgio::tag3(1000.0 * polar_mu)
                          << "  out=" << run_dir.string() << "\n";
            }
            rgmpi::barrier();

            for (double T : cfg.temperature_list) {
                const fs::path T_src = fs::path(cfg.fermi_patch_path) / (std::string("T") + rgio::tag3(T));

                for (double x : X_list) {
                    print_sep(rank);

                    // ---- print point header ----
                    if (rank == 0) {
                        if (cfg.mode == "mu") {
                            std::cout << "[pt] T=" << rgio::tag3(T)
                                      << "  polar_meV=" << rgio::tag3(1000.0 * polar_mu)
                                      << "  mu=" << rgio::tag6(x) << "\n";
                        } else {
                            std::cout << "[pt] T=" << rgio::tag3(T)
                                      << "  polar_meV=" << rgio::tag3(1000.0 * polar_mu)
                                      << "  doping=" << rgio::tag4(x) << "\n";
                        }
                    }

                    // ---- locate base fermiPatch ----
                    const fs::path x_src = (cfg.mode == "mu")
                        ? (T_src / (std::string("mu")     + rgio::tag6(x)))
                        : (T_src / (std::string("doping") + rgio::tag4(x)));

                    const std::string base_file = newest_fermiPatch_in_dir(x_src);
                    if (base_file.empty()) {
                        if (rank == 0) std::cout << "[base] MISS  dir=" << x_src.string() << "\n";
                        rgmpi::barrier();
                        continue;
                    }
                    if (rank == 0) std::cout << "[base] FOUND file=" << base_file << "\n";
                    rgmpi::barrier();

                    // ---- output dirs ----
                    const fs::path out_T = run_dir / (std::string("T") + rgio::tag3(T));
                    const fs::path out_x = (cfg.mode == "doping")
                        ? (out_T / (std::string("doping") + rgio::tag4(x)))
                        : (out_T / (std::string("mu")     + rgio::tag6(x)));

                    if (rank == 0) {
                        std::error_code ec;
                        fs::create_directories(out_x, ec);
                        if (ec) throw std::runtime_error("Failed to create output directory: " + out_x.string());
                    }
                    rgmpi::barrier();

                    // ---- compute chi ----
                    rgio::cal_chi_param param;
                    param.q_range = { cfg.iq_min, cfg.iq_max, cfg.jq_min, cfg.jq_max };
                    param.polar_mu = polar_mu;
                    param.boundary_periodic = periodic;
                    param.eta = cfg.eta;
                    param.T_K = T;
                    if (cfg.mode == "doping") param.doping = x;
                    else                      param.Ef     = x;

                    rgio::cal_chi_result res = rg::cal_chi_grid_from_fermiPatch(base_file, param);
                    core::GridData& qgrid = res.qgrid;

                    if (rank == 0) {
                        const double EF      = param.Ef;
                        const double doping  = res.doping;
                        const double filling = res.filling;

                        const std::string hdr =
                            rgio::write_cal_chi_header(cfg, st, model, param.T_K, EF, doping, filling, polar_mu);

                        const fs::path out_path = out_x / (std::string("chi_") + run_stamp + ".txt");
                        rgio::write_chi_txt(out_path.string(), qgrid, hdr);
                        std::cout << "[out] Wrote: " << out_path.string() << "\n";
                    }

                    rgmpi::barrier();
                }
            }
        }

        if (rank == 0) std::cout << "DONE.\n";

        rgmpi::finalize();
        return 0;

    } catch (const std::exception& ex) {
        if (rank == 0) std::cerr << "[ERROR] " << ex.what() << "\n";
        try { rgmpi::abort_all(std::string("calculate_chi FAILED: ") + ex.what(), 1); }
        catch (...) {}
        rgmpi::finalize();
        return 1;
    }
}