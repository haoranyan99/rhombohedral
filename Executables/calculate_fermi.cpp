// File: calculate_fermi.cpp
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <memory>
#include <limits>
#include <filesystem>
#include <algorithm>

#ifdef USE_MPI
  #include <mpi.h>
#endif

#include "Util/MPI/mpi.h"
#include "PhysStruct/RG_Structure.h"
#include "Source/Models/RG_SKModel.h"
#include "Source/Models/RG_KPModel.h"
#include "Util/Config/config_reader.h"
#include "Util/IO/rg_io.h"
#include "LinearAlgebra/MathFunctions.h"   // la::filling_to_doping


static std::pair<double,double> minmax_or_nan(const std::vector<double>& v) {
    if (v.empty()) {
        const double NaN = std::numeric_limits<double>::quiet_NaN();
        return {NaN, NaN};
    }
    auto [mn_it, mx_it] = std::minmax_element(v.begin(), v.end());
    return {*mn_it, *mx_it};
}

int main(int argc, char** argv) {
    try {

#ifdef USE_MPI
        MPI_Init(&argc, &argv);
#endif

        int rank = 0, nprocs = 1;
        rgmpi::rank_size(rank, nprocs);

        const std::string config_file = (argc > 1) ? std::string(argv[1]) : "config.json";
        const std::string run_stamp   = rgio::make_time_stamp("");

        // -----------------------------
        // 1) read config + parameters
        // -----------------------------
        const auto cfg  = config::read_cal_fermi_config(config_file);
        const auto para = config::read_rg_para(cfg.para_file);

        // -----------------------------
        // 1.5) sanity (mode-dependent)
        // -----------------------------
        if (cfg.temperature_list.empty())
            rgmpi::abort_all("calculate_fermi: empty temperature_list");

        const std::string mode = cfg.mode.empty() ? std::string("doping") : cfg.mode;

        if (mode == "doping") {
            if (cfg.doping_list.empty())
                rgmpi::abort_all("calculate_fermi: empty doping_list (mode=doping)");
            if (cfg.filling_list.empty())
                rgmpi::abort_all("calculate_fermi: empty filling_list (mode=doping)");
            if (cfg.doping_list.size() != cfg.filling_list.size())
                rgmpi::abort_all("calculate_fermi: doping_list.size != filling_list.size (mode=doping)");
        } else if (mode == "mu") {
            if (cfg.mu_list.empty())
                rgmpi::abort_all("calculate_fermi: empty mu_list (mode=mu)");
        } else {
            rgmpi::abort_all("calculate_fermi: unsupported mode = " + mode +
                             " (must be \"doping\" or \"mu\")");
        }

        // -----------------------------
        // 2) build structure
        // -----------------------------
        const double pressure = 0.0;
        rg::RG_Structure st(para.layer_num, pressure, para.a0, para.d0, para.vacuum);

        // -----------------------------
        // 3) build kpatch
        // -----------------------------
        core::GridData kpatch;
        if (cfg.kmesh.type == "localK_hex") {
            kpatch = st.generate_localK_kmesh_hex_b1b2(cfg.kmesh.Nk, cfg.kmesh.dk_frac);
        } else if (cfg.kmesh.type == "BZ") {
            kpatch = st.generate_BZ_kmesh(cfg.kmesh.Nk, cfg.kmesh.Nk);
        } else {
            rgmpi::abort_all("calculate_fermi: unsupported kmesh.type = " + cfg.kmesh.type);
        }

        // -----------------------------
        // 4) build model
        // -----------------------------
        std::unique_ptr<rg::RG_ModelBase> model_ptr;
        if (cfg.model == "sk") {
            model_ptr = std::make_unique<rg::RG_SKModel>(st);
        } else if (cfg.model == "kp") {
            model_ptr = std::make_unique<rg::RG_KPModel>(st);
        } else {
            rgmpi::abort_all("calculate_fermi: model must be \"sk\" or \"kp\"");
        }
        rg::RG_ModelBase& model = *model_ptr;

        if (cfg.model == "sk") {
            if (auto* sk = dynamic_cast<rg::RG_SKModel*>(model_ptr.get())) {
                sk->ensure_hoppings_cached();
            }
        }
        model.set_Dfield(cfg.Dfield_eV);

        // -----------------------------
        // 5) output roots:
        //    {data_dir}/fermi_{model}_{mode}_{suffix}/D.../T...
        // -----------------------------
        namespace fs = std::filesystem;

        const fs::path base_dir = cfg.data_dir.empty() ? fs::path("data") : fs::path(cfg.data_dir);
        const fs::path root_dir = base_dir / (
            std::string("fermi_") + cfg.model + "_" + mode + "_" + cfg.output_suffix
        );
        const fs::path D_dir    = root_dir / (std::string("D") + rgio::tag3(cfg.Dfield_eV));

        const auto [Tmin,  Tmax ] = minmax_or_nan(cfg.temperature_list);
        const auto [dop_min, dop_max] = minmax_or_nan(cfg.doping_list);
        const auto [mu_min,  mu_max ] = minmax_or_nan(cfg.mu_list);

        if (rank == 0) {
            std::error_code ec;
            fs::create_directories(D_dir, ec);
            if (ec) {
                throw std::runtime_error("Failed to create output directory: " +
                                         D_dir.string() + " | " + ec.message());
            }

            std::cout << "=== calculate_fermi ===\n";
            std::cout << "config   = " << config_file << "\n";
            std::cout << "run_stamp= " << run_stamp << " (used only in filenames)\n";
            std::cout << "data_dir = " << base_dir.string() << "\n";
            std::cout << "out_root = " << D_dir.string() << "\n";
            std::cout << "model    = " << cfg.model << "\n";
            std::cout << "mode     = " << mode << "\n";
            std::cout << "suffix   = " << cfg.output_suffix << "\n";
            std::cout << "layers   = " << para.layer_num << "\n";
            std::cout << "pressure = " << pressure << "\n";
            std::cout << "Dfield   = " << cfg.Dfield_eV << "\n";
            std::cout << "kmesh    = " << cfg.kmesh.type
                      << " Nk=" << cfg.kmesh.Nk
                      << " dk_frac=" << cfg.kmesh.dk_frac << "\n";
#ifdef USE_MPI
            std::cout << "MPI      = " << nprocs
                      << " ranks (inited=" << (rgmpi::inited() ? "true" : "false") << ")\n";
#endif
            std::cout << std::fixed;
            std::cout << "scan T: count=" << cfg.temperature_list.size()
                      << " in [" << std::setprecision(6) << Tmin
                      << ", "   << std::setprecision(6) << Tmax << "]\n";

            if (mode == "doping") {
                std::cout << "scan x(doping): count=" << cfg.doping_list.size()
                          << " in [" << std::setprecision(12) << dop_min
                          << ", "   << std::setprecision(12) << dop_max << "]\n";
            } else {
                std::cout << "scan x(mu): count=" << cfg.mu_list.size()
                          << " in [" << std::setprecision(12) << mu_min
                          << ", "   << std::setprecision(12) << mu_max << "] (eV)\n";
            }

            std::cout << "Patch size = " << kpatch.size() << "\n";
        }
        rgmpi::barrier();

        // thresholds for quick diagnostics on occ_k_avg (inside only)
        constexpr double occ_near0 = 1e-3;
        constexpr double occ_near1 = 1.0 - 1e-3;

        // -----------------------------
        // 6) scan
        // -----------------------------
        for (double T_K : cfg.temperature_list) {

            const fs::path T_dir = D_dir / (std::string("T") + rgio::tag3(T_K));

            if (rank == 0) {
                std::error_code ec;
                fs::create_directories(T_dir, ec);
                if (ec) {
                    throw std::runtime_error("Failed to create T directory: " +
                                             T_dir.string() + " | " + ec.message());
                }
                std::cout << "\n[T=" << std::fixed << std::setprecision(3) << T_K
                          << "] " << T_dir.string() << "\n";
            }
            rgmpi::barrier();

            // ---------------- mode=doping ----------------
            if (mode == "doping") {

                for (size_t i = 0; i < cfg.filling_list.size(); ++i) {

                    double doping = cfg.doping_list[i];
                    if (std::abs(doping) < 5e-13) doping = 0.0;

                    const double filling = cfg.filling_list[i];

                    // all ranks compute (may have MPI collectives inside)
                    core::GridData fsgrid = model.cal_fermi_patch_from_filling(filling, T_K, kpatch);

                    if (rank != 0) continue;

                    fsgrid.assert_consistent();

                    // inside optional
                    std::vector<unsigned char> inside_fallback;
                    const std::vector<unsigned char>* inPtr = nullptr;
                    try {
                        inPtr = &fsgrid.get<unsigned char>("inside").v;
                        if (inPtr->size() != fsgrid.size())
                            rgmpi::abort_all("calculate_fermi: inside size mismatch");
                    } catch (...) {
                        inside_fallback.assign(fsgrid.size(), (unsigned char)1);
                        inPtr = &inside_fallback;
                    }
                    const auto& inList = *inPtr;

                    const auto& EFv  = fsgrid.get<double>("EF_used").v;
                    const auto& occv = fsgrid.get<double>("occ_k_avg").v;

                    const size_t N = fsgrid.size();
                    if (EFv.size() != N || occv.size() != N)
                        rgmpi::abort_all("calculate_fermi: field size mismatch in fsgrid");

                    long long N_in = 0, N_occ0 = 0, N_occ1 = 0;
                    double occ_min = +std::numeric_limits<double>::infinity();
                    double occ_max = -std::numeric_limits<double>::infinity();
                    double occ_sum = 0.0;

                    for (size_t j = 0; j < N; ++j) {
                        if (!inList[j]) continue;
                        const double occ = occv[j];
                        if (!std::isfinite(occ)) continue;
                        ++N_in;
                        occ_sum += occ;
                        occ_min = std::min(occ_min, occ);
                        occ_max = std::max(occ_max, occ);
                        if (occ < occ_near0) ++N_occ0;
                        if (occ > occ_near1) ++N_occ1;
                    }

                    const double EF_used  = (N > 0) ? EFv[0] : std::numeric_limits<double>::quiet_NaN();
                    const double occ_mean = (N_in > 0) ? (occ_sum / double(N_in))
                                                       : std::numeric_limits<double>::quiet_NaN();

                    std::cout << std::fixed << std::setprecision(4);
                    std::cout << "[x=doping " << doping << "] "
                              << "EF=" << std::setprecision(8) << EF_used
                              << " occ_mean=" << std::setprecision(6) << occ_mean
                              << " occ[min,max]=[" << std::setprecision(6) << occ_min
                              << "," << std::setprecision(6) << occ_max << "]"
                              << " inside=" << N_in
                              << " near0=" << N_occ0
                              << " near1=" << N_occ1
                              << "\n";

                    // folder: dopingXXXX
                    const fs::path x_dir = T_dir / (std::string("doping") + rgio::tag4(doping));
                    {
                        std::error_code ec;
                        fs::create_directories(x_dir, ec);
                        if (ec) throw std::runtime_error("Failed to create x directory: " +
                                                        x_dir.string() + " | " + ec.message());
                    }

                    const fs::path out_path = x_dir / (std::string("fermiPatch_") + run_stamp + ".txt");

                    const std::string header =
                        rgio::write_cal_fermi_header(cfg, st, model, T_K, EF_used, doping, filling);

                    rgio::write_fermi_patch_txt(out_path.string(), fsgrid, header, /*debug_bands=*/true);
                    std::cout << "Wrote: " << out_path.string() << "\n";
                }
            }

            // ---------------- mode=mu ----------------
            else {

                for (size_t i = 0; i < cfg.mu_list.size(); ++i) {

                    double mu = cfg.mu_list[i];
                    if (std::abs(mu) < 5e-13) mu = 0.0;

                    // all ranks compute
                    core::GridData fsgrid = model.cal_fermi_patch_from_mu(mu, T_K, kpatch);

                    if (rank != 0) continue;

                    fsgrid.assert_consistent();

                    // inside optional
                    std::vector<unsigned char> inside_fallback;
                    const std::vector<unsigned char>* inPtr = nullptr;
                    try {
                        inPtr = &fsgrid.get<unsigned char>("inside").v;
                        if (inPtr->size() != fsgrid.size())
                            rgmpi::abort_all("calculate_fermi: inside size mismatch");
                    } catch (...) {
                        inside_fallback.assign(fsgrid.size(), (unsigned char)1);
                        inPtr = &inside_fallback;
                    }
                    const auto& inList = *inPtr;

                    const auto& EFv  = fsgrid.get<double>("EF_used").v;
                    const auto& occv = fsgrid.get<double>("occ_k_avg").v;

                    const size_t N = fsgrid.size();
                    if (EFv.size() != N || occv.size() != N)
                        rgmpi::abort_all("calculate_fermi: field size mismatch in fsgrid");

                    long long N_in = 0, N_occ0 = 0, N_occ1 = 0;
                    double occ_min = +std::numeric_limits<double>::infinity();
                    double occ_max = -std::numeric_limits<double>::infinity();
                    double occ_sum = 0.0;

                    for (size_t j = 0; j < N; ++j) {
                        if (!inList[j]) continue;
                        const double occ = occv[j];
                        if (!std::isfinite(occ)) continue;
                        ++N_in;
                        occ_sum += occ;
                        occ_min = std::min(occ_min, occ);
                        occ_max = std::max(occ_max, occ);
                        if (occ < occ_near0) ++N_occ0;
                        if (occ > occ_near1) ++N_occ1;
                    }

                    // EF from fsgrid (should equal mu), fallback to mu
                    const double EF_used  = (N > 0) ? EFv[0] : mu;

                    // filling = mean_k_in(occ_k_avg)
                    const double filling_used = (N_in > 0)
                        ? (occ_sum / double(N_in))
                        : std::numeric_limits<double>::quiet_NaN();

                    // doping from filling (unit: 1e12 cm^-2)
                    const double doping_used = la::filling_to_doping(filling_used, st.a());

                    std::cout << std::fixed << std::setprecision(6);
                    std::cout << "[x=mu " << mu << "] "
                              << "EF=" << std::setprecision(8) << EF_used
                              << " filling=" << std::setprecision(10) << filling_used
                              << " doping="  << std::setprecision(10) << doping_used
                              << " occ[min,max]=[" << std::setprecision(6) << occ_min
                              << "," << std::setprecision(6) << occ_max << "]"
                              << " inside=" << N_in
                              << " near0=" << N_occ0
                              << " near1=" << N_occ1
                              << "\n";

                    const fs::path x_dir = T_dir / (std::string("mu") + rgio::tag6(mu));
                    {
                        std::error_code ec;
                        fs::create_directories(x_dir, ec);
                        if (ec) throw std::runtime_error("Failed to create x directory: " +
                                                        x_dir.string() + " | " + ec.message());
                    }

                    const fs::path out_path = x_dir / (std::string("fermiPatch_") + run_stamp + ".txt");

                    const std::string header =
                        rgio::write_cal_fermi_header(cfg, st, model, T_K, EF_used, doping_used, filling_used);

                    rgio::write_fermi_patch_txt(out_path.string(), fsgrid, header, /*debug_bands=*/true);
                    std::cout << "Wrote: " << out_path.string() << "\n";
                }
            }
        }

        if (rank == 0) std::cout << "\nDONE.\n";

        rgmpi::finalize();
        return 0;

    } catch (const std::exception& e) {
        try {
            rgmpi::abort_all(std::string("calculate_fermi FAILED: ") + e.what(), 1);
        } catch (...) {}
        rgmpi::finalize();
        return 1;
    }
}
