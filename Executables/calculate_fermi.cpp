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
#include <fstream>
#include <cstdint>
#include <cstring>

#ifdef USE_MPI
  #include <mpi.h>
#endif

#include "Util/MPI/mpi.h"
#include "PhysStruct/RG_Structure.h"
#include "Source/Models/RG_ModelBase.h"
#include "Source/Models/RG_SKModel.h"
#include "Source/Models/RG_KPModel.h"
#include "Util/Config/config_reader.h"
#include "Util/IO/rg_io.h"
#include "LinearAlgebra/MathFunctions.h"

namespace fs = std::filesystem;

static std::pair<double, double> minmax_or_nan(const std::vector<double>& v)
{
    if (v.empty()) {
        const double NaN = std::numeric_limits<double>::quiet_NaN();
        return {NaN, NaN};
    }

    auto [mn_it, mx_it] = std::minmax_element(v.begin(), v.end());
    return {*mn_it, *mx_it};
}

static std::string compact_double_tag(double x)
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

static std::string auto_fermi_output_suffix(
    const config::CalFermiConfig& cfg,
    double B_T
)
{
    std::ostringstream oss;

    oss << cfg.kmesh.Nk
        << "_"
        << cfg.kmesh.type
        << "_"
        << compact_double_tag(cfg.kmesh.dk_frac)
        << "_D"
        << compact_double_tag(cfg.Dfield_eV);

    if (std::abs(B_T) > 0.0) {
        oss << "_B"
            << compact_double_tag(B_T)
            << "T";
    }

    return oss.str();
}

static std::string fermi_output_suffix(
    const config::CalFermiConfig& cfg,
    double B_T
) {
    std::string suffix =
        auto_fermi_output_suffix(cfg, B_T);

    if (
        !cfg.output_suffix.empty()
     && cfg.output_suffix != "auto"
    ) {
        suffix += "_" + cfg.output_suffix;
    }

    return suffix;
}

static core::GridData make_kmesh(
    const rg::RG_Structure& st,
    const config::MeshCfg& mesh
) {
    if (mesh.type == "hex") {
        return st.generate_localK_kmesh_hex_b1b2(
            mesh.Nk,
            mesh.dk_frac
        );
    }

    if (mesh.type == "square") {
        return st.generate_localK_kmesh_square(
            mesh.Nk,
            mesh.dk_frac
        );
    }

    if (mesh.type == "b1b2") {
        return st.generate_localK_kmesh_b1b2(
            mesh.Nk,
            mesh.dk_frac
        );
    }

    return st.generate_BZ_kmesh(
        mesh.Nk,
        mesh.Nk
    );
}

static double patch_area_fraction(
    const core::GridData& kpatch
) {
    if (kpatch.size() == 0) {
        return 0.0;
    }

    const double frac =
        static_cast<double>(kpatch.size())
      * kpatch.dx
      * kpatch.dy;

    if (!std::isfinite(frac) || frac < 0.0) {
        return 0.0;
    }

    return std::min(frac, 1.0);
}

static double patch_filling(
    const core::GridData& fsgrid
) {
    fsgrid.assert_consistent();

    const auto& occAvg =
        fsgrid.get<double>("occ_k_avg").v;

    if (occAvg.empty()) {
        return 0.5;
    }

    double acc = 0.0;
    int n = 0;

    for (double x : occAvg) {
        if (!std::isfinite(x)) {
            continue;
        }

        acc += x;
        ++n;
    }

    if (n == 0) {
        return 0.5;
    }

    return acc / static_cast<double>(n);
}

static double total_filling_from_local_patch(
    double filling_patch,
    double area_fraction
) {
    return 0.5
         + area_fraction
         * (filling_patch - 0.5);
}

// ============================================================
// Binary writer for fermi patch
//
// Header:
// int32 magic
// int32 version
// int32 NkTot
// int32 dim
// double EF
// double T_K
// double doping_total
// double filling_total
// double dx
// double dy
// char mesh_type[32]
// int32 spin_sign       (v5: +1 up, -1 down, 0 no spin split)
// double doping_spin    (v5)
// double filling_spin   (v5)
//
// Per k:
// int32 iq
// int32 jq
// double kx
// double ky
// double occ_avg
// double evals[dim]
// double occ_band[dim]
// double evec_re[dim*dim]
// double evec_im[dim*dim]
// ============================================================
static void write_fermi_patch_bin(
    const std::string& out_path,
    const core::GridData& fsgrid,
    double EF,
    double T_K,
    double doping_total,
    double filling_total,
    int spin_sign,
    double doping_spin,
    double filling_spin
) {
    fsgrid.assert_consistent();

    std::ofstream ofs(out_path, std::ios::binary);
    if (!ofs) {
        throw std::runtime_error("Cannot open output bin file: " + out_path);
    }

    const auto& kvecF   = fsgrid.get<rg::Vec2>("kvec").v;
    const auto& occAvgF = fsgrid.get<double>("occ_k_avg").v;
    const auto& evalsF  = fsgrid.get<std::vector<double>>("evals").v;
    const auto& occbF   = fsgrid.get<std::vector<double>>("occ_band").v;
    const auto& evecReF = fsgrid.get<std::vector<double>>("evec_re").v;
    const auto& evecImF = fsgrid.get<std::vector<double>>("evec_im").v;

    if (fsgrid.size() == 0) {
        throw std::runtime_error("write_fermi_patch_bin: empty fsgrid");
    }

    const int32_t magic   = 20260510;
    const int32_t version = 5;
    const int32_t NkTot   = static_cast<int32_t>(fsgrid.size());
    const int32_t dim     = static_cast<int32_t>(evalsF[0].size());
    const int32_t spin_sign_i32 = static_cast<int32_t>(spin_sign);

    char mesh_type[32] = {};
    std::snprintf(
        mesh_type,
        sizeof(mesh_type),
        "%s",
        fsgrid.mesh_type.c_str()
    );

    ofs.write(reinterpret_cast<const char*>(&magic),   sizeof(magic));
    ofs.write(reinterpret_cast<const char*>(&version), sizeof(version));
    ofs.write(reinterpret_cast<const char*>(&NkTot),   sizeof(NkTot));
    ofs.write(reinterpret_cast<const char*>(&dim),     sizeof(dim));

    ofs.write(reinterpret_cast<const char*>(&EF),      sizeof(EF));
    ofs.write(reinterpret_cast<const char*>(&T_K),     sizeof(T_K));
    ofs.write(reinterpret_cast<const char*>(&doping_total),  sizeof(doping_total));
    ofs.write(reinterpret_cast<const char*>(&filling_total), sizeof(filling_total));

    ofs.write(reinterpret_cast<const char*>(&fsgrid.dx), sizeof(fsgrid.dx));
    ofs.write(reinterpret_cast<const char*>(&fsgrid.dy), sizeof(fsgrid.dy));
    ofs.write(reinterpret_cast<const char*>(mesh_type), sizeof(mesh_type));
    ofs.write(reinterpret_cast<const char*>(&spin_sign_i32), sizeof(spin_sign_i32));
    ofs.write(reinterpret_cast<const char*>(&doping_spin), sizeof(doping_spin));
    ofs.write(reinterpret_cast<const char*>(&filling_spin), sizeof(filling_spin));

    for (size_t ik = 0; ik < fsgrid.size(); ++ik) {
        const int32_t iq = static_cast<int32_t>(fsgrid.iq[ik]);
        const int32_t jq = static_cast<int32_t>(fsgrid.jq[ik]);

        const double kx = kvecF[ik].x();
        const double ky = kvecF[ik].y();
        const double occ_avg = occAvgF[ik];

        if ((int)evalsF[ik].size()  != dim ||
            (int)occbF[ik].size()   != dim ||
            (int)evecReF[ik].size() != dim * dim ||
            (int)evecImF[ik].size() != dim * dim)
        {
            throw std::runtime_error(
                "write_fermi_patch_bin: inconsistent per-k vector size"
            );
        }

        ofs.write(reinterpret_cast<const char*>(&iq), sizeof(iq));
        ofs.write(reinterpret_cast<const char*>(&jq), sizeof(jq));
        ofs.write(reinterpret_cast<const char*>(&kx), sizeof(kx));
        ofs.write(reinterpret_cast<const char*>(&ky), sizeof(ky));
        ofs.write(reinterpret_cast<const char*>(&occ_avg), sizeof(occ_avg));

        ofs.write(
            reinterpret_cast<const char*>(evalsF[ik].data()),
            sizeof(double) * dim
        );

        ofs.write(
            reinterpret_cast<const char*>(occbF[ik].data()),
            sizeof(double) * dim
        );

        ofs.write(
            reinterpret_cast<const char*>(evecReF[ik].data()),
            sizeof(double) * dim * dim
        );

        ofs.write(
            reinterpret_cast<const char*>(evecImF[ik].data()),
            sizeof(double) * dim * dim
        );

    }
}

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

        const std::string run_stamp =
            rgio::make_time_stamp("");

        const auto cfg =
            config::read_cal_fermi_config(config_file);

        const auto para =
            config::read_rg_para(cfg.para_file);

        rg::RG_Structure st(
            para.layer_num,
            para.pressure,
            para.a0,
            para.d0,
            para.vacuum
        );

        core::GridData kpatch =
            make_kmesh(st, cfg.kmesh);

        const double kmesh_area_fraction =
            patch_area_fraction(kpatch);

        std::unique_ptr<rg::RG_ModelBase> model_ptr;

        if (cfg.model == "sk") {
            model_ptr =
                std::make_unique<rg::RG_SKModel>(
                    st,
                    cfg.para_file
                );
        } else {
            model_ptr =
                std::make_unique<rg::RG_KPModel>(st);
        }

        rg::RG_ModelBase& model = *model_ptr;

        if (cfg.model == "sk") {
            auto* sk =
                dynamic_cast<rg::RG_SKModel*>(model_ptr.get());

            sk->ensure_hoppings_cached();
        }

        model.set_Dfield(cfg.Dfield_eV);

        const fs::path base_dir =
            cfg.data_dir.empty()
            ? fs::path("data")
            : fs::path(cfg.data_dir);

        const auto [Tmin, Tmax] =
            minmax_or_nan(cfg.temperature_list);

        const auto [mu_min, mu_max] =
            minmax_or_nan(cfg.mu_list);

        const auto [B_min, B_max] =
            minmax_or_nan(cfg.Bfield_list_T);

        if (rank == 0) {
            std::cout << "=== calculate_fermi ===\n";
            std::cout << "config     = " << config_file << "\n";
            std::cout << "run_stamp  = " << run_stamp << "\n";
            std::cout << "data_dir   = " << base_dir.string() << "\n";
            std::cout << "model      = " << cfg.model << "\n";
            std::cout << "scan       = mu\n";
            std::cout << "layers     = " << para.layer_num << "\n";
            std::cout << "pressure   = " << para.pressure << "\n";
            std::cout << "Dfield     = " << cfg.Dfield_eV << "\n";
            std::cout << "g_factor   = " << cfg.g_factor << "\n";
            std::cout << "orb_dk     = " << cfg.orbital_derivative_dk << "\n";

            std::cout << "kmesh      = " << cfg.kmesh.type
                      << " Nk=" << cfg.kmesh.Nk
                      << " dk_frac=" << cfg.kmesh.dk_frac
                      << " mesh_type=" << kpatch.mesh_type
                      << " size=" << kpatch.size()
                      << "\n";

            std::cout << "doping     = local-kmesh approximation\n";
            std::cout << "kmesh_frac = " << std::setprecision(12)
                      << kmesh_area_fraction
                      << " (outside patch fixed at filling=0.5, doping=0)\n";

            std::cout << "MPI        = " << nprocs << " ranks\n";

            std::cout << "scan T     = count="
                      << cfg.temperature_list.size()
                      << " in [" << Tmin << ", " << Tmax << "]\n";

            std::cout << "scan mu    = count="
                      << cfg.mu_list.size()
                      << " in [" << mu_min << ", " << mu_max << "]\n";

            std::cout << "scan B_T   = count="
                      << cfg.Bfield_list_T.size()
                      << " in [" << B_min << ", " << B_max << "]\n";
        }

        rgmpi::barrier();

        for (double B_T : cfg.Bfield_list_T) {
            const std::string output_suffix =
                fermi_output_suffix(cfg, B_T);

            const fs::path root_dir =
                base_dir
                / (
                    std::string("fermi_")
                  + cfg.model
                  + "_mu_"
                  + output_suffix
                  + run_stamp
                );

            const fs::path D_dir =
                root_dir
                / (std::string("D") + rgio::tag3(cfg.Dfield_eV));

            if (rank == 0) {
                fs::create_directories(D_dir);
                std::cout << "\n[B="
                          << compact_double_tag(B_T)
                          << " T] suffix="
                          << output_suffix
                          << " out_root="
                          << D_dir.string()
                          << "\n";
            }

            rgmpi::barrier();

            const bool use_magnetic_splitting =
                std::abs(B_T) > 0.0;

            for (double T_K : cfg.temperature_list) {
                const fs::path T_dir =
                    D_dir / (std::string("T") + rgio::tag3(T_K));

                if (rank == 0) {
                    fs::create_directories(T_dir);

                    std::cout << "\n[T="
                              << std::fixed
                              << std::setprecision(3)
                              << T_K
                              << "] "
                              << T_dir.string()
                              << "\n";
                }

                rgmpi::barrier();

                for (double mu_raw : cfg.mu_list) {
                    double mu = mu_raw;

                    if (std::abs(mu) < 5e-13) {
                        mu = 0.0;
                    }

                    // ====================================================
                    // 1. Local/truncated fermi patch with eigenvectors
                    // ====================================================
                    core::GridData fsgrid;
                    core::GridData fsgrid_up;
                    core::GridData fsgrid_dn;

                    if (use_magnetic_splitting) {
                        fsgrid_up =
                            model.cal_fermi_patch_from_mu_magnetic(
                                mu,
                                T_K,
                                kpatch,
                                B_T,
                                cfg.g_factor,
                                cfg.orbital_derivative_dk,
                                +1
                            );

                        fsgrid_dn =
                            model.cal_fermi_patch_from_mu_magnetic(
                                mu,
                                T_K,
                                kpatch,
                                B_T,
                                cfg.g_factor,
                                cfg.orbital_derivative_dk,
                                -1
                            );
                    } else {
                        fsgrid =
                            model.cal_fermi_patch_from_mu(
                                mu,
                                T_K,
                                kpatch
                            );
                    }

                    // ====================================================
                    // 2. Approximate total filling/doping from local patch
                    //    Outside the patch is fixed at half filling.
                    // ====================================================
                    double filling_patch = 0.5;
                    double filling_up_patch = 0.5;
                    double filling_dn_patch = 0.5;

                    double filling_total = 0.5;
                    double doping_total = 0.0;
                    double filling_up_total = 0.5;
                    double filling_dn_total = 0.5;
                    double doping_up_total = 0.0;
                    double doping_dn_total = 0.0;

                    if (use_magnetic_splitting) {
                        filling_up_patch =
                            patch_filling(fsgrid_up);

                        filling_dn_patch =
                            patch_filling(fsgrid_dn);

                        filling_up_total =
                            total_filling_from_local_patch(
                                filling_up_patch,
                                kmesh_area_fraction
                            );

                        filling_dn_total =
                            total_filling_from_local_patch(
                                filling_dn_patch,
                                kmesh_area_fraction
                            );

                        filling_total =
                            0.5 * (filling_up_total + filling_dn_total);

                        doping_up_total =
                            la::filling_to_doping(
                                filling_up_total,
                                st.a()
                            );

                        doping_dn_total =
                            la::filling_to_doping(
                                filling_dn_total,
                                st.a()
                            );

                        doping_total =
                            la::filling_to_doping(
                                filling_total,
                                st.a()
                            );
                    } else {
                        filling_patch =
                            patch_filling(fsgrid);

                        filling_total =
                            total_filling_from_local_patch(
                                filling_patch,
                                kmesh_area_fraction
                            );

                        doping_total =
                            la::filling_to_doping(
                                filling_total,
                                st.a()
                            );
                    }

                    if (rank != 0) {
                        continue;
                    }

                    const double EF_used = mu;

                    std::cout << std::fixed
                              << std::setprecision(6);

                    std::cout << "[mu "
                              << mu
                              << "] EF="
                              << std::setprecision(8)
                              << EF_used
                              << " filling_total="
                              << std::setprecision(10)
                              << filling_total
                              << " doping_total="
                              << std::setprecision(10)
                              << doping_total
                              << " filling_patch="
                              << std::setprecision(10)
                              << (
                                     use_magnetic_splitting
                                   ? 0.5 * (filling_up_patch + filling_dn_patch)
                                   : filling_patch
                                 );

                    if (use_magnetic_splitting) {
                        std::cout
                                  << " filling_up="
                                  << filling_up_total
                                  << " filling_down="
                                  << filling_dn_total
                                  << " filling_up_patch="
                                  << filling_up_patch
                                  << " filling_down_patch="
                                  << filling_dn_patch
                                  << " doping_up="
                                  << doping_up_total
                                  << " doping_down="
                                  << doping_dn_total;
                    }

                    std::cout << "\n";

                    const std::string mu_tag =
                        std::string("mu") + rgio::tag6(mu);

                    const fs::path mu_dir =
                        T_dir / mu_tag;
                    fs::create_directories(mu_dir);

                    if (use_magnetic_splitting) {
                        const fs::path out_path_up =
                            mu_dir
                            / "fermi_spin_up_patch.bin";

                        const fs::path out_path_dn =
                            mu_dir
                            / "fermi_spin_down_patch.bin";

                        write_fermi_patch_bin(
                            out_path_up.string(),
                            fsgrid_up,
                            EF_used,
                            T_K,
                            doping_total,
                            filling_total,
                            +1,
                            doping_up_total,
                            filling_up_total
                        );

                        write_fermi_patch_bin(
                            out_path_dn.string(),
                            fsgrid_dn,
                            EF_used,
                            T_K,
                            doping_total,
                            filling_total,
                            -1,
                            doping_dn_total,
                            filling_dn_total
                        );

                        std::cout << "Wrote: "
                                  << out_path_up.string()
                                  << "\n";

                        std::cout << "Wrote: "
                                  << out_path_dn.string()
                                  << "\n";
                    } else {
                        const fs::path out_path =
                            mu_dir
                            / "fermiPatch.bin";

                        write_fermi_patch_bin(
                            out_path.string(),
                            fsgrid,
                            EF_used,
                            T_K,
                            doping_total,
                            filling_total,
                            0,
                            doping_total,
                            filling_total
                        );

                        std::cout << "Wrote: "
                                  << out_path.string()
                                  << "\n";
                    }
                }
            }
        }

        if (rank == 0) {
            std::cout << "\nDONE.\n";
        }

        rgmpi::finalize();
        return 0;

    } catch (const std::exception& e) {
        try {
            rgmpi::abort_all(
                std::string("calculate_fermi FAILED: ") + e.what(),
                1
            );
        } catch (...) {}

        rgmpi::finalize();
        return 1;
    }
}
