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
// double doping
// double filling
// double dx
// double dy
// char mesh_type[32]
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
    double doping,
    double filling
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
    const int32_t version = 3;
    const int32_t NkTot   = static_cast<int32_t>(fsgrid.size());
    const int32_t dim     = static_cast<int32_t>(evalsF[0].size());

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
    ofs.write(reinterpret_cast<const char*>(&doping),  sizeof(doping));
    ofs.write(reinterpret_cast<const char*>(&filling), sizeof(filling));

    ofs.write(reinterpret_cast<const char*>(&fsgrid.dx), sizeof(fsgrid.dx));
    ofs.write(reinterpret_cast<const char*>(&fsgrid.dy), sizeof(fsgrid.dy));
    ofs.write(reinterpret_cast<const char*>(mesh_type), sizeof(mesh_type));

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

        core::GridData bzpatch =
            st.generate_BZ_kmesh(
                cfg.bzmesh.Nk,
                cfg.bzmesh.Nk
            );

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

        const fs::path root_dir =
            base_dir
            / (
                std::string("fermi_")
              + cfg.model
              + "_mu_"
              + cfg.output_suffix
            );

        const fs::path D_dir =
            root_dir
            / (std::string("D") + rgio::tag3(cfg.Dfield_eV));

        const auto [Tmin, Tmax] =
            minmax_or_nan(cfg.temperature_list);

        const auto [mu_min, mu_max] =
            minmax_or_nan(cfg.mu_list);

        if (rank == 0) {
            fs::create_directories(D_dir);

            std::cout << "=== calculate_fermi ===\n";
            std::cout << "config     = " << config_file << "\n";
            std::cout << "run_stamp  = " << run_stamp << "\n";
            std::cout << "data_dir   = " << base_dir.string() << "\n";
            std::cout << "out_root   = " << D_dir.string() << "\n";
            std::cout << "model      = " << cfg.model << "\n";
            std::cout << "scan       = mu\n";
            std::cout << "suffix     = " << cfg.output_suffix << "\n";
            std::cout << "layers     = " << para.layer_num << "\n";
            std::cout << "pressure   = " << para.pressure << "\n";
            std::cout << "Dfield     = " << cfg.Dfield_eV << "\n";

            std::cout << "kmesh      = " << cfg.kmesh.type
                      << " Nk=" << cfg.kmesh.Nk
                      << " dk_frac=" << cfg.kmesh.dk_frac
                      << " mesh_type=" << kpatch.mesh_type
                      << " size=" << kpatch.size()
                      << "\n";

            std::cout << "bzmesh     = " << cfg.bzmesh.type
                      << " Nk=" << cfg.bzmesh.Nk
                      << " mesh_type=" << bzpatch.mesh_type
                      << " size=" << bzpatch.size()
                      << "\n";

            std::cout << "MPI        = " << nprocs << " ranks\n";

            std::cout << "scan T     = count="
                      << cfg.temperature_list.size()
                      << " in [" << Tmin << ", " << Tmax << "]\n";

            std::cout << "scan mu    = count="
                      << cfg.mu_list.size()
                      << " in [" << mu_min << ", " << mu_max << "]\n";
        }

        rgmpi::barrier();

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
                // 1. Full BZ filling/doping
                // ====================================================
                const double filling_BZ =
                    model.find_filling_from_EF(
                        mu,
                        T_K,
                        bzpatch
                    );

                const double doping_BZ =
                    la::filling_to_doping(
                        filling_BZ,
                        st.a()
                    );

                // ====================================================
                // 2. Local/truncated fermi patch with eigenvectors
                // ====================================================
                core::GridData fsgrid =
                    model.cal_fermi_patch_from_mu(
                        mu,
                        T_K,
                        kpatch
                    );

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
                          << " filling_BZ="
                          << std::setprecision(10)
                          << filling_BZ
                          << " doping_BZ="
                          << std::setprecision(10)
                          << doping_BZ
                          << "\n";

                const fs::path mu_dir =
                    T_dir
                    / (std::string("mu") + rgio::tag6(mu));

                fs::create_directories(mu_dir);

                const fs::path out_path =
                    mu_dir
                    / (std::string("fermiPatch_") + run_stamp + ".bin");

                write_fermi_patch_bin(
                    out_path.string(),
                    fsgrid,
                    EF_used,
                    T_K,
                    doping_BZ,
                    filling_BZ
                );

                std::cout << "Wrote: "
                          << out_path.string()
                          << "\n";
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