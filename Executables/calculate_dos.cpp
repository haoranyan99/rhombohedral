// File: Executables/calculate_dos.cpp
#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <memory>
#include <vector>
#include <filesystem>

#ifdef USE_MPI
  #include <mpi.h>
#endif

#include "Util/MPI/mpi.h"
#include "Util/Config/config_reader.h"
#include "Util/IO/rg_io.h"
#include "LinearAlgebra/MathFunctions.h"

#include "PhysStruct/RG_Structure.h"
#include "Source/Models/RG_SKModel.h"
#include "Source/Models/RG_KPModel.h"

namespace fs = std::filesystem;

static core::GridData make_kmesh(
    const rg::RG_Structure& st,
    const config::MeshCfg& mesh
) {
    if (mesh.type == "hex" || mesh.type == "localK_hex") {
        return st.generate_localK_kmesh_hex_b1b2(mesh.Nk, mesh.dk_frac);
    }

    if (mesh.type == "square" || mesh.type == "localK_square") {
        return st.generate_localK_kmesh_square(mesh.Nk, mesh.dk_frac);
    }

    if (mesh.type == "b1b2" || mesh.type == "localK_b1b2") {
        return st.generate_localK_kmesh_b1b2(mesh.Nk, mesh.dk_frac);
    }

    if (mesh.type == "BZ") {
        return st.generate_BZ_kmesh(mesh.Nk, mesh.Nk);
    }

    rgmpi::abort_all("calculate_dos: unsupported kmesh.type = " + mesh.type);
    return core::GridData{};
}

int main(int argc, char** argv) {
#ifdef USE_MPI
    MPI_Init(&argc, &argv);
#endif

    int rank = 0, nprocs = 1;
    rgmpi::rank_size(rank, nprocs);

    try {
        const std::string config_file =
            (argc >= 2 ? std::string(argv[1]) : "config.json");

        // ============================================================
        // read config / para
        // ============================================================
        const config::CalDOSConfig cfg =
            config::read_cal_dos_config(config_file);

        const config::RGPara para =
            config::read_rg_para(cfg.para_file);

        const fs::path base_dir =
            cfg.data_dir.empty() ? fs::path("data")
                                 : fs::path(cfg.data_dir);

        const std::string run_stamp =
            rgio::make_time_stamp("dos");

        if (rank == 0) {
            std::cout << "=== calculate_dos ===\n";
            std::cout << "config   = " << config_file << "\n";
            std::cout << "run_stamp= " << run_stamp << "\n";
            std::cout << "base_dir = " << base_dir.string() << "\n";
            std::cout << "model    = " << cfg.model << "\n";
            std::cout << "suffix   = " << cfg.output_suffix << "\n";
#ifdef USE_MPI
            std::cout << "MPI      = " << nprocs << " ranks\n";
#endif
        }

        // ============================================================
        // build structure
        // ============================================================
        rg::RG_Structure st(
            para.layer_num,
            para.pressure,
            para.a0,
            para.d0,
            para.vacuum
        );

        // ============================================================
        // build model
        // ============================================================
        std::unique_ptr<rg::RG_ModelBase> model_ptr;
        if (cfg.model == "sk") {
            model_ptr = std::make_unique<rg::RG_SKModel>(st, cfg.para_file);
        } else if (cfg.model == "kp") {
            model_ptr = std::make_unique<rg::RG_KPModel>(st);
        } else {
            rgmpi::abort_all("calculate_dos: model must be \"sk\" or \"kp\"");
        }

        rg::RG_ModelBase& model = *model_ptr;

        if (cfg.model == "sk") {
            if (auto* sk = dynamic_cast<rg::RG_SKModel*>(model_ptr.get())) {
                sk->ensure_hoppings_cached();
            }
        }

        // ============================================================
        // output root: data/dos_{model}_{suffix}/D...
        // ============================================================
        const fs::path out_root =
            base_dir
            / (
                std::string("dos_")
              + cfg.model
              + "_"
              + cfg.output_suffix
            );

        if (rank == 0) {
            std::error_code ec;
            fs::create_directories(out_root, ec);
            if (ec) {
                throw std::runtime_error("Failed to create output directory: " +
                                         out_root.string() + " | " + ec.message());
            }
            std::cout << "out_root = " << out_root.string() << "\n";
        }
        rgmpi::barrier();

        // ============================================================
        // generate k-mesh (GridData)
        // ============================================================
        core::GridData kmesh =
            make_kmesh(st, cfg.kmesh);

        if (rank == 0) {
            std::cout << "kmesh    = " << cfg.kmesh.type
                      << " Nk=" << cfg.kmesh.Nk
                      << " dk_frac=" << cfg.kmesh.dk_frac
                      << " (size=" << kmesh.size() << ")\n";
            std::cout << "DOS(E)   = [" << cfg.dos.e_low << ", " << cfg.dos.e_high
                      << "] num_e=" << cfg.dos.num_e
                      << " eta=" << cfg.dos.eta
                      << " T_K=" << cfg.T_K << "\n";
        }

        for (double Dfield_eV : cfg.Dfield_list_eV) {
            model.set_Dfield(Dfield_eV);

            if (rank == 0) {
                std::cout << "\n[D="
                          << std::fixed
                          << std::setprecision(6)
                          << Dfield_eV
                          << "] calculating DOS\n";
            }

            core::SeriesData dosS =
                model.cal_dos_gaussian(
                    kmesh,
                    cfg.dos.e_low,
                    cfg.dos.e_high,
                    cfg.dos.num_e,
                    cfg.dos.eta,
                    cfg.T_K,
                    /*enforce_hermitian=*/true
                );

            if (rank == 0) {
                const fs::path D_dir =
                    out_root
                    / (std::string("D") + rgio::tag3(Dfield_eV));

                fs::create_directories(D_dir);

                const fs::path out_path =
                    D_dir / (run_stamp + ".txt");

                std::ostringstream hdr;
                st.write_info(hdr, "# ");
                model.write_info(hdr, "# ");

                hdr << "# --- calculate_dos ---\n";
                hdr << "# config_file = " << config_file << "\n";
                hdr << "# model = " << cfg.model << "\n";
                hdr << "# Dfield_eV = " << Dfield_eV << "\n";
                hdr << "# T_K = " << cfg.T_K << " K\n";
                hdr << "# output_suffix = " << cfg.output_suffix << "\n";
                hdr << "# kmesh = " << cfg.kmesh.type
                    << " Nk=" << cfg.kmesh.Nk
                    << " dk_frac=" << cfg.kmesh.dk_frac << "\n";
                hdr << "# E_scan = [" << cfg.dos.e_low << ", " << cfg.dos.e_high
                    << "], num_e=" << cfg.dos.num_e
                    << ", eta=" << cfg.dos.eta << "\n";
#ifdef USE_MPI
                hdr << "# mpi_ranks = " << nprocs << "\n";
#endif

                rgio::write_dos_txt(out_path.string(), dosS, hdr.str());

                std::cout << "Wrote: " << out_path.string() << "\n";
            }
        }

        if (rank == 0) {
            std::cout << "\nDONE.\n";
        }

        rgmpi::finalize();
        return 0;

    } catch (const std::exception& e) {
        try {
            rgmpi::abort_all(std::string("calculate_dos FAILED: ") + e.what(), 1);
        } catch (...) {}
        rgmpi::finalize();
        return 1;
    }
}
