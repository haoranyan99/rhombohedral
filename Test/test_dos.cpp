// File: Test/test_dos.cpp
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

static std::string tag_fixed(double x, int prec) {
    std::ostringstream oss;
    oss.setf(std::ios::fixed);
    oss << std::setprecision(prec) << x;
    return oss.str();
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
        const config::TestDOSConfig cfg =
            config::read_test_dos_config(config_file);

        const config::RGPara para =
            config::read_rg_para(cfg.para_file);

        const std::filesystem::path base_dir =
            cfg.data_dir.empty() ? std::filesystem::path("data")
                                 : std::filesystem::path(cfg.data_dir);

        if (rank == 0) {
            std::cout << "=== test_RG_dos ===\n";
            std::cout << "config   = " << config_file << "\n";
            std::cout << "base_dir = " << base_dir.string() << "\n";
            std::cout << "model    = " << cfg.model << "\n";
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
            model_ptr = std::make_unique<rg::RG_SKModel>(st);
        } else if (cfg.model == "kp") {
            model_ptr = std::make_unique<rg::RG_KPModel>(st);
        } else {
            rgmpi::abort_all("test_RG_dos: model must be \"sk\" or \"kp\"");
        }

        rg::RG_ModelBase& model = *model_ptr;
        model.set_Dfield(cfg.Dfield_eV);

        // (optional) SK cache hoppings
        if (cfg.model == "sk") {
            if (auto* sk = dynamic_cast<rg::RG_SKModel*>(model_ptr.get())) {
                sk->ensure_hoppings_cached();
            }
        }

        // ============================================================
        // output subdir: data/dos_{model}_D{Dfield}eV/
        // ============================================================
        const std::string subdir_name =
            std::string("dos_") + cfg.model + "_D" + tag_fixed(cfg.Dfield_eV, 3) + "eV";
        const std::filesystem::path out_subdir = base_dir / subdir_name;

        if (rank == 0) {
            std::error_code ec;
            std::filesystem::create_directories(out_subdir, ec);
            if (ec) {
                throw std::runtime_error("Failed to create output directory: " +
                                         out_subdir.string() + " | " + ec.message());
            }
            std::cout << "out_subdir = " << out_subdir.string() << "\n";
        }
        rgmpi::barrier();

        // ============================================================
        // generate k-mesh (GridData)
        // ============================================================
        core::GridData kmesh;
        if (cfg.kmesh.type == "localK_hex") {
            kmesh = st.generate_localK_kmesh_hex_b1b2(cfg.kmesh.Nk, cfg.kmesh.dk_frac);
        } else if (cfg.kmesh.type == "BZ") {
            kmesh = st.generate_BZ_kmesh(cfg.kmesh.Nk, cfg.kmesh.Nk);
        } else {
            rgmpi::abort_all("test_RG_dos: unsupported kmesh.type = " + cfg.kmesh.type);
        }

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

        // ============================================================
        // DOS Gaussian (also filling/doping columns if model provides them)
        // ============================================================
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

        // ============================================================
        // output (rank0 only)
        // ============================================================
        if (rank == 0) {
            // filename
            const std::string run_stamp = rgio::make_time_stamp("dos_");
            const std::string out_name = run_stamp + ".txt";

            // header = st.write_info + model.write_info + run info
            std::ostringstream hdr;
            st.write_info(hdr, "# ");
            model.write_info(hdr, "# ");

            hdr << "# --- dos run ---\n";
            hdr << "# config_file = " << config_file << "\n";
            hdr << "# model = " << cfg.model << "\n";
            hdr << "# Dfield_eV = " << cfg.Dfield_eV << "\n";
            hdr << "# T_K = " << cfg.T_K << " K\n";
            hdr << "# kmesh = " << cfg.kmesh.type
                << " Nk=" << cfg.kmesh.Nk
                << " dk_frac=" << cfg.kmesh.dk_frac << "\n";
            hdr << "# E_scan = [" << cfg.dos.e_low << ", " << cfg.dos.e_high
                << "], num_e=" << cfg.dos.num_e
                << ", eta=" << cfg.dos.eta << "\n";
#ifdef USE_MPI
            hdr << "# mpi_ranks = " << nprocs << "\n";
#endif

            rgio::write_dos_txt(out_name, dosS, hdr.str());

            std::cout << "Wrote:\n  " << out_name << "\nDONE.\n";
        }

        rgmpi::finalize();
        return 0;

    } catch (const std::exception& e) {
        try {
            rgmpi::abort_all(std::string("test_RG_dos FAILED: ") + e.what(), 1);
        } catch (...) {}
        rgmpi::finalize();
        return 1;
    }
}
