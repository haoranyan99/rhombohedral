// File: Test/test_bareband.cpp
//
// Bare band structure generator (one global timestamp for the whole run)
//
// Output layout (NO unit suffix in folder names):
//   {data_dir}/band_{suffix}/D{Dfield}/GMKG_{STAMP}.txt
//   {data_dir}/band_{suffix}/D{Dfield}/localK_MKKp_{STAMP}.txt
//
// - STAMP is generated ONCE by rgio::make_time_stamp("") and reused for ALL data files.
// - Directory names do NOT contain STAMP.
//
// Config block (ONLY):
//   tasks.bareband
// including:
//   "output_suffix": "1"
//
// Notes:
// - Only rank 0 writes files.
// - All ranks should call model.bands_along_path() in the same order.

#include <Eigen/Dense>

#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <memory>
#include <limits>
#include <filesystem>
#include <fstream>
#include <algorithm>

#ifdef USE_MPI
  #include <mpi.h>
#endif

#include "Util/MPI/mpi.h"                 // rgmpi helpers
#include "PhysStruct/RG_Structure.h"
#include "Source/Models/RG_SKModel.h"
#include "Source/Models/RG_KPModel.h"
#include "Util/Config/config_reader.h"    // config::read_test_bareband_config + read_rg_para
#include "Util/IO/rg_io.h"                // rgio::write_bands_txt + make_time_stamp

static std::string tag_double(double x, int prec) {
    std::ostringstream oss;
    oss.setf(std::ios::fixed);
    oss << std::setprecision(prec) << x;
    return oss.str();
}

// Write header (truncate) then later rgio::write_bands_txt will append data.
static void write_header_only(const std::string& out_path,
                              const std::string& model_name,
                              const std::string& config_file,
                              const std::string& run_stamp,
                              const std::string& output_suffix,
                              const rg::RG_Structure& st,
                              const rg::RG_ModelBase& model,
                              const core::PathData& path,
                              const std::string& path_name,
                              int Nk_seg,
                              double frac_local /*ignored for GMKG*/)
{
    std::ofstream fout(out_path); // truncate
    if (!fout)
        throw std::runtime_error("test_bareband: cannot write " + out_path);

    fout << "# RG bare band structure\n";
    fout << "# k in 1/Angstrom, E in eV\n";
    fout << "# --- test_bareband ---\n";
    fout << "# config_file = " << config_file << "\n";
    fout << "# run_stamp   = " << run_stamp << "\n";
    fout << "# output_suffix = " << output_suffix << "\n";
    fout << "# model = " << model_name << "\n";
    fout << "# path  = " << path_name << "\n";
    fout << "# Nk_seg = " << Nk_seg << "\n";
    if (path_name != "GMKG") {
        fout << "# frac_local = " << std::setprecision(12) << frac_local << "\n";
    }

    // xticks (optional)
    if (!path.xtick_pos.empty() && !path.xtick_lab.empty()) {
        fout << "# xticks:";
        const size_t n = std::min(path.xtick_pos.size(), path.xtick_lab.size());
        for (size_t i = 0; i < n; ++i) {
            fout << " (" << path.xtick_lab[i] << "," << path.xtick_pos[i] << ")";
        }
        fout << "\n";
    }

    // structure info first
    st.write_info(fout, "# ");

    // model info next
    model.write_info(fout, "# ");

    fout << "# ---- data ----\n";
}

int main(int argc, char** argv) {
    try {

#ifdef USE_MPI
        MPI_Init(&argc, &argv);
#endif

        int rank = 0, nprocs = 1;
        rgmpi::rank_size(rank, nprocs);

        const std::string config_file = (argc > 1) ? std::string(argv[1]) : "config.json";

        // ONE global timestamp for the whole run (used ONLY in data filenames)
        const std::string run_stamp = rgio::make_time_stamp("");

        // -----------------------------
        // 1) read config + parameters
        // -----------------------------
        const auto cfg  = config::read_test_bareband_config(config_file); // reads tasks.bareband
        const auto para = config::read_rg_para(cfg.para_file);

        // -----------------------------
        // 1.5) sanity
        // -----------------------------
        if (cfg.GMKG.Nk_seg <= 0)
            rgmpi::abort_all("test_bareband: GMKG.Nk_seg must be > 0");
        if (cfg.localK_MKKp.Nk_seg <= 0)
            rgmpi::abort_all("test_bareband: localK_MKKp.Nk_seg must be > 0");
        if (!(cfg.localK_MKKp.frac_local > 0.0 && cfg.localK_MKKp.frac_local <= 1.0))
            rgmpi::abort_all("test_bareband: localK_MKKp.frac_local must be in (0,1]");

        // -----------------------------
        // 2) build structure
        // -----------------------------
        const double pressure = 0.0; // keep consistent with your tests
        rg::RG_Structure st(para.layer_num, pressure, para.a0, para.d0, para.vacuum);

        // -----------------------------
        // 3) build model
        // -----------------------------
        std::unique_ptr<rg::RG_ModelBase> model_ptr;
        if (cfg.model == "sk") {
            model_ptr = std::make_unique<rg::RG_SKModel>(st);
        } else if (cfg.model == "kp") {
            model_ptr = std::make_unique<rg::RG_KPModel>(st);
        } else {
            rgmpi::abort_all("test_bareband: model must be \"sk\" or \"kp\"");
        }

        rg::RG_ModelBase& model = *model_ptr;

        if (cfg.model == "sk") {
            if (auto* sk = dynamic_cast<rg::RG_SKModel*>(model_ptr.get())) {
                sk->ensure_hoppings_cached();
            }
        }

        model.set_Dfield(cfg.Dfield_eV);

        // -----------------------------
        // 4) output roots:
        //    {data_dir}/band_{suffix}/D...
        // -----------------------------
        namespace fs = std::filesystem;

        const fs::path base_dir = cfg.data_dir.empty() ? fs::path("data") : fs::path(cfg.data_dir);

        // IMPORTANT: directories do NOT contain stamp; NO unit suffix
        const fs::path root_dir = base_dir / (std::string("band_") + cfg.output_suffix);
        const fs::path D_dir    = root_dir / (std::string("D") + tag_double(cfg.Dfield_eV, 3));

        if (rank == 0) {
            std::error_code ec;
            fs::create_directories(D_dir, ec);
            if (ec) {
                throw std::runtime_error("Failed to create output directory: " +
                                         D_dir.string() + " | " + ec.message());
            }

            std::cout << "=== test_bareband ===\n";
            std::cout << "config   = " << config_file << "\n";
            std::cout << "run_stamp= " << run_stamp << " (used only in filenames)\n";
            std::cout << "data_dir = " << base_dir.string() << "\n";
            std::cout << "suffix   = " << cfg.output_suffix << "\n";
            std::cout << "out_root = " << D_dir.string() << "\n";
            std::cout << "model    = " << cfg.model << "\n";
            std::cout << "layers   = " << para.layer_num << "\n";
            std::cout << "pressure = " << pressure << "\n";
            std::cout << "Dfield   = " << cfg.Dfield_eV << "\n";
#ifdef USE_MPI
            std::cout << "MPI      = " << nprocs
                      << " ranks (inited=" << (rgmpi::inited() ? "true" : "false") << ")\n";
#endif
            std::cout << "GMKG: Nk_seg=" << cfg.GMKG.Nk_seg << "\n";
            std::cout << "localK_MKKp: Nk_seg=" << cfg.localK_MKKp.Nk_seg
                      << " frac_local=" << cfg.localK_MKKp.frac_local << "\n";
        }
        rgmpi::barrier();

        // output filenames: only file name has stamp (same stamp for whole run)
        const fs::path out_gmkg  = D_dir / (std::string("GMKG_") + run_stamp + ".txt");
        const fs::path out_local = D_dir / (std::string("localK_MKKp_") + run_stamp + ".txt");

        // -----------------------------
        // (1) Global Γ–M–K–Γ
        // -----------------------------
        {
            const core::PathData path = st.generate_GMKG(cfg.GMKG.Nk_seg);

            // IMPORTANT: all ranks compute in same order
            const Eigen::MatrixXd E = model.bands_along_path(path);

            if (rank == 0) {
                write_header_only(out_gmkg.string(),
                                  cfg.model,
                                  config_file,
                                  run_stamp,
                                  cfg.output_suffix,
                                  st,
                                  model,
                                  path,
                                  /*path_name=*/"GMKG",
                                  /*Nk_seg=*/cfg.GMKG.Nk_seg,
                                  /*frac_local=*/0.0);
                rgio::write_bands_txt(out_gmkg.string(), path, E);
                std::cout << "Wrote: " << out_gmkg.string() << "\n";
            }
        }
        rgmpi::barrier();

        // -----------------------------
        // (2) Local M–K–K' near K
        // -----------------------------
        {
            const core::PathData path = st.generate_localK_MKKp(cfg.localK_MKKp.Nk_seg,
                                                                cfg.localK_MKKp.frac_local);

            // IMPORTANT: all ranks compute in same order
            const Eigen::MatrixXd E = model.bands_along_path(path);

            if (rank == 0) {
                write_header_only(out_local.string(),
                                  cfg.model,
                                  config_file,
                                  run_stamp,
                                  cfg.output_suffix,
                                  st,
                                  model,
                                  path,
                                  /*path_name=*/"localK_MKKp",
                                  /*Nk_seg=*/cfg.localK_MKKp.Nk_seg,
                                  /*frac_local=*/cfg.localK_MKKp.frac_local);
                rgio::write_bands_txt(out_local.string(), path, E);
                std::cout << "Wrote: " << out_local.string() << "\n";
            }
        }
        rgmpi::barrier();

        if (rank == 0) std::cout << "DONE.\n";

        rgmpi::finalize();
        return 0;

    } catch (const std::exception& e) {
        try {
            rgmpi::abort_all(std::string("test_bareband FAILED: ") + e.what(), 1);
        } catch (...) {}
        rgmpi::finalize();
        return 1;
    }
}
