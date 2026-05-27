// File: Executables/calculate_berry_curvature.cpp
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef USE_MPI
  #include <mpi.h>
#endif

#include "Common/DataContainers.h"
#include "LinearAlgebra/MathFunctions.h"
#include "PhysStruct/RG_Structure.h"
#include "Source/Models/RG_KPModel.h"
#include "Source/Models/RG_ModelBase.h"
#include "Source/Models/RG_SKModel.h"
#include "Util/Config/config_reader.h"
#include "Util/IO/rg_io.h"
#include "Util/MPI/mpi.h"

namespace fs = std::filesystem;

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

    rgmpi::abort_all(
        "calculate_berry_curvature: unsupported kmesh.type = "
        + mesh.type
    );
    return core::GridData{};
}

static core::PathData make_kpath(
    const rg::RG_Structure& st,
    const config::CalBerryCurvatureConfig& cfg
) {
    if (cfg.kpath_type == "GMKG") {
        return st.generate_GMKG(cfg.kpath_Nk_seg);
    }

    if (
        cfg.kpath_type == "localK_MKKp"
     || cfg.kpath_type == "MKKp"
    ) {
        return st.generate_localK_MKKp(
            cfg.kpath_Nk_seg,
            cfg.kpath_frac_local
        );
    }

    rgmpi::abort_all(
        "calculate_berry_curvature: unsupported kpath.type = "
        + cfg.kpath_type
    );
    return core::PathData{};
}

static std::string auto_output_suffix(
    const config::CalBerryCurvatureConfig& cfg
) {
    std::ostringstream oss;
    if (cfg.use_kpath) {
        oss << cfg.kpath_type
            << "_Nkseg"
            << cfg.kpath_Nk_seg;

        if (
            cfg.kpath_type == "localK_MKKp"
         || cfg.kpath_type == "MKKp"
        ) {
            oss << "_frac"
                << compact_double_tag(cfg.kpath_frac_local);
        }
    } else {
        oss << cfg.kmesh.Nk
            << "_"
            << cfg.kmesh.type
            << "_"
            << compact_double_tag(cfg.kmesh.dk_frac);
    }

    return oss.str();
}

static std::string output_suffix(
    const config::CalBerryCurvatureConfig& cfg
) {
    std::string suffix =
        auto_output_suffix(cfg);

    if (
        !cfg.output_suffix.empty()
     && cfg.output_suffix != "auto"
    ) {
        suffix += "_" + cfg.output_suffix;
    }

    return suffix;
}

static std::vector<int> default_band_list(int dim)
{
    std::vector<int> bands(static_cast<size_t>(dim));
    for (int i = 0; i < dim; ++i) {
        bands[static_cast<size_t>(i)] = i;
    }
    return bands;
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

        const auto cfg =
            config::read_cal_berry_curvature_config(config_file);

        const auto para =
            config::read_rg_para(cfg.para_file);

        rg::RG_Structure st(
            para.layer_num,
            para.pressure,
            para.a0,
            para.d0,
            para.vacuum
        );

        core::GridData kmesh;
        core::PathData kpath;
        std::vector<rg::RG_ModelBase::Vec2> kvec_storage;
        const std::vector<rg::RG_ModelBase::Vec2>* kvec_ptr = nullptr;

        if (cfg.use_kpath) {
            kpath =
                make_kpath(st, cfg);
            kpath.assert_consistent();
            kvec_ptr =
                &kpath.k_list;
        } else {
            kmesh =
                make_kmesh(st, cfg.kmesh);
            kmesh.assert_consistent();
            kvec_storage =
                kmesh.get<rg::RG_ModelBase::Vec2>("kvec").v;
            kvec_ptr =
                &kvec_storage;
        }

        const auto& kvec =
            *kvec_ptr;

        std::unique_ptr<rg::RG_ModelBase> model_ptr;

        if (cfg.model == "sk") {
            model_ptr =
                std::make_unique<rg::RG_SKModel>(
                    st,
                    cfg.para_file
                );
        } else if (cfg.model == "kp") {
            model_ptr =
                std::make_unique<rg::RG_KPModel>(st);
        } else {
            rgmpi::abort_all(
                "calculate_berry_curvature: model must be sk or kp"
            );
        }

        rg::RG_ModelBase& model =
            *model_ptr;

        if (cfg.model == "sk") {
            if (auto* sk =
                    dynamic_cast<rg::RG_SKModel*>(model_ptr.get()))
            {
                sk->ensure_hoppings_cached();
            }
        }

        model.set_Dfield(cfg.Dfield_eV);

        Eigen::VectorXd e0 =
            model.bands_at_k(kvec.front(), true);

        const int dim =
            static_cast<int>(e0.size());

        std::vector<int> bands =
            cfg.band_list.empty()
            ? default_band_list(dim)
            : cfg.band_list;

        for (int b : bands) {
            if (b < 0 || b >= dim) {
                rgmpi::abort_all(
                    "calculate_berry_curvature: band index out of range"
                );
            }
        }

        const size_t NkTot =
            kvec.size();

        const auto [i0, i1] =
            rgmpi::block_1d_int(
                static_cast<int>(NkTot),
                rank,
                nprocs
            );

        const size_t nb =
            bands.size();

        std::vector<double> energy_all(
            NkTot * nb,
            0.0
        );

        std::vector<double> berry_all(
            NkTot * nb,
            0.0
        );

        for (int ii = i0; ii < i1; ++ii) {
            const size_t ik =
                static_cast<size_t>(ii);

            Eigen::VectorXd ev =
                model.bands_at_k(kvec[ik], true);

            const std::vector<double> berry =
                model.berry_curvature_A2_at_k(
                    kvec[ik],
                    cfg.derivative_dk,
                    true
                );

            for (size_t ib = 0; ib < nb; ++ib) {
                const int b =
                    bands[ib];

                const size_t pos =
                    ik * nb + ib;

                energy_all[pos] =
                    ev(b);

                berry_all[pos] =
                    berry[static_cast<size_t>(b)];
            }
        }

        rgmpi::allreduce_sum_vector(energy_all);
        rgmpi::allreduce_sum_vector(berry_all);

        const fs::path base_dir =
            cfg.data_dir.empty()
            ? fs::path("data")
            : fs::path(cfg.data_dir);

        const std::string suffix =
            output_suffix(cfg);

        const std::string run_stamp =
            rgio::make_time_stamp("berry_curvature");

        if (rank == 0) {
            const fs::path out_root =
                base_dir
                / (
                    std::string("berry_curvature_")
                  + cfg.model
                  + "_"
                  + suffix
                )
                / (std::string("D") + rgio::tag3(cfg.Dfield_eV));

            fs::create_directories(out_root);

            const fs::path out_path =
                out_root / (run_stamp + ".txt");

            std::ofstream ofs(out_path.string());
            if (!ofs) {
                throw std::runtime_error(
                    "Cannot open output file: "
                    + out_path.string()
                );
            }

            ofs << std::setprecision(15);
            ofs << "# --- calculate_berry_curvature ---\n";
            ofs << "# config_file = " << config_file << "\n";
            ofs << "# model = " << cfg.model << "\n";
            ofs << "# Dfield_eV = " << cfg.Dfield_eV << "\n";
            ofs << "# derivative_dk = " << cfg.derivative_dk << "\n";
            if (cfg.use_kpath) {
                ofs << "# mode = kpath\n";
                ofs << "# kpath = " << cfg.kpath_type
                    << " Nk_seg=" << cfg.kpath_Nk_seg
                    << " frac_local=" << cfg.kpath_frac_local
                    << " size=" << NkTot << "\n";
                for (size_t it = 0; it < kpath.xtick_pos.size(); ++it) {
                    ofs << "# xtick = "
                        << kpath.xtick_pos[it]
                        << " "
                        << kpath.xtick_lab[it]
                        << "\n";
                }
            } else {
                ofs << "# mode = kmesh\n";
                ofs << "# kmesh = " << cfg.kmesh.type
                    << " Nk=" << cfg.kmesh.Nk
                    << " dk_frac=" << cfg.kmesh.dk_frac
                    << " mesh_type=" << kmesh.mesh_type
                    << " size=" << NkTot << "\n";
            }
            ofs << "# band_index_base = 0\n";
            ofs << "# berry_curvature_unit = Angstrom^2\n";
            ofs << "# energy_unit = eV\n";
            ofs << "# k_unit = Angstrom^-1\n";
            if (cfg.use_kpath) {
                ofs << "# kline_unit = Angstrom^-1\n";
                ofs << "# columns = ik s kx ky";
            } else {
                ofs << "# columns = ik iq jq kx ky";
            }
            for (int b : bands) {
                ofs << " E_b" << b
                    << " berry_A2_b" << b;
            }
            ofs << "\n";

            for (size_t ik = 0; ik < NkTot; ++ik) {
                if (cfg.use_kpath) {
                    ofs << ik << " "
                        << kpath.kline[ik] << " "
                        << kvec[ik].x() << " "
                        << kvec[ik].y();
                } else {
                    ofs << ik << " "
                        << kmesh.iq[ik] << " "
                        << kmesh.jq[ik] << " "
                        << kvec[ik].x() << " "
                        << kvec[ik].y();
                }

                for (size_t ib = 0; ib < nb; ++ib) {
                    const size_t pos =
                        ik * nb + ib;

                    ofs << " "
                        << energy_all[pos]
                        << " "
                        << berry_all[pos];
                }
                ofs << "\n";
            }

            std::cout << "=== calculate_berry_curvature ===\n";
            std::cout << "config    = " << config_file << "\n";
            std::cout << "out_path  = " << out_path.string() << "\n";
            std::cout << "model     = " << cfg.model << "\n";
            if (cfg.use_kpath) {
                std::cout << "kpath     = " << cfg.kpath_type
                          << " Nk_seg=" << cfg.kpath_Nk_seg
                          << " frac_local=" << cfg.kpath_frac_local
                          << " size=" << NkTot << "\n";
            } else {
                std::cout << "kmesh     = " << cfg.kmesh.type
                          << " Nk=" << cfg.kmesh.Nk
                          << " dk_frac=" << cfg.kmesh.dk_frac
                          << " size=" << NkTot << "\n";
            }
            std::cout << "bands     =";
            for (int b : bands) {
                std::cout << " " << b;
            }
            std::cout << "\nDONE.\n";
        }

        rgmpi::finalize();
        return 0;

    } catch (const std::exception& e) {
        try {
            rgmpi::abort_all(
                std::string("calculate_berry_curvature FAILED: ")
                + e.what(),
                1
            );
        } catch (...) {}

        rgmpi::finalize();
        return 1;
    }
}
