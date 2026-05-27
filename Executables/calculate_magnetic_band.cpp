// File: Executables/calculate_magnetic_band.cpp
#include <cmath>
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

#include "LinearAlgebra/Constants.h"
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
    if (std::abs(x) < 5e-13) x = 0.0;
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(10) << x;
    std::string s = oss.str();
    while (!s.empty() && s.back() == '0') s.pop_back();
    if (!s.empty() && s.back() == '.') s.pop_back();
    if (s == "-0") s = "0";
    return s;
}

static std::string B_tag(double B_T)
{
    return std::string("B") + compact_double_tag(B_T) + "T";
}

static double magnetic_shift(
    double m_orb_muB,
    double B_T,
    double g_factor,
    int spin_sign
) {
    const int s = (spin_sign >= 0) ? +1 : -1;
    return -m_orb_muB * la::muB_eV_per_T * B_T
         + static_cast<double>(s)
              * g_factor * la::muB_eV_per_T * B_T;
}

static core::PathData make_kpath(
    const rg::RG_Structure& st,
    const config::CalMagneticBandConfig& cfg
) {
    if (cfg.kpath_type == "GMKG") {
        return st.generate_GMKG(cfg.kpath_Nk_seg);
    }

    if (cfg.kpath_type == "localK_MKKp" || cfg.kpath_type == "MKKp") {
        return st.generate_localK_MKKp(
            cfg.kpath_Nk_seg,
            cfg.kpath_frac_local
        );
    }

    rgmpi::abort_all(
        "calculate_magnetic_band: unsupported kpath.type = "
        + cfg.kpath_type
    );
    return core::PathData{};
}

static std::vector<int> default_band_list(int dim)
{
    std::vector<int> out(static_cast<size_t>(dim));
    for (int i = 0; i < dim; ++i) out[static_cast<size_t>(i)] = i;
    return out;
}

static std::string auto_suffix(const config::CalMagneticBandConfig& cfg)
{
    std::ostringstream oss;
    oss << cfg.kpath_type
        << "_Nkseg" << cfg.kpath_Nk_seg;

    if (cfg.kpath_type == "localK_MKKp" || cfg.kpath_type == "MKKp") {
        oss << "_frac" << compact_double_tag(cfg.kpath_frac_local);
    }

    oss << "_D" << compact_double_tag(cfg.Dfield_eV);
    return oss.str();
}

static std::string output_suffix(const config::CalMagneticBandConfig& cfg)
{
    std::string suffix = auto_suffix(cfg);
    if (!cfg.output_suffix.empty() && cfg.output_suffix != "auto") {
        suffix += "_" + cfg.output_suffix;
    }
    return suffix;
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
            config::read_cal_magnetic_band_config(config_file);

        const auto para =
            config::read_rg_para(cfg.para_file);

        rg::RG_Structure st(
            para.layer_num,
            para.pressure,
            para.a0,
            para.d0,
            para.vacuum
        );

        core::PathData kpath = make_kpath(st, cfg);
        kpath.assert_consistent();
        const auto& kvec = kpath.k_list;

        std::unique_ptr<rg::RG_ModelBase> model_ptr;
        if (cfg.model == "sk") {
            model_ptr =
                std::make_unique<rg::RG_SKModel>(st, cfg.para_file);
        } else if (cfg.model == "kp") {
            model_ptr =
                std::make_unique<rg::RG_KPModel>(st);
        } else {
            rgmpi::abort_all("calculate_magnetic_band: model must be sk or kp");
        }

        rg::RG_ModelBase& model = *model_ptr;
        if (cfg.model == "sk") {
            if (auto* sk = dynamic_cast<rg::RG_SKModel*>(model_ptr.get())) {
                sk->ensure_hoppings_cached();
            }
        }
        model.set_Dfield(cfg.Dfield_eV);

        const int dim =
            static_cast<int>(model.bands_at_k(kvec.front(), true).size());

        const std::vector<int> bands =
            cfg.band_list.empty() ? default_band_list(dim) : cfg.band_list;

        for (int b : bands) {
            if (b < 0 || b >= dim) {
                rgmpi::abort_all("calculate_magnetic_band: band index out of range");
            }
        }

        const size_t NkTot = kvec.size();
        const size_t nb = bands.size();

        std::vector<double> E0(NkTot * nb, 0.0);
        std::vector<double> morb(NkTot * nb, 0.0);

        const auto [i0, i1] =
            rgmpi::block_1d_int(static_cast<int>(NkTot), rank, nprocs);

        for (int ii = i0; ii < i1; ++ii) {
            const size_t ik = static_cast<size_t>(ii);
            Eigen::VectorXd ev = model.bands_at_k(kvec[ik], true);
            const std::vector<double> m =
                model.orbital_moment_muB_at_k(
                    kvec[ik],
                    cfg.derivative_dk,
                    true
                );

            for (size_t ib = 0; ib < nb; ++ib) {
                const int b = bands[ib];
                const size_t pos = ik * nb + ib;
                E0[pos] = ev(b);
                morb[pos] = m[static_cast<size_t>(b)];
            }
        }

        rgmpi::allreduce_sum_vector(E0);
        rgmpi::allreduce_sum_vector(morb);

        if (rank == 0) {
            const fs::path base_dir =
                cfg.data_dir.empty() ? fs::path("data") : fs::path(cfg.data_dir);

            const fs::path out_root =
                base_dir
                / (std::string("magnetic_band_") + cfg.model + "_" + output_suffix(cfg))
                / (std::string("D") + rgio::tag3(cfg.Dfield_eV));

            fs::create_directories(out_root);

            const std::string run_stamp =
                rgio::make_time_stamp("magnetic_band");

            for (double B_T : cfg.Bfield_list_T) {
                const fs::path B_dir =
                    out_root / B_tag(B_T);
                fs::create_directories(B_dir);

                const fs::path out_path =
                    B_dir / (run_stamp + ".txt");

                std::ofstream ofs(out_path.string());
                if (!ofs) {
                    throw std::runtime_error(
                        "Cannot open output file: " + out_path.string()
                    );
                }

                ofs << std::setprecision(15);
                ofs << "# --- calculate_magnetic_band ---\n";
                ofs << "# config_file = " << config_file << "\n";
                ofs << "# model = " << cfg.model << "\n";
                ofs << "# Dfield_eV = " << cfg.Dfield_eV << "\n";
                ofs << "# B_T = " << B_T << "\n";
                ofs << "# g_factor = " << cfg.g_factor << "\n";
                ofs << "# orbital_derivative_dk = " << cfg.derivative_dk << "\n";
                ofs << "# energy_unit = eV\n";
                ofs << "# m_orb_unit = mu_B\n";
                ofs << "# k_unit = Angstrom^-1\n";
                ofs << "# kline_unit = Angstrom^-1\n";
                ofs << "# kpath = " << cfg.kpath_type
                    << " Nk_seg=" << cfg.kpath_Nk_seg
                    << " frac_local=" << cfg.kpath_frac_local
                    << " size=" << NkTot << "\n";
                for (size_t it = 0; it < kpath.xtick_pos.size(); ++it) {
                    ofs << "# xtick = " << kpath.xtick_pos[it]
                        << " " << kpath.xtick_lab[it] << "\n";
                }
                ofs << "# columns = ik s kx ky";
                for (int b : bands) {
                    ofs << " E0_b" << b
                        << " m_orb_muB_b" << b
                        << " Eup_b" << b
                        << " Edn_b" << b;
                }
                ofs << "\n";

                for (size_t ik = 0; ik < NkTot; ++ik) {
                    ofs << ik << " "
                        << kpath.kline[ik] << " "
                        << kvec[ik].x() << " "
                        << kvec[ik].y();

                    for (size_t ib = 0; ib < nb; ++ib) {
                        const size_t pos = ik * nb + ib;
                        const double e0 = E0[pos];
                        const double m = morb[pos];
                        ofs << " " << e0
                            << " " << m
                            << " " << e0 + magnetic_shift(m, B_T, cfg.g_factor, +1)
                            << " " << e0 + magnetic_shift(m, B_T, cfg.g_factor, -1);
                    }
                    ofs << "\n";
                }

                std::cout << "Wrote: " << out_path.string() << "\n";
            }
        }

        rgmpi::finalize();
        return 0;

    } catch (const std::exception& e) {
        try {
            rgmpi::abort_all(
                std::string("calculate_magnetic_band FAILED: ") + e.what(),
                1
            );
        } catch (...) {}
        rgmpi::finalize();
        return 1;
    }
}
