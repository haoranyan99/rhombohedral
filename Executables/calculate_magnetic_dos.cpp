// File: Executables/calculate_magnetic_dos.cpp
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

#include <Eigen/Eigenvalues>

#include "LinearAlgebra/Constants.h"
#include "LinearAlgebra/MathFunctions.h"
#include "PhysStruct/RG_Structure.h"
#include "Source/Models/RG_KPModel.h"
#include "Source/Models/RG_ModelBase.h"
#include "Source/Models/RG_SKModel.h"
#include "Util/Config/config_reader.h"
#include "Util/MPI/mpi.h"
#include "Util/IO/rg_io.h"

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
        "calculate_magnetic_dos: unsupported kmesh.type = "
        + mesh.type
    );
    return core::GridData{};
}

static std::string auto_suffix(const config::CalMagneticDOSConfig& cfg)
{
    std::ostringstream oss;
    oss << cfg.kmesh.Nk
        << "_" << cfg.kmesh.type
        << "_" << compact_double_tag(cfg.kmesh.dk_frac);
    return oss.str();
}

static std::string output_suffix(const config::CalMagneticDOSConfig& cfg)
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
            config::read_cal_magnetic_dos_config(config_file);

        const auto para =
            config::read_rg_para(cfg.para_file);

        rg::RG_Structure st(
            para.layer_num,
            para.pressure,
            para.a0,
            para.d0,
            para.vacuum
        );

        core::GridData kmesh =
            make_kmesh(st, cfg.kmesh);
        kmesh.assert_consistent();
        const auto& kvec =
            kmesh.get<rg::RG_ModelBase::Vec2>("kvec").v;

        std::unique_ptr<rg::RG_ModelBase> model_ptr;
        if (cfg.model == "sk") {
            model_ptr =
                std::make_unique<rg::RG_SKModel>(st, cfg.para_file);
        } else if (cfg.model == "kp") {
            model_ptr =
                std::make_unique<rg::RG_KPModel>(st);
        } else {
            rgmpi::abort_all("calculate_magnetic_dos: model must be sk or kp");
        }

        rg::RG_ModelBase& model = *model_ptr;
        if (cfg.model == "sk") {
            if (auto* sk = dynamic_cast<rg::RG_SKModel*>(model_ptr.get())) {
                sk->ensure_hoppings_cached();
            }
        }

        const int num_e = cfg.dos.num_e;
        if (num_e < 2) {
            rgmpi::abort_all("calculate_magnetic_dos: dos.num_e must be >= 2");
        }

        const double e_low = cfg.dos.e_low;
        const double e_high = cfg.dos.e_high;
        if (!(e_high > e_low)) {
            rgmpi::abort_all("calculate_magnetic_dos: dos.e_high must be > e_low");
        }

        const double dE =
            (e_high - e_low) / static_cast<double>(num_e - 1);

        double eta = cfg.dos.eta;
        if (!(eta > 0.0)) eta = 0.6 * dE;
        if (!(eta > 0.0)) {
            rgmpi::abort_all("calculate_magnetic_dos: eta must be > 0");
        }

        std::vector<double> E_axis(static_cast<size_t>(num_e), 0.0);
        for (int ie = 0; ie < num_e; ++ie) {
            E_axis[static_cast<size_t>(ie)] =
                e_low + dE * static_cast<double>(ie);
        }

        const double inv_sqrt2pi =
            1.0 / std::sqrt(2.0 * la::pi);
        const double pref =
            inv_sqrt2pi / eta;
        const double inv2eta2 =
            1.0 / (2.0 * eta * eta);

        const auto [i0, i1] =
            rgmpi::block_1d_int(static_cast<int>(kmesh.size()), rank, nprocs);

        const fs::path base_dir =
            cfg.data_dir.empty() ? fs::path("data") : fs::path(cfg.data_dir);

        const fs::path out_root =
            base_dir
            / (std::string("magnetic_dos_") + cfg.model + "_" + output_suffix(cfg));

        const std::string run_stamp =
            rgio::make_time_stamp("magnetic_dos");

        for (double Dfield_eV : cfg.Dfield_list_eV) {
            model.set_Dfield(Dfield_eV);

            const int dim =
                static_cast<int>(model.bands_at_k(kvec.front(), true).size());

            for (double B_T : cfg.Bfield_list_T) {
                std::vector<double> dos_up_local(static_cast<size_t>(num_e), 0.0);
                std::vector<double> dos_dn_local(static_cast<size_t>(num_e), 0.0);
                std::vector<double> fill_up_local(static_cast<size_t>(num_e), 0.0);
                std::vector<double> fill_dn_local(static_cast<size_t>(num_e), 0.0);

                for (int ii = i0; ii < i1; ++ii) {
                    const size_t ik = static_cast<size_t>(ii);

                    Eigen::VectorXd ev =
                        model.bands_at_k(kvec[ik], true);

                    const std::vector<double> m =
                        model.orbital_moment_muB_at_k(
                            kvec[ik],
                            cfg.orbital_derivative_dk,
                            true
                        );

                    for (int b = 0; b < dim; ++b) {
                        const double e_up =
                            ev(b)
                            + magnetic_shift(
                                m[static_cast<size_t>(b)],
                                B_T,
                                cfg.g_factor,
                                +1
                            );

                        const double e_dn =
                            ev(b)
                            + magnetic_shift(
                                m[static_cast<size_t>(b)],
                                B_T,
                                cfg.g_factor,
                                -1
                            );

                        for (int ie = 0; ie < num_e; ++ie) {
                            const size_t i = static_cast<size_t>(ie);
                            const double E = E_axis[i];
                            const double de_up = E - e_up;
                            const double de_dn = E - e_dn;

                            dos_up_local[i] +=
                                std::exp(-(de_up * de_up) * inv2eta2) * pref;

                            dos_dn_local[i] +=
                                std::exp(-(de_dn * de_dn) * inv2eta2) * pref;

                            fill_up_local[i] +=
                                la::fermi(e_up, E, cfg.T_K);

                            fill_dn_local[i] +=
                                la::fermi(e_dn, E, cfg.T_K);
                        }
                    }
                }

                rgmpi::allreduce_sum_vector(dos_up_local);
                rgmpi::allreduce_sum_vector(dos_dn_local);
                rgmpi::allreduce_sum_vector(fill_up_local);
                rgmpi::allreduce_sum_vector(fill_dn_local);

                if (rank == 0) {
                    const fs::path out_dir =
                        out_root
                        / (std::string("D") + rgio::tag3(Dfield_eV))
                        / B_tag(B_T);

                    fs::create_directories(out_dir);

                    const fs::path out_path =
                        out_dir / (run_stamp + ".txt");

                    std::ofstream ofs(out_path.string());
                    if (!ofs) {
                        throw std::runtime_error(
                            "Cannot open output file: " + out_path.string()
                        );
                    }

                    const double invNk =
                        1.0 / static_cast<double>(kmesh.size());

                    const double invNkDim =
                        1.0 / (
                            static_cast<double>(kmesh.size())
                          * static_cast<double>(dim)
                        );

                    ofs << std::setprecision(15);
                    ofs << "# --- calculate_magnetic_dos ---\n";
                    ofs << "# config_file = " << config_file << "\n";
                    ofs << "# model = " << cfg.model << "\n";
                    ofs << "# Dfield_eV = " << Dfield_eV << "\n";
                    ofs << "# B_T = " << B_T << "\n";
                    ofs << "# g_factor = " << cfg.g_factor << "\n";
                    ofs << "# orbital_derivative_dk = "
                        << cfg.orbital_derivative_dk << "\n";
                    ofs << "# T_K = " << cfg.T_K << "\n";
                    ofs << "# kmesh = " << cfg.kmesh.type
                        << " Nk=" << cfg.kmesh.Nk
                        << " dk_frac=" << cfg.kmesh.dk_frac
                        << " mesh_type=" << kmesh.mesh_type
                        << " size=" << kmesh.size() << "\n";
                    ofs << "# E_scan = [" << e_low << ", " << e_high
                        << "], num_e=" << num_e
                        << ", eta=" << eta << "\n";
                    ofs << "# dos_up/dos_down are per spin; dos_sum is spin-summed; dos_avg is spin-averaged.\n";
                    ofs << "# columns = i E_eV filling_up filling_down filling_avg doping_avg_1e12cm^-2 dos_up dos_down dos_sum dos_avg\n";
                    ofs << "# dE = " << dE << "\n";

                    for (int ie = 0; ie < num_e; ++ie) {
                        const size_t i = static_cast<size_t>(ie);
                        const double dos_up = dos_up_local[i] * invNk;
                        const double dos_dn = dos_dn_local[i] * invNk;
                        const double filling_up =
                            fill_up_local[i] * invNkDim;
                        const double filling_dn =
                            fill_dn_local[i] * invNkDim;
                        const double filling_avg =
                            0.5 * (filling_up + filling_dn);
                        const double doping_avg =
                            la::filling_to_doping(filling_avg, st.a());

                        ofs << ie << " "
                            << E_axis[i] << " "
                            << filling_up << " "
                            << filling_dn << " "
                            << filling_avg << " "
                            << doping_avg << " "
                            << dos_up << " "
                            << dos_dn << " "
                            << dos_up + dos_dn << " "
                            << 0.5 * (dos_up + dos_dn)
                            << "\n";
                    }

                    std::cout << "Wrote: " << out_path.string() << "\n";
                }
            }
        }

        rgmpi::finalize();
        return 0;

    } catch (const std::exception& e) {
        try {
            rgmpi::abort_all(
                std::string("calculate_magnetic_dos FAILED: ") + e.what(),
                1
            );
        } catch (...) {}
        rgmpi::finalize();
        return 1;
    }
}
