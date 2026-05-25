#pragma once

#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>
#include "LinearAlgebra/MathFunctions.h"

namespace config {

using json = nlohmann::json;

// ============================================================
// Small utilities
// ============================================================

inline bool file_exists(const std::string& path)
{
    std::ifstream f(path);
    return static_cast<bool>(f);
}

inline std::vector<double> read_range(const json& j)
{
    const double xmin = j.at("min").get<double>();
    const double xmax = j.at("max").get<double>();
    const int num = j.at("num").get<int>();

    if (num < 1) {
        throw std::runtime_error("read_range: num must be >= 1");
    }

    std::vector<double> x(static_cast<size_t>(num));

    if (num == 1) {
        x[0] = xmin;
        return x;
    }

    for (int i = 0; i < num; ++i) {
        x[static_cast<size_t>(i)] =
            xmin
          + (xmax - xmin)
          * static_cast<double>(i)
          / static_cast<double>(num - 1);
    }

    return x;
}

// ============================================================
// Struct definitions
// ============================================================

struct RGPara {
    int    layer_num = 3;
    double a0        = 2.46;
    double d0        = 3.35;
    double vacuum    = 10.0;
    double pressure  = 0.0;
    double epsilon_r = 3.5;
};

struct MeshCfg {
    std::string type = "localK_hex";
    int    Nk        = 80;
    double dk_frac   = 0.1;
};

struct TestStructureConfig {
    std::string model     = "sk";
    std::string para_file = "Source/Parameters/RG_para.json";
    std::string data_dir  = "data";

    int Nk_seg = 200;
};

struct BarebandPathCfg {
    int Nk_seg = 200;
    double frac_local = 0.20;
};

struct TestBarebandConfig {
    std::string data_dir  = "data";
    std::string para_file = "Source/Parameters/RG_para.json";
    std::string output_suffix = "0";

    std::string model = "sk";
    double Dfield_eV = 0.0;

    BarebandPathCfg GMKG;
    BarebandPathCfg localK_MKKp;
};

struct DOSCfg {
    double e_low  = -0.3;
    double e_high =  0.3;
    int    num_e  = 2001;
    double eta    = 0.001;
};

struct CalDOSConfig {
    std::string data_dir  = "data";
    std::string para_file = "Source/Parameters/RG_para.json";

    std::string model = "sk";
    double T_K = 0.0;

    std::string output_suffix = "0";
    std::vector<double> Dfield_list_eV;

    MeshCfg kmesh;
    DOSCfg dos;
};

struct CalFermiConfig {
    std::string model = "sk";
    std::string data_dir  = "data";
    std::string output_suffix = "0";
    std::string para_file = "Source/Parameters/RG_para.json";

    double Dfield_eV = 0.0;
    double g_factor = 2.0;
    double orbital_derivative_dk = 1e-5;
    std::vector<double> Bfield_list_T;

    MeshCfg kmesh;   // local / truncated output patch

    std::vector<double> temperature_list;
    std::vector<double> mu_list;
};

struct CalOrbitalMomentConfig {
    std::string model = "sk";
    std::string data_dir  = "data";
    std::string output_suffix = "auto";
    std::string para_file = "Source/Parameters/RG_para.json";

    double Dfield_eV = 0.0;
    double derivative_dk = 1e-5;

    MeshCfg kmesh;
    std::string mode = "auto";
    bool use_kpath = false;
    std::string kpath_type = "localK_MKKp";
    int kpath_Nk_seg = 200;
    double kpath_frac_local = 0.20;
    std::vector<int> band_list;
};

using CalBerryCurvatureConfig = CalOrbitalMomentConfig;

inline CalOrbitalMomentConfig read_cal_orbital_moment_like_config(
    const std::string& config_file,
    const std::string& task_name
);

struct CalChiConfig {
    std::string data_dir = "data";
    std::string para_file = "Source/Parameters/RG_para.json";
    std::string output_suffix = "0";

    std::string model = "sk";
    std::string boundary = "open";

    std::string fermi_patch_path;

    double Dfield_eV = 0.0;
    double eta = 1e-5;

    int iq = 0;
    int jq = 0;

    bool use_form_factor = true;
    bool diagonal_band_only = false;

    std::vector<double> temperature_list;
    std::vector<double> polar_list; // stored in eV
    std::vector<double> mu_list;
};

struct CalChiV2Config {
    std::string data_dir = "data";
    std::string para_file = "Source/Parameters/RG_para.json";
    std::string output_suffix = "0";

    std::string model = "sk";
    std::string boundary = "open";

    double Dfield_eV = 0.0;

    int Nk = 400;
    int truncate_Nk = 40;
    int iq = 0;
    int jq = 0;

    double eta = 1e-6;
    bool use_form_factor = true;

    std::vector<double> temperature_list;
    std::vector<double> mu_list;
};

// ============================================================
// Common readers
// ============================================================

inline RGPara read_rg_para(const std::string& para_file)
{
    std::ifstream fin(para_file);

    if (!fin) {
        throw std::runtime_error(
            "Cannot open RG_para file: " + para_file
        );
    }

    json j;
    fin >> j;

    RGPara p;

    p.a0 =
        j.at("lattice").at("a0_A").get<double>();

    p.d0 =
        j.at("lattice").at("d0_A").get<double>();

    p.vacuum =
        j.at("lattice").at("vacuum_A").get<double>();

    p.pressure =
        j.at("lattice").at("pressure").get<double>();

    p.layer_num =
        j.at("lattice").at("layer_num").get<int>();

    p.epsilon_r =
        j.at("lattice").at("epsilon_r").get<double>();

    return p;
}

inline void read_io_common(
    const json& j,
    std::string& data_dir,
    std::string& para_file
) {
    if (!j.contains("io")) {
        return;
    }

    const auto& io = j.at("io");

    if (io.contains("data_dir")) {
        data_dir =
            io.at("data_dir").get<std::string>();
    }

    if (io.contains("model_para")) {
        para_file =
            io.at("model_para").get<std::string>();
    }
}

inline void read_mesh_cfg(
    const json& jmesh,
    MeshCfg& mesh
) {
    if (jmesh.contains("type")) {
        mesh.type =
            jmesh.at("type").get<std::string>();
    }

    if (jmesh.contains("Nk")) {
        mesh.Nk =
            jmesh.at("Nk").get<int>();
    }

    if (jmesh.contains("dk_frac")) {
        mesh.dk_frac =
            jmesh.at("dk_frac").get<double>();
    }
}

// ============================================================
// Task readers
// Top-level task layout, no "tasks" wrapper
// ============================================================

inline TestStructureConfig read_test_structure_config(
    const std::string& config_file
) {
    std::ifstream fin(config_file);

    if (!fin) {
        throw std::runtime_error(
            "Cannot open config file: " + config_file
        );
    }

    json j;
    fin >> j;

    TestStructureConfig c;

    read_io_common(
        j,
        c.data_dir,
        c.para_file
    );

    const auto& task =
        j.at("structure_and_kpath");

    if (task.contains("model")) {
        c.model =
            task.at("model").get<std::string>();
    }

    if (
        task.contains("GMKG")
     && task.at("GMKG").contains("Nk_seg")
    ) {
        c.Nk_seg =
            task.at("GMKG").at("Nk_seg").get<int>();
    }

    return c;
}

inline TestBarebandConfig read_test_bareband_config(
    const std::string& config_file
) {
    std::ifstream fin(config_file);

    if (!fin) {
        throw std::runtime_error(
            "Cannot open config file: " + config_file
        );
    }

    json j;
    fin >> j;

    TestBarebandConfig c;

    read_io_common(
        j,
        c.data_dir,
        c.para_file
    );

    const auto& task =
        j.at("bareband");

    if (task.contains("model")) {
        c.model =
            task.at("model").get<std::string>();
    }

    if (task.contains("Dfield_eV")) {
        c.Dfield_eV =
            task.at("Dfield_eV").get<double>();
    }

    if (task.contains("output_suffix")) {
        c.output_suffix =
            task.at("output_suffix").get<std::string>();
    }

    auto read_path =
        [](const json& jp, BarebandPathCfg& out) {
            if (jp.contains("Nk_seg")) {
                out.Nk_seg =
                    jp.at("Nk_seg").get<int>();
            }

            if (jp.contains("frac_local")) {
                out.frac_local =
                    jp.at("frac_local").get<double>();
            }
        };

    if (task.contains("GMKG")) {
        read_path(task.at("GMKG"), c.GMKG);
    }

    if (task.contains("localK_MKKp")) {
        read_path(task.at("localK_MKKp"), c.localK_MKKp);
    }

    return c;
}

inline CalDOSConfig read_cal_dos_config(
    const std::string& config_file
) {
    std::ifstream fin(config_file);

    if (!fin) {
        throw std::runtime_error(
            "Cannot open config file: " + config_file
        );
    }

    json j;
    fin >> j;

    CalDOSConfig c;

    read_io_common(
        j,
        c.data_dir,
        c.para_file
    );

    const auto& task =
        j.at("calculate_dos");

    if (task.contains("model")) {
        c.model =
            task.at("model").get<std::string>();
    }

    if (task.contains("T_K")) {
        c.T_K =
            task.at("T_K").get<double>();
    }

    if (task.contains("output_suffix")) {
        c.output_suffix =
            task.at("output_suffix").get<std::string>();
    }

    if (task.contains("kmesh")) {
        read_mesh_cfg(
            task.at("kmesh"),
            c.kmesh
        );
    }

    if (task.contains("dos")) {
        const auto& d =
            task.at("dos");

        if (d.contains("e_low")) {
            c.dos.e_low =
                d.at("e_low").get<double>();
        }

        if (d.contains("e_high")) {
            c.dos.e_high =
                d.at("e_high").get<double>();
        }

        if (d.contains("num_e")) {
            c.dos.num_e =
                d.at("num_e").get<int>();
        }

        if (d.contains("eta")) {
            c.dos.eta =
                d.at("eta").get<double>();
        }
    }

    const auto& Dlist =
        task.at("DfieldList_eV");

    const double Dmin =
        Dlist.at("min").get<double>();

    const double Dmax =
        Dlist.at("max").get<double>();

    const int nD =
        Dlist.at("num").get<int>();

    if (nD < 1) {
        throw std::runtime_error(
            "DfieldList_eV.num must be >= 1"
        );
    }

    if (Dmax < Dmin) {
        throw std::runtime_error(
            "DfieldList_eV.max must be >= min"
        );
    }

    c.Dfield_list_eV =
        la::linspace(Dmin, Dmax, nD);

    return c;
}

inline CalFermiConfig read_cal_fermi_config(
    const std::string& config_file
) {
    std::ifstream fin(config_file);

    if (!fin) {
        throw std::runtime_error(
            "Cannot open config file: " + config_file
        );
    }

    json j;
    fin >> j;

    CalFermiConfig c;

    read_io_common(
        j,
        c.data_dir,
        c.para_file
    );

    const auto& task =
        j.at("calculate_fermi");

    if (task.contains("model")) {
        c.model =
            task.at("model").get<std::string>();
    }

    if (task.contains("Dfield_eV")) {
        c.Dfield_eV =
            task.at("Dfield_eV").get<double>();
    }

    if (task.contains("BfieldList_T")) {
        c.Bfield_list_T =
            read_range(task.at("BfieldList_T"));
    }

    if (task.contains("g_factor")) {
        c.g_factor =
            task.at("g_factor").get<double>();
    }

    if (task.contains("orbital_derivative_dk")) {
        c.orbital_derivative_dk =
            task.at("orbital_derivative_dk").get<double>();
    }

    if (task.contains("magnetic_field")) {
        const auto& mf = task.at("magnetic_field");

        if (mf.contains("BList_T")) {
            c.Bfield_list_T =
                read_range(mf.at("BList_T"));
        }

        if (mf.contains("BfieldList_T")) {
            c.Bfield_list_T =
                read_range(mf.at("BfieldList_T"));
        }

        if (mf.contains("g_factor")) {
            c.g_factor =
                mf.at("g_factor").get<double>();
        }

        if (mf.contains("orbital_derivative_dk")) {
            c.orbital_derivative_dk =
                mf.at("orbital_derivative_dk").get<double>();
        }
    }

    if (task.contains("output_suffix")) {
        c.output_suffix =
            task.at("output_suffix").get<std::string>();
    }

    if (task.contains("kmesh")) {
        read_mesh_cfg(
            task.at("kmesh"),
            c.kmesh
        );
    }

    c.temperature_list =
        read_range(task.at("temperatureList_K"));

    c.mu_list =
        read_range(task.at("muList_eV"));

    if (c.Bfield_list_T.empty()) {
        c.Bfield_list_T.push_back(0.0);
    }

    return c;
}

inline CalOrbitalMomentConfig read_cal_orbital_moment_config(
    const std::string& config_file
) {
    return read_cal_orbital_moment_like_config(
        config_file,
        "calculate_orbital_moment"
    );
}

inline CalBerryCurvatureConfig read_cal_berry_curvature_config(
    const std::string& config_file
) {
    return read_cal_orbital_moment_like_config(
        config_file,
        "calculate_berry_curvature"
    );
}

inline CalOrbitalMomentConfig read_cal_orbital_moment_like_config(
    const std::string& config_file,
    const std::string& task_name
) {
    std::ifstream fin(config_file);

    if (!fin) {
        throw std::runtime_error(
            "Cannot open config file: " + config_file
        );
    }

    json j;
    fin >> j;

    CalOrbitalMomentConfig c;

    read_io_common(
        j,
        c.data_dir,
        c.para_file
    );

    const auto& task =
        j.at(task_name);

    if (task.contains("model")) {
        c.model =
            task.at("model").get<std::string>();
    }

    if (task.contains("Dfield_eV")) {
        c.Dfield_eV =
            task.at("Dfield_eV").get<double>();
    }

    if (task.contains("derivative_dk")) {
        c.derivative_dk =
            task.at("derivative_dk").get<double>();
    }

    if (task.contains("orbital_derivative_dk")) {
        c.derivative_dk =
            task.at("orbital_derivative_dk").get<double>();
    }

    if (task.contains("output_suffix")) {
        c.output_suffix =
            task.at("output_suffix").get<std::string>();
    }

    if (task.contains("mode")) {
        c.mode =
            task.at("mode").get<std::string>();
    }

    if (task.contains("k_mode")) {
        c.mode =
            task.at("k_mode").get<std::string>();
    }

    if (task.contains("kmesh")) {
        read_mesh_cfg(
            task.at("kmesh"),
            c.kmesh
        );
    }

    if (task.contains("kpath")) {
        const auto& kp =
            task.at("kpath");

        c.use_kpath = true;

        if (kp.contains("type")) {
            c.kpath_type =
                kp.at("type").get<std::string>();
        }

        if (kp.contains("Nk_seg")) {
            c.kpath_Nk_seg =
                kp.at("Nk_seg").get<int>();
        }

        if (kp.contains("frac_local")) {
            c.kpath_frac_local =
                kp.at("frac_local").get<double>();
        }
    }

    if (c.mode == "kpath" || c.mode == "path") {
        c.use_kpath = true;
    } else if (c.mode == "kmesh" || c.mode == "mesh") {
        c.use_kpath = false;
    } else if (c.mode == "auto") {
        // Legacy behavior: if kpath exists, use it.
    } else {
        throw std::runtime_error(
            task_name + ".mode must be kmesh, kpath, or auto"
        );
    }

    if (task.contains("bandList")) {
        c.band_list =
            task.at("bandList").get<std::vector<int>>();
    }

    if (task.contains("band_list")) {
        c.band_list =
            task.at("band_list").get<std::vector<int>>();
    }

    return c;
}

inline CalChiConfig read_cal_chi_config(
    const std::string& config_file
) {
    std::ifstream fin(config_file);

    if (!fin) {
        throw std::runtime_error(
            "Cannot open config file: " + config_file
        );
    }

    json j;
    fin >> j;

    CalChiConfig cfg;

    read_io_common(
        j,
        cfg.data_dir,
        cfg.para_file
    );

    const auto& c =
        j.at("calculate_chi");

    cfg.model =
        c.value("model", cfg.model);

    cfg.boundary =
        c.value("boundary", cfg.boundary);

    cfg.Dfield_eV =
        c.value("Dfield_eV", cfg.Dfield_eV);

    cfg.eta =
        c.value("eta", cfg.eta);

    cfg.iq =
        c.value("iq", cfg.iq);

    cfg.jq =
        c.value("jq", cfg.jq);

    cfg.fermi_patch_path =
        c.value("fermi_patch_path", cfg.fermi_patch_path);

    cfg.output_suffix =
        c.value("output_suffix", cfg.output_suffix);

    cfg.use_form_factor =
        c.value("use_form_factor", cfg.use_form_factor);

    cfg.diagonal_band_only =
        c.value("diagonal_band_only", cfg.diagonal_band_only);

    cfg.diagonal_band_only =
        c.value("diagonal_overlap_only", cfg.diagonal_band_only);

    cfg.polar_list =
        read_range(c.at("polarList_meV"));

    for (double& x : cfg.polar_list) {
        x *= 1e-3; // meV -> eV
    }

    cfg.temperature_list =
        read_range(c.at("temperatureList_K"));

    cfg.mu_list =
        read_range(c.at("muList_eV"));

    return cfg;
}

inline CalChiV2Config read_cal_chi_v2_config(
    const std::string& config_file
) {
    std::ifstream fin(config_file);

    if (!fin) {
        throw std::runtime_error(
            "Cannot open config file: " + config_file
        );
    }

    json j;
    fin >> j;

    CalChiV2Config cfg;

    read_io_common(
        j,
        cfg.data_dir,
        cfg.para_file
    );

    const auto& c =
        j.at("calculate_chi_v2");

    cfg.model =
        c.value("model", cfg.model);

    cfg.boundary =
        c.value("boundary", cfg.boundary);

    cfg.Dfield_eV =
        c.value("Dfield_eV", cfg.Dfield_eV);

    cfg.Nk =
        c.value("Nk", cfg.Nk);

    cfg.truncate_Nk =
        c.value("truncate_Nk", cfg.truncate_Nk);

    cfg.iq =
        c.value("iq", cfg.iq);

    cfg.jq =
        c.value("jq", cfg.jq);

    cfg.eta =
        c.value("eta", cfg.eta);

    cfg.use_form_factor =
        c.value("use_form_factor", cfg.use_form_factor);

    cfg.output_suffix =
        c.value("output_suffix", cfg.output_suffix);

    cfg.temperature_list =
        read_range(c.at("temperatureList_K"));

    cfg.mu_list =
        read_range(c.at("muList_eV"));

    return cfg;
}

} // namespace config
