// File: Util/Config/config_reader.h
#pragma once

#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <sstream>

#include <nlohmann/json.hpp>
#include "LinearAlgebra/MathFunctions.h" // la::linspace, la::doping_to_filling

namespace config {

using json = nlohmann::json;

#pragma region ===================== Small utilities =====================

inline bool file_exists(const std::string& path) {
    std::ifstream f(path);
    return (bool)f;
}

#pragma endregion
// ============================================================================



#pragma region ===================== Struct definitions =====================

// -------------------- RG parameters --------------------
struct RGPara {
    int    layer_num = 3;
    double a0        = 2.46;
    double d0        = 3.35;
    double vacuum    = 10.0;
    double pressure  = 0.0;
};


// -------------------- mesh --------------------
struct MeshCfg {
    std::string type = "localK_hex";
    int    Nk        = 80;
    double dk_frac   = 0.1;
};

// -------------------- structure / band --------------------
struct TestStructureConfig {
    std::string model     = "sk";
    std::string para_file = "RG_para.json";
    std::string data_dir  = "data";
    int Nk_seg            = 200;
};

struct BarebandPathCfg {
    int Nk_seg = 200;
    double frac_local = 0.20; // for local path; GMKG can ignore but keep consistent
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

// -------------------- DOS --------------------
struct DOSCfg {
    double e_low  = -0.3;
    double e_high =  0.3;
    int    num_e  = 2001;
    double eta    = 0.001; // positive by default
};

struct TestDOSConfig {
    std::string data_dir  = "data";
    std::string para_file = "Source/Parameters/RG_para.json";
    std::string model = "sk";

    double Dfield_eV  = 0.0;
    double T_K = 0.0;

    MeshCfg kmesh;
    DOSCfg dos;
};


struct CalFermiConfig {
    std::string model = "sk";
    std::string mode = "doping";
    std::string data_dir  = "data";
    std::string output_suffix = "0";
    std::string para_file = "Source/Parameters/RG_para.json";
    

    double Dfield_eV  = 0.0;

    MeshCfg kmesh;
    std::vector<double> temperature_list;
    std::vector<double> doping_list;
    std::vector<double> filling_list;
    std::vector<double> mu_list;
};

// -------------------- batch chi --------------------
struct CalChiConfig {
    std::string data_dir  = "data";
    std::string para_file = "Source/Parameters/RG_para.json";
    std::string output_suffix = "0";
    std::string model = "sk";
    std::string mode = "doping";
    std::string boundary = "periodic";
    std::string fermi_patch_path; 

    double Dfield_eV  = 0.0;
    int iq_min = -10;
    int iq_max = 10;
    int jq_min = -10;
    int jq_max = 10;
    double eta = 0.0001;

    MeshCfg kmesh;

    std::vector<double> temperature_list;
    std::vector<double> polar_list;
    std::vector<double> doping_list;
    std::vector<double> filling_list;
    std::vector<double> mu_list;

};



#pragma endregion
// ============================================================================





#pragma region ===================== Readers (JSON) =====================

inline RGPara read_rg_para(const std::string& para_file)
{
    std::ifstream fin(para_file);
    if (!fin)
        throw std::runtime_error("Cannot open RG_para file: " + para_file);

    json j;
    fin >> j;

    RGPara p;

    // lattice
    p.a0        = j.at("lattice").at("a0_A").get<double>();
    p.d0        = j.at("lattice").at("d0_A").get<double>();
    p.vacuum    = j.at("lattice").at("vacuum_A").get<double>();
    p.pressure  = j.at("lattice").at("pressure").get<double>();
    p.layer_num = j.at("lattice").at("layer_num").get<int>();

    return p;
}

inline TestStructureConfig read_test_structure_config(const std::string& config_file)
{
    std::ifstream fin(config_file);
    if (!fin)
        throw std::runtime_error("Cannot open config file: " + config_file);

    json j;
    fin >> j;

    TestStructureConfig c;

    // io
    if (j.contains("io")) {
        const auto& io = j.at("io");
        if (io.contains("model_para"))
            c.para_file = io.at("model_para").get<std::string>();
        if (io.contains("data_dir"))
            c.data_dir = io.at("data_dir").get<std::string>();
    }

    // task-specific
    const auto& task = j.at("tasks").at("structure_and_kpath");

    if (task.contains("model"))
        c.model = task.at("model").get<std::string>();

    if (task.contains("GMKG") && task.at("GMKG").contains("Nk_seg"))
        c.Nk_seg = task.at("GMKG").at("Nk_seg").get<int>();

    return c;
}

inline TestBarebandConfig read_test_bareband_config(const std::string& config_file)
{
    std::ifstream fin(config_file);
    if (!fin) throw std::runtime_error("Cannot open config file: " + config_file);

    json j;
    fin >> j;

    TestBarebandConfig c;

    // io
    if (j.contains("io")) {
        const auto& io = j.at("io");
        if (io.contains("data_dir"))   c.data_dir  = io.at("data_dir").get<std::string>();
        if (io.contains("model_para")) c.para_file = io.at("model_para").get<std::string>();
    }

    // tasks.bareband
    const auto& task = j.at("tasks").at("bareband");

    if (task.contains("model"))     c.model     = task.at("model").get<std::string>();
    if (task.contains("output_suffix")) c.output_suffix = task.at("output_suffix").get<std::string>();
    if (task.contains("Dfield_eV")) c.Dfield_eV = task.at("Dfield_eV").get<double>();

    auto read_path = [&](const json& jp, BarebandPathCfg& out) {
        if (jp.contains("Nk_seg"))     out.Nk_seg = jp.at("Nk_seg").get<int>();
        if (jp.contains("frac_local")) out.frac_local = jp.at("frac_local").get<double>();
    };

    if (task.contains("GMKG"))        read_path(task.at("GMKG"), c.GMKG);
    if (task.contains("localK_MKKp")) read_path(task.at("localK_MKKp"), c.localK_MKKp);

    return c;
}

inline TestDOSConfig read_test_dos_config(const std::string& config_file)
{
    std::ifstream fin(config_file);
    if (!fin) throw std::runtime_error("Cannot open config file: " + config_file);

    json j;
    fin >> j;

    TestDOSConfig c;

    // io
    if (j.contains("io")) {
        const auto& io = j.at("io");
        if (io.contains("data_dir"))   c.data_dir  = io.at("data_dir").get<std::string>();
        if (io.contains("model_para")) c.para_file = io.at("model_para").get<std::string>();
    }

    // tasks.test_dos
    const auto& task = j.at("tasks").at("test_dos");

    if (task.contains("model"))     c.model     = task.at("model").get<std::string>();
    if (task.contains("Dfield_eV")) c.Dfield_eV = task.at("Dfield_eV").get<double>();
    if (task.contains("T_K"))       c.T_K       = task.at("T_K").get<double>();

    // kmesh
    if (task.contains("kmesh")) {
        const auto& km = task.at("kmesh");
        if (km.contains("type"))    c.kmesh.type = km.at("type").get<std::string>();
        if (km.contains("Nk"))      c.kmesh.Nk   = km.at("Nk").get<int>();
        if (km.contains("dk_frac")) c.kmesh.dk_frac = km.at("dk_frac").get<double>();
    }

    // dos
    if (task.contains("dos")) {
        const auto& d = task.at("dos");
        if (d.contains("e_low"))  c.dos.e_low  = d.at("e_low").get<double>();
        if (d.contains("e_high")) c.dos.e_high = d.at("e_high").get<double>();
        if (d.contains("num_e"))  c.dos.num_e  = d.at("num_e").get<int>();
        if (d.contains("eta"))    c.dos.eta    = d.at("eta").get<double>();
    }

    return c;
}

inline CalFermiConfig read_cal_fermi_config(const std::string& config_file)
{
    std::ifstream fin(config_file);
    if (!fin)
        throw std::runtime_error("Cannot open config file: " + config_file);

    json j;
    fin >> j;

    CalFermiConfig c;

    // ---------------- io ----------------
    if (j.contains("io")) {
        const auto& io = j.at("io");
        if (io.contains("data_dir"))   c.data_dir  = io.at("data_dir").get<std::string>();
        if (io.contains("model_para")) c.para_file = io.at("model_para").get<std::string>();
    }

    // ---------------- task ----------------
    const auto& task = j.at("tasks").at("calculate_fermi");

    if (task.contains("model"))         c.model         = task.at("model").get<std::string>();
    if (task.contains("mode"))          c.mode          = task.at("mode").get<std::string>();
    if (task.contains("Dfield_eV"))     c.Dfield_eV     = task.at("Dfield_eV").get<double>();
    if (task.contains("output_suffix")) c.output_suffix = task.at("output_suffix").get<std::string>();

    // ---------------- kmesh ----------------
    if (task.contains("kmesh")) {
        const auto& km = task.at("kmesh");
        if (km.contains("type"))    c.kmesh.type = km.at("type").get<std::string>();
        if (km.contains("Nk"))      c.kmesh.Nk   = km.at("Nk").get<int>();
        if (km.contains("dk_frac")) c.kmesh.dk_frac = km.at("dk_frac").get<double>();
    }

    // ---------------- temperature (always required) ----------------
    const auto& tlist = task.at("temperatureList_K");
    const double tmin = tlist.at("min").get<double>();
    const double tmax = tlist.at("max").get<double>();
    const int    nt   = tlist.at("num").get<int>();

    if (nt < 1)
        throw std::runtime_error("temperatureList_K.num must be >= 1");
    if (tmax < tmin)
        throw std::runtime_error("temperatureList_K.max must be >= min");

    c.temperature_list = la::linspace(tmin, tmax, nt);

    // ============================================================
    // mode = doping
    // ============================================================
    if (c.mode == "doping") {

        const auto& dlist = task.at("dopingList_10^12cm^-2");
        const double dmin = dlist.at("min").get<double>();
        const double dmax = dlist.at("max").get<double>();
        const int    nd   = dlist.at("num").get<int>();

        if (nd < 1)
            throw std::runtime_error("dopingList_10^12cm^-2.num must be >= 1");
        if (dmax < dmin)
            throw std::runtime_error("dopingList_10^12cm^-2.max must be >= min");

        c.doping_list = la::linspace(dmin, dmax, nd);

        // ---- doping -> filling ----
        const auto para = read_rg_para(c.para_file);
        const double a  = para.a0 - 0.00197 * para.pressure;

        c.filling_list.resize(c.doping_list.size());
        for (size_t i = 0; i < c.doping_list.size(); ++i) {
            c.filling_list[i] = la::doping_to_filling(c.doping_list[i], a);
        }

        // mu_list ignored
        c.mu_list.clear();
    }

    // ============================================================
    // mode = mu
    // ============================================================
    else if (c.mode == "mu") {

        const auto& mlist = task.at("muList_eV");
        const double mmin = mlist.at("min").get<double>();
        const double mmax = mlist.at("max").get<double>();
        const int    nm   = mlist.at("num").get<int>();

        if (nm < 1)
            throw std::runtime_error("muList_eV.num must be >= 1");
        if (mmax < mmin)
            throw std::runtime_error("muList_eV.max must be >= min");

        c.mu_list = la::linspace(mmin, mmax, nm);

        // doping / filling ignored
        c.doping_list.clear();
        c.filling_list.clear();
    }

    // ============================================================
    // unknown mode
    // ============================================================
    else {
        throw std::runtime_error(
            "calculate_fermi: unknown mode = \"" + c.mode +
            "\" (must be \"doping\" or \"mu\")"
        );
    }

    return c;
}


inline CalChiConfig read_cal_chi_config(const std::string& config_file)
{
    std::ifstream fin(config_file);
    if (!fin) throw std::runtime_error("Cannot open config file: " + config_file);

    json j;
    fin >> j;

    CalChiConfig c;

    // -------------------------
    // io
    // -------------------------
    if (j.contains("io")) {
        const auto& io = j.at("io");
        if (io.contains("data_dir"))   c.data_dir  = io.at("data_dir").get<std::string>();
        if (io.contains("model_para")) c.para_file = io.at("model_para").get<std::string>();
    }

    // -------------------------
    // tasks.calculate_chi
    // -------------------------
    const auto& task = j.at("tasks").at("calculate_chi");

    if (task.contains("model"))        c.model        = task.at("model").get<std::string>();
    if (task.contains("mode"))         c.mode         = task.at("mode").get<std::string>();
    if (task.contains("boundary"))     c.boundary     = task.at("boundary").get<std::string>();
    if (task.contains("Dfield_eV"))    c.Dfield_eV    = task.at("Dfield_eV").get<double>();
    if (task.contains("iq_min"))       c.iq_min       = task.at("iq_min").get<int>();
    if (task.contains("iq_max"))       c.iq_max       = task.at("iq_max").get<int>();
    if (task.contains("jq_min"))       c.jq_min       = task.at("jq_min").get<int>();
    if (task.contains("jq_max"))       c.jq_max       = task.at("jq_max").get<int>();
    if (task.contains("eta"))          c.eta          = task.at("eta").get<double>();
    if (task.contains("output_suffix"))
        c.output_suffix= task.at("output_suffix").get<std::string>();
    if (task.contains("fermi_patch_path")) 
        c.fermi_patch_path = task.at("fermi_patch_path").get<std::string>();

    // -------------------------
    // kmesh
    // -------------------------
    if (task.contains("kmesh")) {
        const auto& km = task.at("kmesh");
        if (km.contains("type"))       c.kmesh.type    = km.at("type").get<std::string>();
        if (km.contains("Nk"))         c.kmesh.Nk      = km.at("Nk").get<int>();
        if (km.contains("dk_frac"))    c.kmesh.dk_frac = km.at("dk_frac").get<double>();
    }

    if (!task.contains("temperatureList_K"))
        throw std::runtime_error("calculate_chi: missing temperatureList_K");

    {
        const auto& tlist = task.at("temperatureList_K");
        const double tmin = tlist.at("min").get<double>();
        const double tmax = tlist.at("max").get<double>();
        const int    nt   = tlist.at("num").get<int>();
        if (nt < 1) throw std::runtime_error("temperatureList_K.num must be >= 1");
        c.temperature_list = la::linspace(tmin, tmax, nt); // allow decreasing
    }

    if (!task.contains("polarList_meV"))
        throw std::runtime_error("calculate_chi: missing polarList_meV");

    {
        const auto& plist = task.at("polarList_meV");
        const double pmin = 0.001*plist.at("min").get<double>();
        const double pmax = 0.001*plist.at("max").get<double>();
        const int    np   = plist.at("num").get<int>();
        if (np < 1) throw std::runtime_error("polarList_meV.num must be >= 1");
        c.polar_list = la::linspace(pmin, pmax, np); // allow decreasing
    }


    // normalize mode
    if (c.mode.empty()) c.mode = "doping";

    // -------------------------
    // mode-dependent lists
    // -------------------------
    if (c.mode == "doping") {

        if (!task.contains("dopingList_10^12cm^-2"))
            throw std::runtime_error("calculate_chi: mode=doping but missing dopingList_10^12cm^-2");

        const auto& dlist = task.at("dopingList_10^12cm^-2");
        const double dmin = dlist.at("min").get<double>();
        const double dmax = dlist.at("max").get<double>();
        const int    nd   = dlist.at("num").get<int>();
        if (nd < 1) throw std::runtime_error("dopingList_10^12cm^-2.num must be >= 1");

        c.doping_list = la::linspace(dmin, dmax, nd); // allow decreasing

        // convert doping -> filling (metadata / fallback)
        const auto para = read_rg_para(c.para_file);
        const double a = para.a0 - 0.00197 * para.pressure;

        c.filling_list.resize(c.doping_list.size());
        for (size_t i = 0; i < c.doping_list.size(); ++i) {
            c.filling_list[i] = la::doping_to_filling(c.doping_list[i], a);
        }

        // mu_list in this mode: ignore (leave empty)
        c.mu_list.clear();

    } else if (c.mode == "mu") {

        if (!task.contains("muList_eV"))
            throw std::runtime_error("calculate_chi: mode=mu but missing muList_eV");

        const auto& mlist = task.at("muList_eV");
        const double mmin = mlist.at("min").get<double>();
        const double mmax = mlist.at("max").get<double>();
        const int    nm   = mlist.at("num").get<int>();
        if (nm < 1) throw std::runtime_error("muList_eV.num must be >= 1");

        c.mu_list = la::linspace(mmin, mmax, nm); // allow decreasing

        // doping/filling in this mode: not required (leave empty)
        c.doping_list.clear();
        c.filling_list.clear();

    } else {
        throw std::runtime_error("calculate_chi: unsupported mode = " + c.mode +
                                 " (must be \"doping\" or \"mu\")");
    }

    return c;
}


#pragma endregion
// ============================================================================

} // namespace config
