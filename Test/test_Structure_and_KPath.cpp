#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

#include "PhysStruct/RG_Structure.h"
#include "Util/IO/rg_io.h"
#include "Util/Config/config_reader.h"

int main(int argc, char** argv) {
    try {
        std::cout << std::fixed << std::setprecision(8);

        const std::string config_file =
            (argc > 1) ? std::string(argv[1]) : "config.json";

        // ====================================================
        // Read config for this test
        // ====================================================
        const auto cfg = config::read_test_structure_config(config_file);

        if (cfg.model != "sk") {
            throw std::runtime_error(
                "test_Structure_and_KPath only supports model = \"sk\""
            );
        }

        // ====================================================
        // Read RG_para.json (SK model constants)
        // ====================================================
        const auto para = config::read_rg_para(cfg.para_file);

        // ====================================================
        // Build structure
        // ====================================================
        const double pressure = 0.0;  // fixed for this test
        rg::RG_Structure st(
            para.layer_num,
            pressure,
            para.a0,
            para.d0,
            para.vacuum
        );

        // ====================================================
        // Console summary
        // ====================================================
        std::cout << "=== test_Structure_and_KPath ===\n";
        std::cout << "config   = " << config_file << "\n";
        std::cout << "para     = " << cfg.para_file << "\n";
        std::cout << "model    = " << cfg.model << "\n";
        std::cout << "layers   = " << para.layer_num << "\n";
        std::cout << "Nk_seg   = " << cfg.Nk_seg << "\n";

        rgio::cout_vec2("a1", st.a1());
        rgio::cout_vec2("a2", st.a2());
        rgio::cout_vec2("b1", st.b1());
        rgio::cout_vec2("b2", st.b2());
        rgio::cout_vec2("K",  st.K());
        std::cout << "\n";

        // ====================================================
        // Output (single dump file)
        // ====================================================
        const std::string out_file =
            cfg.data_dir + "/test_Structure_and_KPath_dump.txt";

        std::ofstream f(out_file);
        if (!f)
            throw std::runtime_error("Cannot write to " + out_file);

        f << std::setprecision(15);

        f << "\n# === STRUCTURE_META ===\n";
        f << "# nLayer a0 d0 vacuum\n";
        f << para.layer_num << " "
          << para.a0 << " "
          << para.d0 << " "
          << para.vacuum << "\n";

        f << "\n# === ATOMS ===\n";
        f << "# x y z layer sub\n";
        for (const auto& s : st.atoms()) {
            int sub = (s.sub == rg::Sublattice::A ? 0 : 1);
            f << s.r.x() << " " << s.r.y() << " " << s.r.z() << " "
              << s.layer << " " << sub << "\n";
        }

        f << "\n# === REAL_CELL ===\n";
        f << st.a1().x() << " " << st.a1().y() << " "
          << st.a2().x() << " " << st.a2().y() << "\n";

        f << "\n# === KPATH_GMKG ===\n";
        const auto path = st.generate_GMKG(cfg.Nk_seg);
        for (size_t i = 0; i < path.k_list.size(); ++i) {
            f << i << " "
              << path.k_list[i].x() << " "
              << path.k_list[i].y() << " "
              << path.kline[i] << "\n";
        }

        f.close();

        std::cout << "Saved: " << out_file << "\n";
        std::cout << "DONE.\n";
        return 0;

    } catch (const std::exception& e) {
        std::cerr << "\n[ERROR] " << e.what() << "\n";
        return 1;
    }
}
