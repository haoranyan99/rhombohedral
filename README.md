# Rhombohedral Graphene Project

This project provides a modular C++ framework for band structure, density of states (DOS),
and valley susceptibility calculations in **rhombohedral (ABC-stacked) multilayer graphene**.

The code is designed for **theoretical and numerical studies of low-energy electronic
structure**, with a focus on **valley physics, Fermi surface geometry, and susceptibility
χ(q)**. It is portable across **macOS** and **Windows (MinGW)**, and supports both
production calculations and systematic validation tests.

---

## Features

- Slater–Koster (SK) tight-binding model for multilayer graphene
- 10-band SWMcC \(k\cdot p\) (KP10) model
- Band structure along high-symmetry paths and local-K patches
- Density of states (DOS) calculations
- Valley susceptibility χ(q) on **local-K hexagonal patches**
- Patch-based Fermi energy solving and Fermi-surface diagnostics
- Unified JSON-based parameter control

---

## Project Structure

```text
rhombohedral_project/
│
├─ Makefile
├─ config.json                    # ⭐ Unified runtime configuration file
│
├─ include/
│  ├─ eigen-3.4.0/                # Eigen linear algebra library
│  └─ nlohmann/
│     └─ json.hpp                 # JSON parsing
│
├─ Executables/                   # ===== Main programs =====
│  ├─ rg_band_dos_valley.cpp      # Band structure + valley-resolved DOS
│  └─ rg_valley_chi.cpp           # Valley susceptibility χ(q)
│
├─ Test/                          # ===== Test & validation programs =====
│  ├─ test_all.cpp                # Run / collect all tests
│  ├─ test_RG_KP_bareband.cpp     # Band structure test for KP model
│  ├─ test_RG_SK_bareband.cpp     # Band structure test for SK model
│  ├─ test_RG_chi_hexpatch.cpp    # χ(q) on local-K hex patch (SK model)
│  ├─ test_dos_total_triangle.cpp # Global DOS (triangle method) regression
│  ├─ test_gen_hopp.cpp           # Real-space hopping generation test
│  ├─ test_read_vhs_from_dos.cpp  # VHS peak extraction from DOS
│  ├─ test_Structure_and_KPath.cpp# Geometry / stacking / lattice checks
│  └─ test_twodos.cpp             # Analytic 2D DOS (parabolic band) test
│
├─ data/                          # ===== Output data =====
│  └─ *.txt                       # Band / DOS / χ / FS diagnostics
│
├─ Common/                        # ===== Shared data containers =====
│  └─ DataContainers.h            # GridData, field storage, indexing
│
├─ PhysStruct/                    # ===== Geometry & k-space =====
│  └─ RG_Structure.h              # Lattice, stacking, k-mesh, K-paths
│
├─ Source/
│  ├─ Models/                     # ===== Electronic models =====
│  │  ├─ RG_KPModel.h             # 10-band SWMcC k·p model
│  │  └─ RG_SKModel.h             # Slater–Koster tight-binding model
│  │
│  └─ Parameters/                 # ===== Model parameters =====
│     └─ RG_ModelParams.json      # SK / KP parameters + χ numerical settings
│
├─ Util/
│  ├─ IO/
│  │  └─ rg_io.h
│  │
│  ├─ DOS/                        # DOS-related utilities (future extension)
│  └─ Config/                     # Runtime / JSON helpers (future extension)
│
├─ MeasureEngines/                # ===== Numerical engines =====
│  └─ (reserved for solvers / integrators)
