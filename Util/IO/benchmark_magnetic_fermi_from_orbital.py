#!/usr/bin/env python3
import argparse
import math
import re
import struct
from pathlib import Path


MU_B_EV_PER_T = 5.7883818060e-5
KB_EV_PER_K = 8.617333262145e-5


def read_exact(f, n, what):
    b = f.read(n)
    if len(b) != n:
        raise EOFError(f"unexpected EOF while reading {what}")
    return b


def read_fermi_patch(path):
    path = Path(path)
    rows = {}

    with path.open("rb") as f:
        magic, version, nktot, dim = struct.unpack(
            "<iiii", read_exact(f, 16, "header ints")
        )

        if magic != 20260510:
            raise ValueError(f"{path}: bad magic {magic}")

        ef, t_k, doping, filling, dx, dy = struct.unpack(
            "<dddddd", read_exact(f, 48, "header doubles")
        )
        mesh_type = read_exact(f, 32, "mesh_type").split(b"\0", 1)[0].decode(
            errors="replace"
        )

        spin_sign = None
        doping_spin = None
        filling_spin = None

        if version == 4:
            read_exact(f, 4, "v4 shear flag")
        elif version >= 5:
            (spin_sign,) = struct.unpack("<i", read_exact(f, 4, "v5 spin_sign"))
            doping_spin, filling_spin = struct.unpack(
                "<dd", read_exact(f, 16, "v5 spin metadata")
            )

        for _ in range(nktot):
            iq, jq = struct.unpack("<ii", read_exact(f, 8, "iq jq"))
            kx, ky, occ_avg = struct.unpack("<ddd", read_exact(f, 24, "k/occ"))
            evals = struct.unpack("<" + "d" * dim, read_exact(f, 8 * dim, "evals"))
            occ = struct.unpack("<" + "d" * dim, read_exact(f, 8 * dim, "occ"))
            read_exact(f, 16 * dim * dim, "evecs")

            rows[(iq, jq)] = {
                "kx": kx,
                "ky": ky,
                "occ_avg": occ_avg,
                "evals": evals,
                "occ": occ,
            }

    return {
        "path": str(path),
        "version": version,
        "NkTot": nktot,
        "dim": dim,
        "EF": ef,
        "T_K": t_k,
        "doping": doping,
        "filling": filling,
        "dx": dx,
        "dy": dy,
        "mesh_type": mesh_type,
        "spin_sign": spin_sign,
        "doping_spin": doping_spin,
        "filling_spin": filling_spin,
        "rows": rows,
    }


def parse_orbital_moment(path):
    path = Path(path)
    bands = []
    rows = {}

    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith("# columns ="):
                cols = line.split("=", 1)[1].strip().split()
                for col in cols:
                    m = re.fullmatch(r"m_orb_muB_b(\d+)", col)
                    if m:
                        bands.append(int(m.group(1)))
                continue

            if line.startswith("#"):
                continue

            vals = line.split()
            iq = int(vals[1])
            jq = int(vals[2])
            kx = float(vals[3])
            ky = float(vals[4])

            idx = 5
            m_by_band = {}
            for b in bands:
                idx += 1  # E_b*
                m_by_band[b] = float(vals[idx])
                idx += 1

            rows[(iq, jq)] = {
                "kx": kx,
                "ky": ky,
                "m_muB": m_by_band,
            }

    if not bands:
        raise ValueError(f"{path}: cannot find m_orb_muB_b* columns")

    return {
        "path": str(path),
        "bands": bands,
        "rows": rows,
    }


def fermi_occ(E, mu, T_K):
    if T_K <= 0.0:
        return 1.0 if E <= mu else 0.0

    x = (E - mu) / (KB_EV_PER_K * T_K)

    if x > 50.0:
        return 0.0
    if x < -50.0:
        return 1.0

    return 1.0 / (math.exp(x) + 1.0)


def stats(xs):
    if not xs:
        return {
            "n": 0,
            "mean": float("nan"),
            "mean_abs": float("nan"),
            "rms": float("nan"),
            "max_abs": float("nan"),
        }

    n = len(xs)
    mean = sum(xs) / n
    mean_abs = sum(abs(x) for x in xs) / n
    rms = math.sqrt(sum(x * x for x in xs) / n)
    max_abs = max(abs(x) for x in xs)

    return {
        "n": n,
        "mean": mean,
        "mean_abs": mean_abs,
        "rms": rms,
        "max_abs": max_abs,
    }


def compare_spin(b0, direct, orbital, spin_sign, B_T, g_factor):
    dE = []
    dOcc = []
    missing = 0

    for key, row0 in b0["rows"].items():
        row_direct = direct["rows"].get(key)
        row_m = orbital["rows"].get(key)

        if row_direct is None or row_m is None:
            missing += 1
            continue

        for band in orbital["bands"]:
            if band >= b0["dim"] or band >= direct["dim"]:
                continue

            m_muB = row_m["m_muB"][band]
            E0 = row0["evals"][band]
            E_approx = (
                E0
                - m_muB * MU_B_EV_PER_T * B_T
                + spin_sign * g_factor * MU_B_EV_PER_T * B_T
            )
            E_direct = row_direct["evals"][band]

            occ_approx = fermi_occ(E_approx, b0["EF"], b0["T_K"])
            occ_direct = row_direct["occ"][band]

            dE.append(E_approx - E_direct)
            dOcc.append(occ_approx - occ_direct)

    return {
        "missing_k": missing,
        "dE_eV": stats(dE),
        "dE_meV": stats([x * 1000.0 for x in dE]),
        "dOcc": stats(dOcc),
    }


def print_stats(label, result):
    print(f"\n=== {label} ===")
    print(f"missing_k = {result['missing_k']}")

    for name in ["dE_meV", "dOcc"]:
        s = result[name]
        print(
            f"{name}: n={s['n']} mean={s['mean']:.6g} "
            f"mean_abs={s['mean_abs']:.6g} rms={s['rms']:.6g} "
            f"max_abs={s['max_abs']:.6g}"
        )


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Benchmark direct magnetic fermiPatch bins against applying "
            "orbital moment + Zeeman shifts to a B=0 fermiPatch."
        )
    )
    ap.add_argument("--b0", required=True, help="B=0 fermiPatch.bin")
    ap.add_argument("--up", required=True, help="direct B fermi_spin_up_patch.bin")
    ap.add_argument("--down", required=True, help="direct B fermi_spin_down_patch.bin")
    ap.add_argument("--orbital", required=True, help="orbital_moment txt file")
    ap.add_argument("--B", type=float, default=1.0, help="magnetic field in Tesla")
    ap.add_argument("--g", type=float, default=2.0, help="Zeeman g factor")
    args = ap.parse_args()

    b0 = read_fermi_patch(args.b0)
    up = read_fermi_patch(args.up)
    down = read_fermi_patch(args.down)
    orbital = parse_orbital_moment(args.orbital)

    print("B0     :", b0["path"])
    print("up     :", up["path"])
    print("down   :", down["path"])
    print("orbital:", orbital["path"])
    print("bands  :", orbital["bands"])
    print(f"B_T={args.B} g={args.g} muB={MU_B_EV_PER_T} eV/T")
    print(f"mu={b0['EF']} eV T={b0['T_K']} K Nk={b0['NkTot']} dim={b0['dim']}")

    up_result = compare_spin(b0, up, orbital, +1, args.B, args.g)
    down_result = compare_spin(b0, down, orbital, -1, args.B, args.g)

    print_stats("spin up: approx(B0 + orbital + Zeeman) - direct B", up_result)
    print_stats("spin down: approx(B0 + orbital + Zeeman) - direct B", down_result)


if __name__ == "__main__":
    main()
