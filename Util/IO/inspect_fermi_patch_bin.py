#!/usr/bin/env python3
import argparse
import struct
from pathlib import Path


def read_exact(f, n, what):
    b = f.read(n)
    if len(b) != n:
        raise EOFError(f"unexpected EOF while reading {what}")
    return b


def inspect(path, nrows):
    path = Path(path)
    with path.open("rb") as f:
        magic, version, NkTot, dim = struct.unpack("<iiii", read_exact(f, 16, "header ints"))
        EF, T_K, doping, filling, dx, dy = struct.unpack("<dddddd", read_exact(f, 48, "header doubles"))
        mesh_type = read_exact(f, 32, "mesh_type").split(b"\0", 1)[0].decode(errors="replace")

        out = {
            "file": str(path),
            "magic": magic,
            "version": version,
            "NkTot": NkTot,
            "dim": dim,
            "EF": EF,
            "T_K": T_K,
            "doping_header": doping,
            "filling_header": filling,
            "dx": dx,
            "dy": dy,
            "mesh_type": mesh_type,
        }

        if version == 4:
            has_shear, = struct.unpack("<i", read_exact(f, 4, "v4 shear flag"))
            out["has_legacy_shear_projection"] = has_shear
        elif version >= 5:
            spin_sign, = struct.unpack("<i", read_exact(f, 4, "v5 spin_sign"))
            doping_spin, filling_spin = struct.unpack("<dd", read_exact(f, 16, "v5 spin metadata"))
            out["spin_sign"] = spin_sign
            out["doping_total"] = doping
            out["filling_total"] = filling
            out["doping_spin"] = doping_spin
            out["filling_spin"] = filling_spin

        print("=== fermiPatch header ===")
        for k, v in out.items():
            print(f"{k}: {v}")

        if nrows <= 0:
            return

        print("\n=== first k rows ===")
        print("ik iq jq kx ky occ_avg evals[0:min(4,dim)] occ[0:min(4,dim)]")
        for ik in range(min(nrows, NkTot)):
            iq, jq = struct.unpack("<ii", read_exact(f, 8, "iq jq"))
            kx, ky, occ_avg = struct.unpack("<ddd", read_exact(f, 24, "k/occ"))
            evals = struct.unpack("<" + "d" * dim, read_exact(f, 8 * dim, "evals"))
            occ = struct.unpack("<" + "d" * dim, read_exact(f, 8 * dim, "occ_band"))
            f.seek(16 * dim * dim, 1)
            print(
                ik,
                iq,
                jq,
                f"{kx:.12g}",
                f"{ky:.12g}",
                f"{occ_avg:.12g}",
                [round(x, 10) for x in evals[: min(4, dim)]],
                [round(x, 10) for x in occ[: min(4, dim)]],
            )


def main():
    ap = argparse.ArgumentParser(description="Inspect fermiPatch binary headers and first rows.")
    ap.add_argument("bin_file")
    ap.add_argument("--rows", type=int, default=0, help="number of k rows to print")
    args = ap.parse_args()
    inspect(args.bin_file, args.rows)


if __name__ == "__main__":
    main()
