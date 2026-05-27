#!/bin/bash
set -euo pipefail

# ===================== FIXED PROJECT ROOT =====================
PROJECT_ROOT="/pscratch/sd/h/hryan/rg_main"
cd "$PROJECT_ROOT"

# ===================== USER SETTINGS =====================
QOS="regular"
WALLTIME_REGULAR="2:00:00"
WALLTIME_DEBUG="0:30:00"

# One Slurm job is staged for each (T, B) pair.
# If T_VALUES is non-empty, it overrides T_MIN/T_MAX/T_NUM.
T_VALUES=""
T_MIN=0
T_MAX=10
T_NUM=51

MU_MIN=0.5
MU_MAX=0.8
MU_NUM=301

B_VALUES=""
B_MIN=1
B_MAX=1
B_NUM=1
G_FACTOR=2.0
ORBITAL_DK=1e-5

DFIELD=-0.01
MODEL="sk"

KMESH_TYPE="b1b2"
NK=60
DK_FRAC=0.00125
OUTPUT_SUFFIX="auto"

EXE_ABS="${PROJECT_ROOT}/bin/calculate_fermi"
CONFIG_TEMPLATE="${PROJECT_ROOT}/config.json"
RG_PARA_SRC="${PROJECT_ROOT}/Source/Parameters/RG_para.json"
SBATCH_TEMPLATE="${PROJECT_ROOT}/submit_fermi_oneT.slurm"

DO_SUBMIT=1
SUBMIT_PARALLEL=1
MAX_SUBMIT_JOBS=8

RUN_TAG="$(date +"%Y-%m-%d_%H%M%S")"
OUT_PREFIX="fermi"
[ "$QOS" = "debug" ] && OUT_PREFIX="debug_fermi"

RUN_ROOT="${PROJECT_ROOT}/data/${OUT_PREFIX}_${RUN_TAG}"
JOBS_DIR="${RUN_ROOT}/jobs"
FERMI_DATA_BASE="${RUN_ROOT}"

# ===================== helpers =====================
linspace () {
  awk -v min="$1" -v max="$2" -v num="$3" '
    BEGIN {
      if (num == 1) { printf "%.6f\n", min; exit }
      for (i = 0; i < num; i++)
        printf "%.6f\n", min + (max - min) * i / (num - 1)
    }'
}

tag_T () {
  awk -v T="$1" 'BEGIN { printf "T%.3f", T }'
}

tag_B () {
  awk -v B="$1" 'BEGIN {
    s = sprintf("%.6f", B)
    sub(/0+$/, "", s)
    sub(/\.$/, "", s)
    if (s == "-0") s = "0"
    printf "B%sT", s
  }'
}

throttle_bg () {
  local max_jobs="$1"
  while true; do
    local nj
    nj=$(jobs -rp | wc -l | awk '{print $1}')
    [ "$nj" -lt "$max_jobs" ] && break
    sleep 0.2
  done
}

# ===================== sanity =====================
for f in "$EXE_ABS" "$CONFIG_TEMPLATE" "$RG_PARA_SRC" "$SBATCH_TEMPLATE"; do
  [ -f "$f" ] || { echo "[ERROR] missing: $f"; exit 1; }
done

WALLTIME="$WALLTIME_REGULAR"
[ "$QOS" = "debug" ] && WALLTIME="$WALLTIME_DEBUG"

mkdir -p "$JOBS_DIR" "$FERMI_DATA_BASE"

if [ -n "$T_VALUES" ]; then
  T_LIST="$T_VALUES"
  T_JOB_NUM="$(printf '%s\n' $T_VALUES | wc -l | awk '{print $1}')"
else
  T_LIST="$(linspace "$T_MIN" "$T_MAX" "$T_NUM")"
  T_JOB_NUM="$T_NUM"
fi

if [ -n "$B_VALUES" ]; then
  B_LIST="$B_VALUES"
  B_JOB_NUM="$(printf '%s\n' $B_VALUES | wc -l | awk '{print $1}')"
else
  B_LIST="$(linspace "$B_MIN" "$B_MAX" "$B_NUM")"
  B_JOB_NUM="$B_NUM"
fi

echo "[INFO] RUN_ROOT        : $RUN_ROOT"
echo "[INFO] JOBS_DIR        : $JOBS_DIR"
echo "[INFO] FERMI_DATA_BASE : $FERMI_DATA_BASE"
echo "[INFO] T values        : $T_LIST"
echo "[INFO] B values        : $B_LIST"
echo "[INFO] Jobs            : $((T_JOB_NUM * B_JOB_NUM)) = T_count * B_count"
echo "[INFO] Expect outputs  : ${FERMI_DATA_BASE}/D*/B*/T*/mu*/fermi*.bin"

for T in $T_LIST; do
  TAGT="$(tag_T "$T")"
  for B in $B_LIST; do
    TAGB="$(tag_B "$B")"
    JOBDIR="${JOBS_DIR}/${TAGT}/${TAGB}"
    mkdir -p "$JOBDIR"

    echo "  -> stage $JOBDIR"

    cp "$EXE_ABS"         "$JOBDIR/"
    cp "$CONFIG_TEMPLATE" "$JOBDIR/config.json"
    cp "$RG_PARA_SRC"     "$JOBDIR/RG_para.json"

    python3 - <<PY
import json

cfg_path = r"$JOBDIR/config.json"

with open(cfg_path, "r") as f:
    cfg = json.load(f)

cfg.setdefault("io", {})
cfg["io"]["data_dir"] = r"$FERMI_DATA_BASE"
cfg["io"]["model_para"] = "RG_para.json"

fermi = cfg.setdefault("calculate_fermi", {})

fermi["model"] = r"$MODEL"
fermi["Dfield_eV"] = float(r"$DFIELD")
fermi["magnetic_field"] = {
    "BList_T": {
        "min": float(r"$B"),
        "max": float(r"$B"),
        "num": 1
    },
    "g_factor": float(r"$G_FACTOR"),
    "orbital_derivative_dk": float(r"$ORBITAL_DK")
}

fermi["kmesh"] = {
    "type": r"$KMESH_TYPE",
    "Nk": int(r"$NK"),
    "dk_frac": float(r"$DK_FRAC")
}

fermi["temperatureList_K"] = {
    "min": float(r"$T"),
    "max": float(r"$T"),
    "num": 1
}

fermi["muList_eV"] = {
    "min": float(r"$MU_MIN"),
    "max": float(r"$MU_MAX"),
    "num": int(r"$MU_NUM")
}

fermi["output_suffix"] = r"$OUTPUT_SUFFIX"

with open(cfg_path, "w") as f:
    json.dump(cfg, f, indent=2)
    f.write("\n")
PY

    sed \
      -e "s/{{JOB_NAME}}/fermi_${TAGT}_${TAGB}/g" \
      -e "s|{{EXE_ABS}}|${JOBDIR}/$(basename "$EXE_ABS")|g" \
      -e "s|{{CFG_ABS}}|${JOBDIR}/config.json|g" \
      -e "s/^#SBATCH[[:space:]]\+-q[[:space:]]\+.*/#SBATCH -q ${QOS}/" \
      -e "s/^#SBATCH[[:space:]]\+-t[[:space:]]\+.*/#SBATCH -t ${WALLTIME}/" \
      "$SBATCH_TEMPLATE" > "$JOBDIR/run.sh"

    chmod +x "$JOBDIR/run.sh"

    if [ "$DO_SUBMIT" -eq 1 ]; then
      if [ "$SUBMIT_PARALLEL" -eq 1 ]; then
        throttle_bg "$MAX_SUBMIT_JOBS"
        ( cd "$JOBDIR" && sbatch run.sh ) &
      else
        ( cd "$JOBDIR" && sbatch run.sh )
      fi
    fi
  done
done

[ "$DO_SUBMIT" -eq 1 ] && [ "$SUBMIT_PARALLEL" -eq 1 ] && wait

echo "[DONE] Jobs staged under      : $JOBS_DIR"
echo "[DONE] All fermi outputs under: $FERMI_DATA_BASE"
