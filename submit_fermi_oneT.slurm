#!/bin/bash
#SBATCH -N 4
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -J {{JOB_NAME}}
#SBATCH --mail-type=ALL
#SBATCH -A m3706
#SBATCH -t 0:30:00
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1

module purge
module load craype-x86-milan
module load craype-network-ofi
module load libfabric/1.22.0
module load PrgEnv-gnu/8.5.0
module load cray-mpich/8.1.30
module load cudatoolkit/12.4

echo "[INFO] Job start : $(date)"
echo "[INFO] Job name  : $SLURM_JOB_NAME"
echo "[INFO] Work dir  : $(pwd)"
echo "[INFO] Nodes     : ${SLURM_JOB_NUM_NODES:-unknown}"
echo "[INFO] Tasks     : ${SLURM_NTASKS:-unknown}"
module list 2>&1

export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

t0=$SECONDS
srun -n "${SLURM_NTASKS}" -c "${SLURM_CPUS_PER_TASK}" --cpu_bind=cores {{EXE_ABS}} {{CFG_ABS}}
dt=$((SECONDS - t0))

echo "[INFO] Job end   : $(date)"
printf "[INFO] Elapsed   : %02d:%02d:%02d (hh:mm:ss)\n" $((dt/3600)) $(((dt%3600)/60)) $((dt%60))
