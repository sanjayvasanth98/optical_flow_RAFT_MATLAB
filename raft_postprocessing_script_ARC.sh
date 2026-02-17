#!/bin/bash
#SBATCH --output=RAFT_SideView_%j.out
#SBATCH --error=RAFT_SideView_%j.err
#SBATCH --time=02:00:00
#SBATCH --account=cavitation
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --mem=6G
#SBATCH --partition=normal_q
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kbsanjayvasanth@vt.edu
#SBATCH --job-name=test

set -euo pipefail

module reset
module load MATLAB

cd "$SLURM_SUBMIT_DIR"

echo "============================================================"
echo " SLURM Job ID:   $SLURM_JOB_ID"
echo " Node:           $HOSTNAME"
echo " Workdir:        $SLURM_SUBMIT_DIR"
echo " Start time:     $(date)"
echo " Cores:          $SLURM_CPUS_PER_TASK"
echo "============================================================"

# Keep MATLAB temp files in your job space (optional but recommended)
export TMPDIR="$SLURM_SUBMIT_DIR/tmp_${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"

# Prevent oversubscription / ensure MATLAB uses the allocated cores (important)
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Run under srun so Slurm properly tracks the step
srun matlab -nodisplay -nosplash -nodesktop -batch "\
try, run('/home/kbsanjayvasanth/Inception_raft_test/data_analysis/RAFT_Data_analysis.m'); \
catch ME, disp(getReport(ME,'extended')); exit(1); end; exit(0);"

echo "============================================================"
echo " Job finished at: $(date)"
echo "============================================================"
