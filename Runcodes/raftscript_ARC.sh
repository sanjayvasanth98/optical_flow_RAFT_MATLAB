#!/bin/bash
#SBATCH --output=RAFT_SideView_%j.out
#SBATCH --error=RAFT_SideView_%j.err
#SBATCH --time=02:00:00
#SBATCH --account=cavitation
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:1
#SBATCH --mem=8G
#SBATCH --partition=a30_normal_q
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=youremailid@vt.edu
#SBATCH --job-name=jobname

module reset
module load MATLAB

# Ensure job starts from submission directory
cd $SLURM_SUBMIT_DIR

echo "============================================================"
echo " SLURM Job ID:        $SLURM_JOB_ID"
echo " Node:                $HOSTNAME"
echo " Submitted from:      $SLURM_SUBMIT_DIR"
echo " Start time:          $(date)"
echo "============================================================"

# GPU monitoring (every 300 seconds)
nvidia-smi --query-gpu=timestamp,name,pci.bus_id,driver_version,temperature.gpu,utilization.gpu,utilization.memory,memory.total,memory.free,memory.used --format=csv -l 300 > gpu.perform.$SLURM_JOBID.log &

# -------------------------------------------------------------------------
# Launch MATLAB (non-interactive batch mode)
# -------------------------------------------------------------------------
matlab -nodisplay -nosplash -nodesktop -batch "run('/home/kbsanjayvasanth/Inception_raft_test/smooth_naga/raftmatlabsideview_ARC.m')"

echo "============================================================"
echo " Job finished at:     $(date)"
echo "============================================================"
