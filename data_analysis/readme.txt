RAFT Post-Processing & Data Analysis (ARC/Owl Cluster)
======================================================

This README provides instructions for running the RAFT post-processing scripts on the Virginia Tech ARC (Owl) cluster. The workflow generates cluster-safe, memory-efficient plots for velocity profiles and Reynolds stresses.

FILES OVERVIEW
--------------

1. RAFT_Data_analysis.m
   - Purpose: The main MATLAB computation engine.
   - Features:
     * Generates tiled velocity profiles (Throat + N downstream stations).
     * Calculates and plots Reynolds Stresses (uu, vv, uv).
     * Uses 'matfile' chunking to process large datasets without crashing RAM.
   - Outputs: .png plots and .mat data files.

2. raft_postprocessing_script_ARC.sh
   - Purpose: The Slurm submission script to run the analysis on the ARC cluster.
   - Features: Handles module loading, temporary directory management, and compute node allocation.


CONFIGURATION STEPS
-------------------

Before submitting the job, you MUST edit both files to match your specific case paths.

1. Edit the MATLAB Script (RAFT_Data_analysis.m)
   Open the file and locate the "USER INPUTS" section at the top.

   * matPaths: Provide the full absolute paths to your *_velocity.mat files.
   * caseLabels: Provide a short string name for each case (must match the order of matPaths).
   * userOutDir: The directory where you want the resulting images and data to be saved.
   * dx_mm: The spacing (in mm) between downstream profile stations.
   * nStations: How many profiles to plot after the throat.

   Example:
     matPaths = { "/groups/cavitation/project_A/Run1_velocity.mat" };
     caseLabels = { "Run 1 (Baseline)" };
     userOutDir = "/groups/cavitation/project_A/results";

2. Edit the Slurm Script (raft_postprocessing_script_ARC.sh)
   Open the file and update the following:

   * Email: #SBATCH --mail-user=yourPID@vt.edu
   * Jobname: #SBATCH --job-name=jobname
   * MATLAB Execution Path (CRITICAL):
     Look for the "srun matlab ..." line near the bottom. Update the path inside the run(...) command:
     srun matlab ... -batch "try, run('/path/to/your/RAFT_Data_analysis.m'); ..."


HOW TO RUN
----------

1. Upload both files to your working directory on the cluster.
2. Open 'Owl' cluster from the working directory.
3. Make the shell script executable (optional, can do step 3 directly):
   chmod +x raft_postprocessing_script_ARC.sh
3. Submit the job:
   sbatch raft_postprocessing_script_ARC.sh
4. Check status:
   squeue   
   (or)
   squeue -u yourPID  %specific line in queue


OUTPUTS
-------

The script creates the following structure inside your defined 'userOutDir':

* Velocity profiles/
  - VelProfile_Vmag_*.png: Comparison of magnitude profiles.
  - VelProfile_U_*.png: Comparison of horizontal velocity (u).
  - VelProfile_V_*.png: Comparison of vertical velocity (v).

* Reynolds Stresses/
  - [CaseName]_ReStress_uu.png: Normal stress map (uu).
  - [CaseName]_ReStress_vv.png: Normal stress map (vv).
  - [CaseName]_ReStress_uv.png: Shear stress map (uv).
  - [CaseName]_ReynoldsStresses.mat: Raw calculated stress arrays.


IMPORTANT NOTES
---------------

* Prerequisites: Input .mat files must contain 'u_all', 'v_all', 'mm_per_pixel', 'maskROI', and 'x_throat_mm'.
* Memory Efficiency: The script processes data in chunks (default 100 frames) to handle large 3D arrays.
* Paths: Always use absolute paths (starting with /home/ or /groups/) to avoid file-not-found errors.