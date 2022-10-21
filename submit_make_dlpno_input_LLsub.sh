#!/bin/bash -l

echo "============================================================"
echo "Job ID : $SLURM_JOB_ID"
echo "Job Name : $SLURM_JOB_NAME"
echo "Starting on : $(date)"
echo "Running on node : $SLURMD_NODENAME"
echo "Current directory : $(pwd)"
echo "============================================================"

#QMD
QMD_PATH=/home/gridsan/groups/RMG/Software/QM_descriptors_calculation-radical_workflow
input_smiles=/home/gridsan/groups/RMG/Projects/Hao-Wei-Oscar-Yunsie/production_run/HAbs/inputs/TS_cosmo_dlpno/wb97xd_and_xtb_opted_ts_combo_results_hashed_sep1a_ts_input.csv
xyz_DFT_opt=/home/gridsan/groups/RMG/Projects/Hao-Wei-Oscar-Yunsie/production_run/HAbs/inputs/TS_cosmo_dlpno/wb97xd_and_xtb_opted_ts_combo_results_hashed_sep1a_ts_dft_xyz.pkl

python $QMD_PATH/main_make_input.py --input_smiles $input_smiles --xyz_DFT_opt $xyz_DFT_opt --DLPNO_sp_n_procs 22 --DLPNO_sp_job_ram 3900
