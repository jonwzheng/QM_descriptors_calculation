#!/bin/bash -l
## SBATCH -J COSMO
## SBATCH -n 1
## SBATCH -c 1
## SBATCH --array=0-1

echo "============================================================"
echo "Job ID : $SLURM_JOB_ID"
echo "Job Name : $SLURM_JOB_NAME"
echo "Starting on : $(date)"
echo "Running on node : $SLURMD_NODENAME"
echo "Current directory : $(pwd)"
echo "============================================================"

conda activate rdmc_env

#COSMO
TURBODIR=/home/gridsan/groups/RMG/Software/TmoleX19/TURBOMOLE
source $TURBODIR/Config_turbo_env
COSMOTHERM_PATH=/home/gridsan/groups/RMG/Software/COSMOtherm2021
COSMO_DATABASE_PATH=/home/gridsan/groups/RMG/COSMO_database/COSMObase2021

#QMD
QMD_PATH=/home/gridsan/groups/RMG/Software/QM_descriptors_calculation-radical_workflow
input_smiles=/home/gridsan/groups/RMG/Projects/Hao-Wei-Oscar-Yunsie/production_run/HAbs/inputs/TS_cosmo_dlpno/wb97xd_and_xtb_opted_ts_combo_results_hashed_sep1a_ts_input.csv
xyz_DFT_opt_dict=/home/gridsan/groups/RMG/Projects/Hao-Wei-Oscar-Yunsie/production_run/HAbs/inputs/TS_cosmo_dlpno/wb97xd_and_xtb_opted_ts_combo_results_hashed_sep1a_ts_dft_xyz.pkl

scratch_dir=$TMPDIR/$USER/$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID
mkdir -p $scratch_dir

python $QMD_PATH/main_COSMO_calc.py --input_smiles $input_smiles --xyz_DFT_opt_dict $xyz_DFT_opt_dict --scratch_dir $scratch_dir --COSMO_input_pure_solvents $QMD_PATH/common_solvent_list_final.csv --COSMOtherm_path $COSMOTHERM_PATH --COSMO_database_path $COSMO_DATABASE_PATH