#!/bin/bash -l
#SBATCH --array=0-49

echo "============================================================"
echo "Job ID : $SLURM_JOB_ID"
echo "Job Name : $SLURM_JOB_NAME"
echo "Starting on : $(date)"
echo "Running on node : $SLURMD_NODENAME"
echo "Current directory : $(pwd)"
echo "============================================================"

# conda activate rdmc_env

#RDMC for gaussian-xtb
RDMC_PATH=/home/gridsan/groups/RMG/Software/RDMC-main
export PATH=$RDMC_PATH:$PATH
export PYTHONPATH=$RDMC_PATH:$PYTHONPATH

#QMD
QMD_PATH=/home/gridsan/groups/RMG/Software/QM_descriptors_calculation-radical_workflow
input_smiles=/home/gridsan/groups/RMG/Projects/Hao-Wei-Oscar-Yunsie/production_run/HAbs/inputs/TS_sep1a_all/wb97xd_and_xtb_opted_ts_combo_results_hashed_sep1a_ts_input.csv
xyz_DFT_opt_dict=/home/gridsan/groups/RMG/Projects/Hao-Wei-Oscar-Yunsie/production_run/HAbs/inputs/TS_sep1a_all/wb97xd_and_xtb_opted_ts_combo_results_hashed_sep1a_ts_dft_xyz.pkl

#r p complex semi opt
python $QMD_PATH/main_make_r_p_complex_FF_opt_input.py --input_smiles $input_smiles --RDMC_path $RDMC_PATH --xyz_DFT_opt_dict $xyz_DFT_opt_dict

