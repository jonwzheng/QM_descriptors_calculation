#!/bin/bash -l

echo "============================================================"
echo "Job ID : $SLURM_JOB_ID"
echo "Job Name : $SLURM_JOB_NAME"
echo "Starting on : $(date)"
echo "Running on node : $SLURMD_NODENAME"
echo "Current directory : $(pwd)"
echo "============================================================"

conda activate rdmc_env

#xtb
source /home/gridsan/groups/RMG/Software/xtb-6.4.1/share/xtb/config_env.bash
XTB_PATH=/home/gridsan/groups/RMG/Software/xtb-6.4.1/bin
export PATH=$XTB_PATH:$PATH

#RDMC for gaussian-xtb
RDMC_PATH=/home/gridsan/groups/RMG/Software/RDMC-main

#gaussian
export PATH=$PATH:/home/gridsan/groups/RMG/Software/gaussian/g16
export PATH=$PATH:/home/gridsan/groups/RMG/Software/gaussian/gv
export g16root=/home/gridsan/groups/RMG/Software/gaussian
source /home/gridsan/groups/RMG/Software/gaussian/g16/bsd/g16.profile
export GAUSS_SCRDIR=""

#QMD
QMD_PATH=/home/gridsan/groups/RMG/Software/QM_descriptors_calculation-radical_workflow
#input_smiles=/home/gridsan/groups/RMG/Projects/Hao-Wei-Oscar-Yunsie/production_run/HAbs/inputs/TS_cosmo_dlpno/wb97xd_and_xtb_opted_ts_combo_results_hashed_sep1a_ts_input.csv

# #ff conf
# python $QMD_PATH/main_make_FF_conf_input.py --input_smiles $input_smiles --XTB_path $XTB_PATH 

#semi opt
python $QMD_PATH/main_make_semiempirical_opt_input.py --input_smiles $input_smiles --XTB_path $XTB_PATH --RDMC_path $RDMC_PATH --G16_path $g16root