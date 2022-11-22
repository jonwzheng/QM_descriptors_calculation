#!/bin/bash -l
#SBATCH --array=0-49

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

#inputs
input_smiles=reactants_products_wb97xd_and_xtb_opted_ts_combo_results_hashed_lookup_table_sep1a_filtered.csv
output_name=reactants_products_sep1a_filtered
xyz_FF_dict=${output_name}_ff_opted_results_xyz.pkl
xyz_semiempirical_opt_dict=${output_name}_semiempirical_opted_results_xyz.pkl

#ff conf
python $QMD_PATH/main_make_FF_conf_input.py --input_smiles $input_smiles --XTB_path $XTB_PATH 

# #semi opt
# python $QMD_PATH/main_make_semiempirical_opt_input.py --input_smiles $input_smiles --XTB_path $XTB_PATH --RDMC_path $RDMC_PATH --G16_path $g16root --xyz_FF_dict $xyz_FF_dict

# #dft opt freq
# python $QMD_PATH/main_make_DFT_opt_freq_input.py --input_smiles $input_smiles --G16_path $g16root --xyz_semiempirical_opt_dict $xyz_semiempirical_opt_dict
