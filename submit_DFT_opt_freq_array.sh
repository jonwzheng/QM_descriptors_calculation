#!/bin/bash -l
##SBATCH -J DFT
##SBATCH --ntasks=16
##SBATCH --mem-per-cpu=3900
##SBATCH --array=0-1

echo "============================================================"
echo "Job ID : $SLURM_JOB_ID"
echo "Job Name : $SLURM_JOB_NAME"
echo "Starting on : $(date)"
echo "Running on node : $SLURMD_NODENAME"
echo "Current directory : $(pwd)"
echo "============================================================"

#conda activate rdmc_env

#xtb
source /home/gridsan/groups/RMG/Software/xtb-6.4.1/share/xtb/config_env.bash
XTB_PATH=/home/gridsan/groups/RMG/Software/xtb-6.4.1/bin
export PATH=$XTB_PATH:$PATH

#RDMC for gaussian-xtb
RDMC_PATH=/home/gridsan/groups/RMG/Software/RDMC-main
export PATH=$RDMC_PATH:$PATH
export PYTHONPATH=$PYTHONPATH:$RDMC_PATH

#gaussian
export PATH=$PATH:/home/gridsan/groups/RMG/Software/gaussian/g16
export PATH=$PATH:/home/gridsan/groups/RMG/Software/gaussian/gv
export g16root=/home/gridsan/groups/RMG/Software/gaussian
source /home/gridsan/groups/RMG/Software/gaussian/g16/bsd/g16.profile
export GAUSS_SCRDIR=""

#QMD
QMD_PATH=/home/gridsan/groups/RMG/Software/QM_descriptors_calculation-radical_workflow
input_smiles=reactants_products_wb97xd_and_xtb_opted_ts_combo_results_hashed_chart_aug11b.csv
xyz_semiempirical_opt_dict=reactants_products_aug11b_semiempirical_opted_results_xyz.pkl

scratch_dir=$TMPDIR/$USER/$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID-$LLSUB_RANK-$LLSUB_SIZE
mkdir -p $scratch_dir
echo $scratch_dir

#semiempirical
python $QMD_PATH/main_DFT_opt_freq.py --input_smiles $input_smiles --G16_path $g16root/g16 --xyz_semiempirical_opt_dict $xyz_semiempirical_opt_dict --scratch_dir $scratch_dir

rm -rf $scratch_dir
