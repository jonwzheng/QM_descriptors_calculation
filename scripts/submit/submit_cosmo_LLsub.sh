#!/bin/bash -l

echo "============================================================"
echo "Job ID : $SLURM_JOB_ID"
echo "Job Name : $SLURM_JOB_NAME"
echo "Starting on : $(date)"
echo "Running on node : $SLURMD_NODENAME"
echo "Current directory : $(pwd)"
echo "============================================================"

# conda activate rdmc_env
which python

#RDMC
export PYTHONPATH=/home/gridsan/groups/RMG/Software/RDMC:$PYTHONPATH

#COSMO
TURBODIR=/home/gridsan/groups/RMG/Software/TmoleX19/TURBOMOLE
source $TURBODIR/Config_turbo_env
COSMOTHERM_PATH=/home/gridsan/groups/RMG/Software/COSMOtherm2021
COSMO_DATABASE_PATH=/home/gridsan/groups/RMG/COSMO_database/COSMObase2021

#QMD
QMD_PATH=/home/gridsan/groups/RMG/Software/QM_descriptors_calculation-radical_workflow
export PYTHONPATH=$QMD_PATH:$PYTHONPATH

input_smiles=inputs/reactants_products_sep1a_filtered_inputs.csv
xyz_DFT_opt_dict=reactants_products_sep1a_filtered_dft_opted_results_xyz.pkl

scratch_dir=$TMPDIR/$USER/$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID-$LLSUB_RANK-$LLSUB_SIZE
mkdir -p $scratch_dir
echo $scratch_dir

python -u $QMD_PATH/scripts/calculation/main_COSMO_calc.py --input_smiles $input_smiles --xyz_DFT_opt_dict $xyz_DFT_opt_dict --scratch_dir $scratch_dir --COSMO_input_pure_solvents $QMD_PATH/common_solvent_list_final.csv --COSMOtherm_path $COSMOTHERM_PATH --COSMO_database_path $COSMO_DATABASE_PATH --task_id $LLSUB_RANK --num_tasks $LLSUB_SIZE

rm -rf $scratch_dir




