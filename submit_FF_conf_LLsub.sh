#!/bin/bash -l

echo "============================================================"
echo "Job ID : $SLURM_JOB_ID"
echo "Job Name : $SLURM_JOB_NAME"
echo "Starting on : $(date)"
echo "Running on node : $SLURMD_NODENAME"
echo "Current directory : $(pwd)"
echo "============================================================"

conda activate rdmc_env

#RDMC for gaussian-xtb
RDMC_PATH=/home/gridsan/groups/RMG/Software/RDMC-main
export PATH=$RDMC_PATH:$PATH
export PYTHONPATH=$PYTHONPATH:$RDMC_PATH

#xtb
source /home/gridsan/groups/RMG/Software/xtb-6.4.1/share/xtb/config_env.bash
XTB_PATH=/home/gridsan/groups/RMG/Software/xtb-6.4.1/bin
export PATH=$XTB_PATH:$PATH

#QMD
QMD_PATH=/home/gridsan/groups/RMG/Software/QM_descriptors_calculation-radical_workflow
input_smiles=inputs/reactants_products_aug11b_inputs.csv

scratch_dir=$TMPDIR/$USER/$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID-$LLSUB_RANK-$LLSUB_SIZE
mkdir -p $scratch_dir
echo $scratch_dir

python $QMD_PATH/main_FF_conf.py --input_smiles $input_smiles --XTB_path $XTB_PATH --scratch_dir $scratch_dir

rm -rf $scratch_dir