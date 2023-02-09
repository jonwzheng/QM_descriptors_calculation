#!/bin/bash -l

echo "============================================================"
echo "Job ID : $SLURM_JOB_ID"
echo "Job Name : $SLURM_JOB_NAME"
echo "Starting on : $(date)"
echo "Running on node : $SLURMD_NODENAME"
echo "Current directory : $(pwd)"
echo "============================================================"

source /state/partition1/llgrid/pkg/anaconda/anaconda3-2022b/etc/profile.d/conda.sh
conda activate rdmc_env
which python

#QMD
QMD_PATH=/home/gridsan/groups/RMG/Software/QM_descriptors_calculation-radical_workflow
export PYTHONPATH=$QMD_PATH:$PYTHONPATH
input_smiles=../../../ts_calculations/calculations/feb9a/inputs/ts_feb9a_inputs.csv

scratch_dir=$TMPDIR/$USER/$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID-$LLSUB_RANK-$LLSUB_SIZE
mkdir -p $scratch_dir
echo $scratch_dir

#r p complex semi opt
python $QMD_PATH/scripts/calculation/ff_opt_r_p_complex.py.py --input_smiles $input_smiles --RDMC_path $RDMC_PATH --scratch_dir $scratch_dir

rm -rf $scratch_dir
