#!/bin/bash -l

echo "============================================================"
echo "Job ID : $SLURM_JOB_ID"
echo "Job Name : $SLURM_JOB_NAME"
echo "Starting on : $(date)"
echo "Running on node : $SLURMD_NODENAME"
echo "Current directory : $(pwd)"
echo "============================================================"

source /home/gridsan/groups/RMG/Software/xtb-6.4.1/share/xtb/config_env.bas
conda activate rdmc_env

#QMD
QMD_PATH=/home/gridsan/groups/RMG/Software/QM_descriptors_calculation-radical_workflow
export PYTHONPATH=$QMD_PATH:$PYTHONPATH
input_smiles=/home/gridsan/hwpang/RMG_shared/Projects/Hao-Wei-Oscar-Yunsie/HAbs_calculations/ts_results/roo_ts_checked_feb9a.pkl

scratch_dir=$TMPDIR/$USER/$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID-$LLSUB_RANK-$LLSUB_SIZE
mkdir -p $scratch_dir
echo $scratch_dir

#r p complex semi opt
python $QMD_PATH/scripts/calculation/ff_opt_r_p_complex.py.py --input_smiles $input_smiles --RDMC_path $RDMC_PATH --xyz_DFT_opt_dict $xyz_DFT_opt_dict --scratch_dir $scratch_dir

rm -rf $scratch_dir
