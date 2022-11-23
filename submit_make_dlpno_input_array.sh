#!/bin/bash -l
#SBATCH --array=0-49

echo "============================================================"
echo "Job ID : $SLURM_JOB_ID"
echo "Job Name : $SLURM_JOB_NAME"
echo "Starting on : $(date)"
echo "Running on node : $SLURMD_NODENAME"
echo "Current directory : $(pwd)"
echo "============================================================"

conda activate QM_descriptors

#openmpi
PATH=/home/gridsan/groups/RMG/Software/ompi/bin:$PATH
LD_LIBRARY_PATH=/home/gridsan/groups/RMG/Software/ompi/lib:/home/gridsan/groups/RMG/Software/ompi/etc:$LD_LIBRARY_PATH

#Orca
orcadir=/home/gridsan/groups/RMG/Software/orca
export PATH=$PATH:$orcadir
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$orcadir

#QMD
QMD_PATH=/home/gridsan/groups/RMG/Software/QM_descriptors_calculation-radical_workflow
input_smiles=reactants_products_wb97xd_and_xtb_opted_ts_combo_results_hashed_chart_aug11b.csv
xyz_DFT_opt=reactants_products_aug11b_dft_opted_results_xyz.pkl

python $QMD_PATH/main_make_dlpno_sp_input.py --input_smiles $input_smiles --xyz_DFT_opt $xyz_DFT_opt --DLPNO_sp_n_procs 22 --DLPNO_sp_job_ram 3900 --ORCA_path $orcadir
