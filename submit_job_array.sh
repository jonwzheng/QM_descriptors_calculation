#!/bin/bash -l
#SBATCH -n 4
#SBATCH --array=0-1

source ~/.bashrc
conda activate QM_descriptors

#nbo
NBO_PATH=~/RMG_shared/Software/NBO7/nbo7/bin
export PATH=$PATH:$NBO_PATH

#xtb
source ~/RMG_shared/Software/xtb-6.4.1/share/xtb/config_env.bash
XTB_PATH=/home/gridsan/groups/RMG/Software/xtb-6.4.1/bin
export PATH=$XTB_PATH:$PATH
chmod 750 ~/Software/RDMC/rdmc/external/xtb_tools/xtb_gaussian.pl

#RDMC for gaussian-xtb
RDMC_PATH=~/Software/RDMC

#gaussian
export PATH=$PATH:/home/gridsan/groups/RMG/Software/gaussian/g16
export PATH=$PATH:/home/gridsan/groups/RMG/Software/gaussian/gv
export g16root=/home/gridsan/groups/RMG/Software/gaussian
source /home/gridsan/groups/RMG/Software/gaussian/g16/bsd/g16.profile

#COSMO
TURBODIR=/home/gridsan/groups/RMG/Software/TmoleX19/TURBOMOLE
source $TURBODIR/Config_turbo_env
COSMOTHERMO_PATH=/home/gridsan/groups/RMG/Software/COSMOtherm2020
COSMO_DATABASE_PATH=/home/gridsan/groups/RMG/COSMO_database

#openmpi
PATH=~/RMG_shared/Software/ompi/bin:$PATH
LD_LIBRARY_PATH=~/RMG_shared/Software/ompi/lib:~/RMG_shared/Software/ompi/etc:$LD_LIBRARY_PATH

#Orca
orcadir=~/RMG_shared/Software/orca
export PATH=$PATH:$orcadir
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$orcadir

echo "My task ID: " $SLURM_ARRAY_TASK_ID
echo "Number of Tasks: " $SLURM_ARRAY_TASK_COUNT

#QMD
QMD_PATH=~/Software/QM_descriptors_calculation
python $QMD_PATH/main.py --input_smiles test.csv --task_id $SLURM_ARRAY_TASK_ID --num_tasks $SLURM_ARRAY_TASK_COUNT --XTB_path $XTB_PATH --G16_path $g16root/g16 --ORCA_path $orcadir --RDMC_path $RDMC_PATH --COSMOtherm_path $COSMOTHERMO_PATH --COSMO_database_path $COSMO_DATABASE_PATH --conf_search_FF GFNFF --nconf 5 --COSMO_solvents h2o methanol --COSMO_solvent_solute_ratios 30 70 0 --COSMO_temperature 25 

