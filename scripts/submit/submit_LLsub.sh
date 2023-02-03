#!/bin/bash -l

# source /etc/profile
# module load anaconda
# conda init
# source /home/gridsan/$USER/.bashrc
# conda activate QM_descriptors

#xtb
source /home/gridsan/groups/RMG/Software/xtb-6.4.1/share/xtb/config_env.bash
XTB_PATH=/home/gridsan/groups/RMG/Software/xtb-6.4.1/bin
export PATH=$XTB_PATH:$PATH
# chmod 750 /home/gridsan/groups/RMG/Software/RDMC-main/rdmc/external/xtb_tools/xtb_gaussian.pl

#RDMC for gaussian-xtb
RDMC_PATH=/home/gridsan/groups/RMG/Software/RDMC-main

#gaussian
export PATH=$PATH:/home/gridsan/groups/RMG/Software/gaussian/g16
export PATH=$PATH:/home/gridsan/groups/RMG/Software/gaussian/gv
export g16root=/home/gridsan/groups/RMG/Software/gaussian
source /home/gridsan/groups/RMG/Software/gaussian/g16/bsd/g16.profile
export GAUSS_SCRDIR=""

#COSMO
TURBODIR=/home/gridsan/groups/RMG/Software/TmoleX19/TURBOMOLE
source $TURBODIR/Config_turbo_env
COSMOTHERMO_PATH=/home/gridsan/groups/RMG/Software/COSMOtherm2021
COSMO_DATABASE_PATH=/home/gridsan/groups/RMG/COSMO_database/COSMObase2021

#openmpi
PATH=/home/gridsan/groups/RMG/Software/ompi/bin:$PATH
LD_LIBRARY_PATH=/home/gridsan/groups/RMG/Software/ompi/lib:/home/gridsan/groups/RMG/Software/ompi/etc:$LD_LIBRARY_PATH

#Orca
orcadir=/home/gridsan/groups/RMG/Software/orca
export PATH=$PATH:$orcadir
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$orcadir

echo "My task ID: " $LLSUB_RANK
echo "Number of Tasks: " $LLSUB_SIZE

#QMD
QMD_PATH=/home/gridsan/groups/RMG/Software/QM_descriptors_calculation-radical_workflow
input_smiles=reactants_products_split_0_wb97xd_and_xtb_opted_ts_combo_results_hashed_lookup_table_sep1a_filtered.csv

python $QMD_PATH/main.py --input_smiles $input_smiles --task_id $LLSUB_RANK --num_tasks $LLSUB_SIZE --XTB_path $XTB_PATH --G16_path $g16root/g16 --ORCA_path $orcadir --RDMC_path $RDMC_PATH --COSMOtherm_path $COSMOTHERMO_PATH --COSMO_database_path $COSMO_DATABASE_PATH --COSMO_input_pure_solvents $QMD_PATH/common_solvent_list_final.csv
