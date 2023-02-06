#!/bin/bash -l
#SBATCH -J process
#SBATCH -n 48
#SBATCH -o slurm-process.out

conda activate rdmc_env

#QMD
QMD_PATH=/home/gridsan/groups/RMG/Software/QM_descriptors_calculation-radical_workflow
input_path=inputs/reactants_products_aug11b_inputs.csv 
output_name=reactants_products_aug11b

#FF conf xyzs
python -u $QMD_PATH/scripts/parsing/main_process_FF_conf_result_parallel.py $input_path ${output_name}_ff_opted_results 48

#semiempirical
python -u $QMD_PATH/scripts/parsing/main_process_semi_opt_result_parallel.py $input_path ${output_name}_semiempirical_opted_results 48

#dft
python -u $QMD_PATH/scripts/parsing/main_process_dft_opt_freq_result_parallel.py $input_path ${output_name}_dft_opted_results 48

#dlpno
python -u $QMD_PATH/scripts/parsing/main_process_dlpno_sp_result_parallel.py $input_path ${output_name}_dlpno_sp_results 48

#cosmo
python -u $QMD_PATH/scripts/parsing/main_process_cosmo_result_parallel.py $input_path ${output_name}_cosmo_results 48 $QMD_PATH/common_solvent_list_final.csv 
