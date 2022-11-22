#!/bin/bash -l
#SBATCH -J process
#SBATCH -n 48

#conda activate rdmc_env

#QMD
QMD_PATH=/home/gridsan/groups/RMG/Software/QM_descriptors_calculation-radical_workflow
input_path=reactants_products_wb97xd_and_xtb_opted_ts_combo_results_hashed_chart_aug11b.csv 
output_name=reactants_products_aug11b

# #FF conf xyzs
# python $QMD_PATH/main_process_FF_conf_result_parallel.py $input_path ${output_name}_ff_opted_results 48

#semiempirical
python $QMD_PATH/main_process_semi_opt_result_parallel_aug11a.py $input_path ${output_name}_semiempirical_opted_results 48

# ##dft
# python $QMD_PATH/main_process_dft_opt_result_parallel_july27d.py $input_path ${output_name}_dft_opted_results 48

# #dlpno
# python $QMD_PATH/main_process_dlpno_sp_result_parallel.py $input_path ${output_name}_dlpno_sp_results 48

# #cosmo
# python $QMD_PATH/main_process_cosmo_result_parallel.py $input_path ${output_name}_cosmo_results 48
