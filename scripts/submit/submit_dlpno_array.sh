#!/bin/bash -l
#SBATCH -J orca
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000
#SBATCH --array=0-31
#SBATCH -o slurm-orca_%a.out

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

#openmpi
ompi=/home/gridsan/groups/RMG/Software/ompi-3.1.4
PATH=$ompi/bin:$PATH
LD_LIBRARY_PATH=$ompi/lib:$ompi/etc:$LD_LIBRARY_PATH

#Orca
orcadir=/home/gridsan/groups/RMG/Software/orca_4_2_1_linux_x86-64_openmpi314
export PATH=$PATH:$orcadir
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$orcadir
echo "orcaversion"
which orca
echo $PATH

#QMD
QMD_PATH=/home/gridsan/groups/RMG/Software/QM_descriptors_calculation-radical_workflow
export PYTHONPATH=$QMD_PATH:$PYTHONPATH

#inputs
input_smiles=inputs/reactants_products_aug11b_inputs.csv
output_name=reactants_products_aug11b
xyz_DFT_opt_dict=${output_name}_dft_opted_results_xyz.pkl
DLPNO_sp_n_procs=24
DLPNO_sp_job_ram=4000

#svp calculations

DLPNO_folder='DLPNO_sp'
DLPNO_level_of_theory='uHF dlpno-ccsd(t) def2-svp def2-svp/c TightSCF NormalPNO'
echo $DLPNO_folder
echo $DLPNO_level_of_theory

python -u $QMD_PATH/scripts/calculation/main_make_dlpno_sp_input.py \
--input_smiles $input_smiles \
--xyz_DFT_opt_dict $xyz_DFT_opt_dict \
--DLPNO_folder $DLPNO_folder \
--DLPNO_level_of_theory $DLPNO_level_of_theory \
--DLPNO_sp_n_procs $DLPNO_sp_n_procs \
--DLPNO_sp_job_ram $DLPNO_sp_job_ram \
--ORCA_path $orcadir \
--task_id $SLURM_ARRAY_TASK_ID \
--num_tasks $SLURM_ARRAY_TASK_COUNT

SubmitDir=`pwd`

echo "Starting DLPNO calculations..."

for i in {1..5}; do
    cd $SubmitDir
    for folder in $SubmitDir/output/$DLPNO_folder/inputs/inputs_*; do

        folderind="${folder##$SubmitDir/output/$DLPNO_folder/inputs/inputs_}"
        echo "folderind $folderind"
        
        for filename in $folder/*.in; do
            
            input=`basename $filename .in`
            
            echo "input $input"
            if [ -e $folder/$input.in ]
            then
                mv $folder/$input.in $folder/$input.tmp
                ScratchDir=$HOME/scratch/$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID-$input
                echo "ScratchDir $ScratchDir"
                mkdir -p $ScratchDir

                cd $ScratchDir
                cp $folder/$input.tmp $input.in
                $orcadir/orca $input.in > $input.log
                cp $input.log $SubmitDir/output/$DLPNO_folder/outputs/outputs_$folderind/
                if [ -e $input.log ]
                then
                    if grep -Fq "aborting the run" $input.log
                    then
                        if grep -Fq "ORCA finished by" $input.log
                        then
                            echo "done but failed"
                            rm $folder/$input.tmp
                        else
                            echo "failed"
                            mv $folder/$input.tmp $folder/$input.in
                        fi
                    else
                        echo "done"
                        rm $folder/$input.tmp
                    fi
                else
                    echo "failed"
                    mv $folder/$input.tmp $folder/$input.in
                fi
                cd $SubmitDir
                rm -rf $ScratchDir
            fi
        done
    done
done

# f12 calculations
DLPNO_folder='DLPNO_sp_f12'
DLPNO_level_of_theory='uHF UNO DLPNO-CCSD(T)-F12D cc-pvtz-f12 def2/J cc-pvqz/c cc-pvqz-f12-cabs RIJCOSX VeryTightSCF NormalPNO'
echo $DLPNO_folder
echo $DLPNO_level_of_theory

python -u $QMD_PATH/scripts/calculation/main_make_dlpno_sp_input.py \
--input_smiles $input_smiles \
--xyz_DFT_opt_dict $xyz_DFT_opt_dict \
--DLPNO_folder $DLPNO_folder \
--DLPNO_level_of_theory $DLPNO_level_of_theory \
--DLPNO_sp_n_procs $DLPNO_sp_n_procs \
--DLPNO_sp_job_ram $DLPNO_sp_job_ram \
--ORCA_path $orcadir \
--task_id $SLURM_ARRAY_TASK_ID \
--num_tasks $SLURM_ARRAY_TASK_COUNT

SubmitDir=`pwd`

echo "Starting DLPNO calculations..."

for i in {1..5}; do
    cd $SubmitDir
    for folder in $SubmitDir/output/$DLPNO_folder/inputs/inputs_*; do

        folderind="${folder##$SubmitDir/output/$DLPNO_folder/inputs/inputs_}"
        echo "folderind $folderind"
        
        for filename in $folder/*.in; do
            
            input=`basename $filename .in`
            
            echo "input $input"
            if [ -e $folder/$input.in ]
            then
                mv $folder/$input.in $folder/$input.tmp
                ScratchDir=$HOME/scratch/$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID-$input
                echo "ScratchDir $ScratchDir"
                mkdir -p $ScratchDir

                cd $ScratchDir
                cp $folder/$input.tmp $input.in
                $orcadir/orca $input.in > $input.log
                cp $input.log $SubmitDir/output/$DLPNO_folder/outputs/outputs_$folderind/
                if [ -e $input.log ]
                then
                    if grep -Fq "aborting the run" $input.log
                    then
                        if grep -Fq "ORCA finished by" $input.log
                        then
                            echo "done but failed"
                            rm $folder/$input.tmp
                        else
                            echo "failed"
                            mv $folder/$input.tmp $folder/$input.in
                        fi
                    else
                        echo "done"
                        rm $folder/$input.tmp
                    fi
                else
                    echo "failed"
                    mv $folder/$input.tmp $folder/$input.in
                fi
                cd $SubmitDir
                rm -rf $ScratchDir
            fi
        done
    done
done

