#!/bin/bash -l
#SBATCH -J restart
#SBATCH --array=0-49

echo "============================================================"
echo "Job ID : $SLURM_JOB_ID"
echo "Job Name : $SLURM_JOB_NAME"
echo "Starting on : $(date)"
echo "Running on node : $SLURMD_NODENAME"
echo "Current directory : $(pwd)"
echo "============================================================"

SubmitDir=`pwd`

for i in {1..5}; do
    cd $SubmitDir
    for folder in $SubmitDir/output/FF_conf/inputs/inputs_*; do

        folderind="${folder##$SubmitDir/output/FF_conf/inputs/inputs_}"
        echo "folderind $folderind"
        
        for filename in $folder/*.tmp; do
            
            input=`basename $filename .tmp`
            
            echo "input $input"
            if [ -e $folder/$input.tmp ]
            then
                mv $folder/$input.tmp $folder/$input.in
            fi
        done
    done
done


# for i in {1..5}; do
#     cd $SubmitDir
#     for folder in $SubmitDir/output/DLPNO_sp/inputs/inputs_*; do

#         folderind="${folder##$SubmitDir/output/DLPNO_sp/inputs/inputs_}"
#         echo "folderind $folderind"
        
#         for filename in $folder/*.tmp; do
            
#             input=`basename $filename .tmp`
            
#             echo "input $input"
#             if [ -e $folder/$input.tmp ]
#             then
#                 mv $folder/$input.tmp $folder/$input.in
#             fi
#         done
#     done
# done

