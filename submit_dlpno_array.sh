#!/bin/bash -l
#SBATCH -J orca
#SBATCH --ntasks=22
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3900
#SBATCH --array=0-31
#SBATCH --exclude=c-16-6-4,c-16-11-4,c-16-10-3,c-16-10-2,c-16-7-1,c-16-12-2,d-17-12-2,c-16-12-3,c-16-10-4,c-16-12-3

echo "============================================================"
echo "Job ID : $SLURM_JOB_ID"
echo "Job Name : $SLURM_JOB_NAME"
echo "Starting on : $(date)"
echo "Running on node : $SLURMD_NODENAME"
echo "Current directory : $(pwd)"
echo "============================================================"

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

SubmitDir=`pwd`

for i in {1..5}; do
    for folder in $SubmitDir/output/DLPNO_sp/inputs/inputs_*; do

        folderind="${folder##$SubmitDir/output/DLPNO_sp/inputs/inputs_}"
        echo "folderind $folderind"
        
        for filename in $folder/*.in; do
            
            input=`basename $filename .in`
            
            echo "input $input"
            if [ -e $folder/$input.in ]
            then
                mv $folder/$input.in $folder/$input.tmp

                ScratchDir=$TMPDIR/$USER/scratch/$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID-$input
                echo "ScratchDir $ScratchDir"
                mkdir -p $ScratchDir

                cd $ScratchDir
                cp $folder/$input.tmp $input.in
                $orcadir/orca $input.in > $input.log
                cp $input.log $SubmitDir/output/DLPNO_sp/outputs/outputs_$folderind/
                if [ -e $input.log ]
                then
                    if grep -Fq "aborting the run" $input.log
                    then
                        echo "failed"
                        mv $folder/$input.tmp $folder/$input.in
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
