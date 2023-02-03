#!/bin/bash -l
#SBATCH -J orca
#SBATCH --ntasks=22
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3900
#SBATCH --array=0-31
#SBATCH --exclude=c-16-6-4,c-16-11-4,c-16-10-3,c-16-10-2,c-16-7-1,c-16-12-2,d-17-12-2,c-16-12-3,c-16-10-4,c-16-12-3,d-6-7-1,c-16-13-2,c-16-8-4,c-16-6-1,d-6-6-3

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
    cd $SubmitDir
    for folder in $SubmitDir/output/DLPNO_sp/inputs/inputs_*; do

        folderind="${folder##$SubmitDir/output/DLPNO_sp/inputs/inputs_}"
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
