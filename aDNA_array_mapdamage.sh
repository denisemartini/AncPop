#!/bin/bash
#SBATCH --job-name=aDNA_array_mapdamage      # job name (shows up in the queue)
#SBATCH --array=1-130            # array job list
#SBATCH --account=uoo02327     # Project Account
#SBATCH --time=2:00:00         # Walltime (HH:MM:SS)
#SBATCH --cpus-per-task=8      # number of cores per task
#SBATCH --mem-per-cpu=1500      # memory/cpu (in MB)
#SBATCH --ntasks=1              # number of tasks (e.g. MPI)
#SBATCH --partition=large       # specify a partition
#SBATCH --hint=nomultithread    # don't use hyperthreading
#SBATCH --chdir=/nesi/nobackup/uoo02327/denise/AncPop_analysis   # directory where you run the job
#SBATCH --output=arraylogs/%x-%A_%a.out      # %x is job name, %A is job number and %a is array task id
#SBATCH --mail-type=ALL         # Optional: Send email notifications
#SBATCH --mail-user=marde569@student.otago.ac.nz     # Use with --mail-type option

#29-11-2018 # [Denise] This script was adapted to run as an array in NeSI, from:
#17-10-2018, 19-06-2018 # some updates and modifications by Hugh Cross
#13-10-2016 #protocol by Catherine Collins #scripting by Sophia Cameron-Christie

# loading necessary modules and environments:
module load Miniconda3/4.4.10
source activate /nesi/project/uoo02327/programs/miniconda_envs/mapdamage/
module load SAMtools
module load picard
module load R
module load GSL/2.3-gimkl-2017a
picard=/opt/nesi/mahuika/picard/2.1.0/picard.jar

#variables
#path to the reference file
ref=/nesi/nobackup/uoo02327/denise/AncPop_analysis/kaka_composite_contamination_file.fasta
#the name of the target mt sequence in the FASTA file (everything after > and before the first space)
mtchr="Nestor_meridionalis"
# a filelist.txt in the working directory should contain all the file names to process
samp=$(awk "NR==$SLURM_ARRAY_TASK_ID" filelist.txt)

# sequencing platform
rgpl=illumina

##########################
# Splitting the info from the file name into multiple fields, for read group setting
flowcell=`echo $samp | cut -d '-' -f1`
lane=`echo $samp | cut -d '_' -f3`
number=`echo $samp | cut -d '_' -f2`
sampname=`echo $samp | cut -d '_' -f1 | cut -d '-' -f3`
library=`echo $samp | cut -d '_' -f1 | cut -d '-' -f2`
#work that info into some read group fields
rgid="${flowcell}_${lane}_${number}"
rgpu="${flowcell}.${lane}"
rglb="${sampname}_${library}"
rgsm="${sampname}"
##########################

cd $samp\_process

#mapDamage: quantifies DNA damage patterns aDNA NGS sequencing reads
for tmpname in ${samp} ${samp}_collapsed
do
mapDamage -i ${tmpname}_remdup.bam -r $ref --rescale
done

#merge the two bam files after mapDamage
if [ -e results_${samp}_remdup/${samp}_remdup.rescaled.bam ] && [ -e results_${samp}_collapsed_remdup/${samp}_collapsed_remdup.rescaled.bam ]; then
echo "found both rescaled bam files, merging"
samtools merge ${samp}_merged_rescale.bam results_${samp}_remdup/${samp}_remdup.rescaled.bam results_${samp}_collapsed_remdup/${samp}_collapsed_remdup.rescaled.bam
bamForReadGroup=${samp}_merged_rescale.bam
#samtools -f merge ${samp}_merged.bam ${samp}_remdup.bam ${samp}_collapsed_remdup.bam
#IF ONLY THE COLLAPSED FILE EXISTS
elif [ ! -e results_${samp}_remdup/${samp}_remdup.rescaled.bam ] && [ -e results_${samp}_collapsed_remdup/${samp}_collapsed_remdup.rescaled.bam ]; then
echo "found the collapsed rescaled bam files, not merging"
bamForReadGroup=results_${samp}_collapsed_remdup/${samp}_collapsed_remdup.rescaled.bam
#IF ONLY THE PAIRED-END FILE EXISTS
elif [ -e results_${samp}_remdup/${samp}_remdup.rescaled.bam ] && [ ! -e results_${samp}_collapsed_remdup/${samp}_collapsed_remdup.rescaled.bam ]; then
echo "found the paired-end rescaled bam files, not merging"
bamForReadGroup=results_${samp}_collapsed_remdup/${samp}_collapsed_remdup.rescaled.bam
#IF NEITHER FILE EXISTS
elif [ ! -e results_${samp}_remdup/${samp}_remdup.rescaled.bam ] && [ ! -e results_${samp}_collapsed_remdup/${samp}_collapsed_remdup.rescaled.bam ]; then
echo "cannot find the rescaled bams from MapDamage, check MapDamage has run on either the paired-end or collapsed files"
fi

#add the read groups (required for GATK)
# picard tools format changed slightly for new version (Hugh's note)
java -jar $picard AddOrReplaceReadGroups \
I=$bamForReadGroup \
O=${samp}_merged_rescale_rdgrps.bam \
RGID=$rgid \
RGLB=$rglb \
RGPL=$rgpl \
RGPU=$rgpu \
RGSM=$rgsm

samtools index ${samp}_merged_rescale_rdgrps.bam

echo "finished processing $samp"
