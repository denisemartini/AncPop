#!/bin/bash -e
#SBATCH --job-name=aDNA_array_plot      # job name (shows up in the queue)
#SBATCH --array=1-65            # array job list
#SBATCH --account=uoo02327     # Project Account
#SBATCH --time=00:15:00         # Walltime (HH:MM:SS)
#SBATCH --cpus-per-task=4      # number of cores per task
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
module load SAMtools
module load R

#variables
#path to the reference file
ref=/nesi/nobackup/uoo02327/denise/AncPop_analysis/kaka_composite_contamination_file.fasta
#the name of the target mt sequence in the FASTA file (everything after > and before the first space)
mtchr="Nestor_meridionalis"
# a samplenames.txt in the working directory should contain all the sample names to process
samp=$(awk "NR==$SLURM_ARRAY_TASK_ID" samplenames.txt)
# working directory path
workdir=/nesi/nobackup/uoo02327/denise/AncPop_analysis/
covdir=/nesi/nobackup/uoo02327/denise/AncPop_analysis/coverage/
lengthdir=/nesi/nobackup/uoo02327/denise/AncPop_analysis/length/

#variables for R scripts
covver="make_coverage_plots1.2"
plotreadver="plot_readlengths2.0"
dna="ancient"

###################################
# checking that directories are in place
mkdir -p $covdir
mkdir -p $lengthdir

###################################
# merging all the rescaled and readgrouped bam files for the same sample
mkdir ${samp}_process/

cp H*-${samp}_*_process/H*-${samp}_*_rdgrps.bam ${samp}_process/
cp H*-${samp}_*_process/H*-${samp}_*_sort_stats.txt ${samp}_process/

cd ${samp}_process/

echo "merging $samp"
ls *.bam > bamlist
samtools merge -f ${samp}_final_merged.bam -b bamlist
samtools index ${samp}_final_merged.bam

# creating files for later stats
samtools depth -r $mtchr -a ${samp}_final_merged.bam > $covdir${samp}_coverage.txt
samtools view ${samp}_final_merged.bam | awk '{print length($10)}' | sort -n | uniq -c > $lengthdir${samp}_length_distribution.txt

###################################

echo "making plots and stats for $samp"

ls *stats.txt | cut -d '_' -f 1-3 | sort -u > statslist

echo "Samp PairedTot CollapsedTot" > ${samp}_contamination_stats
cat statslist | while read names
do
  echo -n $names$'\t' >> ${samp}_contamination_stats
  awk '{print $2, $3}' $workdir${names}_process/total_reads.txt >> ${samp}_contamination_stats
	echo "Ref PairedMapped PairedUnmapped" >> ${samp}_contamination_stats
  awk '{print $1,$3,$4}' ${names}_sort_stats.txt >> ${samp}_contamination_stats
  echo "Ref CollapsedMapped CollapsedUnmapped" >> ${samp}_contamination_stats
  awk '{print $1,$3,$4}' ${names}_collapsed_sort_stats.txt >> ${samp}_contamination_stats
done

###################################
#create coverage plot

echo "Creating coverage plot"

cd $covdir

Rscript $workdir${covver}.R $samp $mtchr

#create readlength plot

echo "Creating readlength plot"

cd $lengthdir

Rscript $workdir${plotreadver}.R $samp $PWD

##############################
