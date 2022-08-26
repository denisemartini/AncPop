#!/bin/bash -e
#SBATCH --job-name=aDNA_array      # job name (shows up in the queue)
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
module load BWA
module load SAMtools

#OBLIGATORY CHANGES
#change variables
#list sample names
# enter directory for raw data
datadir=/nesi/nobackup/uoo02327/denise/AncPop_analysis/data/
# set directory for trimmed data
trimdir=/nesi/nobackup/uoo02327/denise/AncPop_analysis/trimmed/
# set output directory
outdir=/nesi/nobackup/uoo02327/denise/AncPop_analysis/
#path to the reference file
ref=/nesi/nobackup/uoo02327/denise/AncPop_analysis/kaka_composite_contamination_file.fasta
#the name of the target mt sequence in the FASTA file (everything after > and before the first space)
mtchr="Nestor_meridionalis"

#suffix pattern of the first fastq file
fq1=_R1_001.fastq
#suffix pattern of the second fastq file
fq2=_R2_001.fastq
# a filelist.txt in the working directory should contain all the file names to process
samp=$(awk "NR==$SLURM_ARRAY_TASK_ID" filelist.txt)

##########################
# Run preprocessing
cd $trimdir

#unzip data if necessary
if [ -e $datadir${samp}${fq1}.gz ] || [ -e $datadir${samp}${fq2}.gz ]; then
echo "Unzipping fastq files"
gunzip $datadir${samp}${fq1}.gz
gunzip $datadir${samp}${fq2}.gz
fi

#adaptor removal
AdapterRemoval \
--file1 ${datadir}$samp$fq1 \
--file2 ${datadir}$samp$fq2  \
--collapse  \
--trimns  \
--trimqualities  \
--minlength 25  \
--mm 3  \
--minquality 20  \
--basename $samp

cd $outdir

# rezip the raw files

gzip $datadir${samp}${fq1}
gzip $datadir${samp}${fq2}

##########################

# Prepare reference file
#$bwa index $ref
#index the reference fasta file
if [ ! -e $ref.amb ]; then
echo "Index file of reference does not exist: creating index with BWA"
bwa index $ref
else
echo "BWA Index file found"
fi

##########################

# Running the alignments

#mkdir for sample
mkdir ${samp}_process
cd ${samp}_process

#unzip fq file if necessary
if [ -e $trimdir${samp}.collapsed.gz ]; then
echo "Unzipping fastq file"
gunzip $trimdir${samp}.collapsed.gz
fi

#script changed to use .collapsed AdapterRemoval output file
#Find the SA coordinates of the input reads
bwa aln -n 0.03 -o 2 -l 1024 $ref $trimdir${samp}.pair1.truncated > ${samp}_p1.sai
bwa aln -n 0.03 -o 2 -l 1024 $ref $trimdir${samp}.pair2.truncated > ${samp}_p2.sai
bwa aln -n 0.03 -o 2 -l 1024 $ref $trimdir${samp}.collapsed > ${samp}_collapsed.sai

## Repeat the above line for pair2.truncated and collapsed.truncated
#Generate alignments in the SAM format given PAIRED-END reads

bwa sampe $ref ${samp}_p1.sai ${samp}_p2.sai ${trimdir}${samp}.pair1.truncated ${trimdir}${samp}.pair2.truncated > ${samp}.sam

#Generate alignments in the SAM format given SINGLE-END reads
bwa samse $ref ${samp}_collapsed.sai ${trimdir}${samp}.collapsed > ${samp}_collapsed.sam

#get the total reads for calculating endogenous percent
collapsedtotal=$(samtools view -c ${samp}_collapsed.sam)
pairedtotal=$(samtools view -c ${samp}.sam)

echo "$samp $pairedtotal $collapsedtotal" > total_reads.txt

#You will need to do the following for both merged and unmerged reads

#begin loop
for tmpname in ${samp} ${samp}_collapsed
do

#output the alignments as bam, ignoring alignment with quality score lower than 20
samtools view -b -S -q 20 ${tmpname}.sam > ${tmpname}.bam
#sort the alignments by coordinates
samtools sort ${tmpname}.bam -o ${tmpname}_sort.bam
# try to pipe some commands together

#samtools view -bh -S -q 20 ${tmpname}.sam | samtools sort - -o ${tmpname}_sort.bam

#index the sorted BAM-formatted alignments
samtools index ${tmpname}_sort.bam

#get and print the stats from the indexed BAM file (outputs to a text file)
samtools idxstats ${tmpname}_sort.bam > ${tmpname}_sort_stats.txt

#this removes unmapped reads, they should have been removed with q -20 in the first samtools view step
samtools view -bh -F 0x0004 ${tmpname}_sort.bam $mtchr > ${tmpname}_maponly.bam

# Anna added this to remove duplicates, as samtools will remove duplicates from collapsed and uncollapsed reads in a single step these days
samtools rmdup -S --reference $ref ${tmpname}_maponly.bam ${tmpname}_remdup.bam

#end loop
done
