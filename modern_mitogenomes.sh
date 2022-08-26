#!/bin/bash -e
#SBATCH --job-name=modDNA_mito      # job name (shows up in the queue)
#SBATCH --array=1-2            # array job list
#SBATCH --account=uoo02327     # Project Account
#SBATCH --time=6:00:00         # Walltime (HH:MM:SS)
#SBATCH --cpus-per-task=8      # number of cores per task
#SBATCH --mem-per-cpu=1500      # memory/cpu (in MB)
#SBATCH --ntasks=1              # number of tasks (e.g. MPI)
#SBATCH --partition=large       # specify a partition
#SBATCH --hint=nomultithread    # don't use hyperthreading
#SBATCH --chdir=/nesi/nobackup/uoo02327/denise/AncPop_analysis/modern   # directory where you run the job
#SBATCH --output=%x-%A_%a.out      # %x is job name, %A is job number and %a is array task id
#SBATCH --mail-type=ALL         # Optional: Send email notifications
#SBATCH --mail-user=marde569@student.otago.ac.nz     # Use with --mail-type option

# 18.12.18, Denise
# adapting the modern and ancient DNA pipelines to call mitogenomes for my modern samples.
# Specifically, I will be using the modern PE realignment script for the alignment and the aDNA methods for the vcf calling.

module load BWA
module load SAMtools
module load picard
module load BCFtools
module load BEDTools

# variables ##need to list here all of the files and directories that I am going to use
# example to fill:
ref=/nesi/nobackup/uoo02327/denise/AncPop_analysis/modern/Nestor_meridionalis.fasta
datadir=/nesi/nobackup/uoo02327/denise/AncPop_analysis/modern/
#next two lines are optional, they depend on the sample name, it can include only the extension of the sample files
#but they are used below to get the correct sample name in the readgroup
fq1=_R1_001.fastq
fq2=_R2_001.fastq

#the whole list needs to be in "" and the samples need to be separated by a white space
platform="Illumina"
picard=/opt/nesi/mahuika/picard/2.1.0/picard.jar

#index the reference fasta file
if [ ! -e $ref.amb ]; then
echo "Index file of reference does not exist: creating index with BWA"
bwa index $ref
else
echo "BWA Index file found"
fi

#create dictionary file of the reference if it doesn't exist
dictfile=${ref%.*}.dict
if [ ! -e "${ref%.*}.dict" ]; then
echo "Dictionary file of reference does not exist: creating dictionary"
java -jar $picard CreateSequenceDictionary R=$ref O=${ref%.*}.dict
else
echo "Dictionary file found"
fi

#index the reference if it is not indexed already
idxfile=${ref%.*}.fai
if [ ! -e "$ref.fai" ]; then
echo "Reference is not indexed: indexing with SAMTOOLS"
samtools faidx $ref
else
echo "Index for reference file found"
fi
#####################################################

samp=$(awk "NR==$SLURM_ARRAY_TASK_ID" filelist.txt)


mkdir -p ${samp}_process/

cd ${samp}_process/

echo "processing data for $samp"

#unzip data if necessary

if [ -e $datadir${samp}${fq1}.gz ] || [ -e $datadir${samp}${fq2}.gz ]; then
echo "Unzipping fastq files"
gunzip $datadir${samp}${fq1}.gz
gunzip $datadir${samp}${fq2}.gz
fi

#get readgroup info
##split the machine/lane info from the fq file
infoline=$(head -n 1 ${datadir}$samp$fq1)
instrument=`echo $infoline | cut -d ':' -f1`
instrument=$(echo $instrument | sed -e 's/ /_/g') # this is to fix issues with spaces in the instrument name
instrumentrun=`echo $infoline | cut -d ':' -f2`
flowcell=`echo $infoline | cut -d ':' -f3`
lane=`echo $infoline | cut -d ':' -f4`
index=`echo $infoline | cut -d ':' -f10`
##the next two lines come from the sample name itself (sample and library in my case)
sampname=`echo $samp | cut -d '_' -f1`
library=`echo $samp | cut -d '_' -f2`
#work that info into some
rgid="ID:${instrument}_${instrumentrun}_${flowcell}_${lane}_${index}"
rgpl="PL:${platform}"
rgpu="PU:${flowcell}.${lane}"
rglb="LB:${sampname}_${library}"
rgsm="SM:${sampname}"

#specify the paired-end reads
fastq1=${datadir}$samp$fq1
fastq2=${datadir}$samp$fq2

#transforming fastq to unsortedBAM
echo "producing uBAM for $samp"
java -jar $picard FastqToSam \
FASTQ=$fastq1 \
FASTQ2=$fastq2 \
OUTPUT=${samp}_fastqtosam.bam \
READ_GROUP_NAME=$rgid \
SAMPLE_NAME=$rgsm \
LIBRARY_NAME=$rglb \
PLATFORM_UNIT=$rgpu \
PLATFORM=$rgpl

echo "Finished preprocessing $samp"

echo "sorting $samp"
java -jar $picard SortSam \
I=${samp}_fastqtosam.bam \
O=${samp}_sorted_fqtosam.bam \
SORT_ORDER=queryname \
TMP_DIR=./tmp

#marking IlluminaAdapters
echo "marking IlluminaAdapters for $samp"
java -jar $picard MarkIlluminaAdapters \
I=${samp}_sorted_fqtosam.bam \
O=${samp}_markilluminaadapters.bam \
M=${samp}_markilluminaadapters_metrics.txt \
TMP_DIR=./tmp

#cleaning up
rm ${samp}_fastqtosam.bam

#Aligning and merging with uBAM all piped
echo "Aligning and merging with uBAM for $samp"
java -jar $picard SamToFastq \
I=${samp}_markilluminaadapters.bam \
FASTQ=/dev/stdout QUIET=true \
CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true \
TMP_DIR=./tmp | bwa mem -M -t 8 -p $ref /dev/stdin | java -jar $picard MergeBamAlignment \
ALIGNED_BAM=/dev/stdin \
UNMAPPED_BAM=${samp}_markilluminaadapters.bam \
OUTPUT=${samp}_piped.bam \
R=$ref CREATE_INDEX=true ADD_MATE_CIGAR=true \
CLIP_OVERLAPPING_READS=true \
INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS  \
TMP_DIR=./tmp

#cleaning up
rm ${samp}_markilluminaadapters.bam

#Getting stats for the alignment
echo "Getting stats for $samp"
samtools stats ${samp}_piped.bam > ${samp}_samtools_stats.txt

echo "Finished alignment $samp"

# sorting .bam file and moving it out
echo "Sorting BAM file for $samp"
java -jar $picard SortSam \
I=${samp}_piped.bam \
O=${samp}_sorted.bam \
SORT_ORDER=coordinate \
TMP_DIR=./tmp

# removing duplicates from .bam file and moving it out
echo "Remove duplicates from BAM file for $samp"
java -jar $picard MarkDuplicates \
I=${samp}_sorted.bam \
O=${samp}_sorted_rmdup.bam \
M=${samp}_marked_dup_metrics.txt \
REMOVE_DUPLICATES=TRUE \
ASSUME_SORTED=TRUE \
TMP_DIR=./tmp

mv ${samp}_sorted_rmdup.bam ..

cd ..

# index bam file and check coverage
samtools index ${samp}_sorted_rmdup.bam
samtools depth -r Nestor_meridionalis -a ${samp}_sorted_rmdup.bam > ${samp}_coverage.txt
Rscript ../make_coverage_plots1.2.R ${samp} Nestor_meridionalis
# prepare a mask file
bedtools genomecov -ibam ${samp}_sorted_rmdup.bam -bga | awk '$4 < 5' > ${samp}_mask.bed


# call consensus
bcftools mpileup --min-MQ 20 --min-BQ 20 \
--output-type u \
--fasta-ref $ref ${samp}_sorted_rmdup.bam | bcftools call \
--ploidy 1 \
--multiallelic-caller \
--output-type v - --output ${samp}.vcf

bgzip ${samp}.vcf
tabix -p vcf ${samp}.vcf.gz

bcftools consensus --fasta-ref $ref \
--haplotype 1 --missing "N" --mask ${samp}_mask.bed \
${samp}.vcf.gz > ${samp}.fa
