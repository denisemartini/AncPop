#!/bin/bash -e
#SBATCH --job-name=aDNA_call_variants      # job name (shows up in the queue)
#SBATCH --account=uoo02327     # Project Account
#SBATCH --time=01:00:00         # Walltime (HH:MM:SS)
#SBATCH --cpus-per-task=6      # number of cores per task
#SBATCH --mem-per-cpu=1500      # memory/cpu (in MB)
#SBATCH --ntasks=1              # number of tasks (e.g. MPI)
#SBATCH --partition=large       # specify a partition
#SBATCH --hint=nomultithread    # don't use hyperthreading
#SBATCH --chdir=/nesi/nobackup/uoo02327/denise/AncPop_analysis   # directory where you run the job
#SBATCH --output=arraylogs/%x-%A_%a.out      # %x is job name, %A is job number and %a is array task id
#SBATCH --mail-type=ALL         # Optional: Send email notifications
#SBATCH --mail-user=marde569@student.otago.ac.nz     # Use with --mail-type option

# 03.12.18 Denise, adapted for NeSI from the below (from Catherine and Sophia?)
#call VCF files from BAM files produced by fastq_pipeline_loop[version].sh
#Uses GATK's HaploCaller to call all variants, combine into a group VCF and

# modules to load in NeSI
module load GATK/3.8-1
gatk=/opt/nesi/mahuika/GATK/3.8-1/GenomeAnalysisTK.jar
module load SAMtools
module load picard
picard=/opt/nesi/mahuika/picard/2.1.0/picard.jar
module load snpEff
SnpSift=/opt/nesi/mahuika/snpEff/4.2/SnpSift.jar


#variables to change (description below each)

reffile=kaka_composite_contamination_file.fasta
#name of the FASTA reference file
refdir=/nesi/nobackup/uoo02327/denise/AncPop_analysis/
#directory that contains the FASTA reference file
grpname=highcov
#name of this group of samples
projectdir=/nesi/nobackup/uoo02327/denise/AncPop_analysis/
#project directory, that contains data processed by the fastq_pipeline_loop[version].sh script
ploidy=1
#the ploidy of the data (eg, 1 for mitochondria, 2 for autosomal)
callingtype="split"
#whether to call the variants as "split" or "grouped". "split" is the method recommended by GATK, but the "grouped" method takes significant amounts of time but make work well for aDNA.
#the "split" method may not be recommended if:
#if sample number is about 20 this method is not recommended
#if sequence target is anything larger than the mitochondria this method may not be recommended
#if the depth of sequencing is high (>100 average, as in modern mtDNA) this method may not be recommended
filtertype="hard"
#how to filter the variants; either "VQSR" or "hard"
#for genome-wide data, the recommendation is Variant Quality Score Recalibration (VQSR) but this is not recommended for when only mitochondrial data is available
#because the number of variants are too few; therefore "hard" filtering is used for these variants
#hard filters SHOULD BE EXAMINED AND ADAPTED according to the dataset; further discussion can be found on GATK's site: software.broadinstitute.org/gatk/guide/article?id=6925
dna="ancient"
#the prefix of the .bam file that the coverage should be called from: usuall "_remdup.bam" for modern DNA and "_merged_rescale_rdgrps.bam" for ancient
MQthreshold="2.0"
#the quality threshold for filtering - should be specific to the group
call="yes"
#either yes or [anything else]. If "yes", variants will be called from sample BAMs. If it is not yes, existing sample VCFs will be used.

#below: a list of sample names (in any order) as they were processed by the fastq_pipeline_loop[version].sh
#these sample names currently should be in double quotes (around all) and separated by whitespace (spaces or tabs)
#samplist="ABC843"
# list of samples provided in a separate file
samplist=${projectdir}samplenames.txt



###########################################################################################################################################################################

#variables that shouldn't need changing
ref=$refdir$reffile
grpcalldir="call_group_genotypes2.0"

###########################################################################################################################################################################

# "13" error in reference file solved with dos2unix
# dos2unix $ref
sed -i 's/-/N/g' $ref
sed -i 's/n/N/g' $ref

#specify the name of the BAM file to call from
if [ $dna = "ancient" ]; then
echo "You have specified ANCIENT DNA"
bamprefix="_final_merged.bam"
elif [ $dna = "modern" ]; then
echo "You have specified MODERN DNA"
bamprefix="_remdup.bam"
fi

if [ ! -d "$grpcalldir" ]; then
echo "creating new directory for group variant calling"
mkdir $grpcalldir
fi

idxfile=${reffile%.*}.fai
if [ ! -e "$refdir$reffile.fai" ]; then
echo "Reference is not indexed: indexing with SAMTOOLS"
samtools faidx $ref
else
echo "Index for reference file found"
fi

dictfile=${reffile%.*}.dict
if [ ! -e "$refdir${reffile%.*}.dict" ]; then
echo "Dictionary file of reference does not exist: creating dictionary"
java -jar $picard CreateSequenceDictionary R=$ref O=$refdir${reffile%.*}.dict
else
echo "Dictionary file found"
fi

cd $grpcalldir

if [ $callingtype = "split" ]; then

#begin the calling on each sample
if [ $call = "yes" ]; then
echo "variants will be called from individual sample BAMs"
cat $samplist | while read samp
do

# #call sample
java -Xmx8G -jar $gatk \
-T HaplotypeCaller \
-R $ref \
-I $projectdir${samp}_process/${samp}${bamprefix} \
-ploidy 1 \
--emitRefConfidence BP_RESOLUTION \
-bamout ${samp}.gatk.bam \
-o $samp.g.forbam.vcf
#
# #

done
else

echo "variants will NOT be called from individual sample BAMs"

fi
#remove the existing list of vcfs

if [ -e "vcf.list" ]; then
echo "removing the existing list of VCFs"
rm vcf.list
fi

#make a list of the gvcfs to compare group genotypes
cat $samplist | while read samp
do
echo $samp.g.forbam.vcf >> vcf.list
done

#joint genotyping

java -jar $gatk \
-T GenotypeGVCFs \
-R $ref \
--variant vcf.list \
-o ${grpname}_unfiltered_variants.vcf

#variant recalibration (not recommended when only using mtDNA vars) /filtering

#filter the SNPs

echo "pulling out SNPs and INDELs"

#select snps
java -jar $gatk \
-T SelectVariants \
-R $ref \
--variant ${grpname}_unfiltered_variants.vcf \
-selectType SNP \
-o raw_snps.vcf

#select the indels
java -jar $gatk \
-T SelectVariants \
-R $ref \
--variant ${grpname}_unfiltered_variants.vcf \
-selectType INDEL \
-o raw_indels.vcf

echo "filtering the SNPs and INDELs"

if [ $filtertype = "hard" ]; then

#filter snps
java -jar $gatk \
-T VariantFiltration \
-R $ref \
--variant raw_snps.vcf \
--filterExpression "QD < 2.0 || FS > 60.0 || MQ < $MQthreshold" \
--filterName "gatk_snp_filter" \
-o ${grpname}_filtered_snps.vcf


#filter indels
java -jar $gatk \
-T VariantFiltration \
-R $ref \
-V raw_indels.vcf \
--filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
--filterName "gatk_indel_filter" \
-o ${grpname}_filtered_indels.vcf

#recombine indels + snps

#recombine indels + snps

java -jar $gatk \
-T CombineVariants \
-R $ref \
--variant:snp ${grpname}_filtered_snps.vcf \
--variant:ind ${grpname}_filtered_indels.vcf \
-genotypeMergeOptions PRIORITIZE \
-priority snp,ind \
-o ${grpname}_samples_include-filtered.vcf

#remove variants that have failed filtering

java -jar $gatk \
-T SelectVariants \
-R $ref \
--variant ${grpname}_samples_include-filtered.vcf \
-ef \
-o ${grpname}_samples.vcf


else

echo "currently not set up to do VQSR filtering, please perform manually"

fi

########################################################################################

else

#create bamlist
if [ ! -e "$projectdir$grpname.bam.list" ]; then
echo "Could not find bamlist! Creating one from the samplist"
ls ${projectdir}*_process/*${bamprefix} > $projectdir$grpname.TEST.bam.list
else
echo "Found bamlist for group $grpname"
fi


#call the variants
if [ ! -e "$projectdir${grpname}_samples.vcf" ]; then
echo "No VCF of called variants for this group! Calling a new VCF with HaplotypeCaller. This make take some time"
java -jar $gatk \
-T HaplotypeCaller \
-R $ref \
-I ${projectdir}${bamlist} \
-ploidy $ploidy \
-o $projectdir${grpcalldir}${grpname}_samples.vcf
else
echo "VCF of called variants exists for this group. Please delete the existing VCF if you wish to overwrite"
fi
fi

cd ..
