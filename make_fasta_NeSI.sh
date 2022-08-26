#!/bin/bash -e
#SBATCH --job-name=aDNA_make_fasta      # job name (shows up in the queue)
#SBATCH --account=uoo02327     # Project Account
#SBATCH --time=00:30:00         # Walltime (HH:MM:SS)
#SBATCH --cpus-per-task=4      # number of cores per task
#SBATCH --mem-per-cpu=1500      # memory/cpu (in MB)
#SBATCH --ntasks=1              # number of tasks (e.g. MPI)
#SBATCH --partition=large       # specify a partition
#SBATCH --hint=nomultithread    # don't use hyperthreading
#SBATCH --chdir=/nesi/nobackup/uoo02327/denise/AncPop_analysis   # directory where you run the job
#SBATCH --output=arraylogs/%x-%A_%a.out      # %x is job name, %A is job number and %a is array task id
#SBATCH --mail-type=ALL         # Optional: Send email notifications
#SBATCH --mail-user=marde569@student.otago.ac.nz     # Use with --mail-type option

# 03.12.18 Denise, adapted for NeSI from the below (from Catherine and Sophia?)
#call FASTA files from BAM files produced by fastq_pipeline_loop[version].sh
#Uses GATK's HaploCaller to call best VARIANTS at any depth
#Nullifies any bases (ref or var) at or less than depth specified by the 'mincoverage' variable

# modules to load in NeSI
module load GATK/3.8-1
gatk=/opt/nesi/mahuika/GATK/3.8-1/GenomeAnalysisTK.jar
module load SAMtools
module load picard
picard=/opt/nesi/mahuika/picard/2.1.0/picard.jar
module load snpEff
SnpSift=/opt/nesi/mahuika/snpEff/4.2/SnpSift.jar
module load BCFtools
module load R


######## variables to change ############
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
dna="ancient"
#which BAM file to use
mincoverage="3 5"
#MT
mtchr="Nestor_meridionalis"
#the MINIMUM coverage allowed to call a base for these samples (eg, if set to 1, all bases covered by 1 or more reads will be in FASTA)
# list of samples provided in a separate file
samplist=${projectdir}samplenames.txt
#below: a list of sample names (in any order) as they were processed by the fastq_pipeline_loop[version].sh
#these sample names currently should be in double quotes (around all) and separated by whitespace (spaces or tabs)
#samplist="MS10765 MS10768"


###########################################################################################################################################################################

#variables that shouldn't need changing
ref=$refdir$reffile
bamlist=$grpname.bam.list
#oneunder=$((mincoverage - 1))
gencalldir="call_genotypes/"
grpcalldir="call_group_genotypes2.0/"

###########################################################################################################################################################################

#enter the project directory

if [[ ! -f ${projectdir}filter_cover_for_vcfs.R ]] ; then
    echo "Could not find filter_cover_for_vcfs.R in main processing directory, aborting."
    exit
fi

cd $projectdir

#specify the name of the BAM file to call from
if [ $dna = "ancient" ]; then
echo "You have specified ANCIENT DNA"
bamprefix="_final_merged.bam"
elif [ $dna = "modern" ]; then
echo "You have specified MODERN DNA"
bamprefix="_remdup.bam"
fi

#create dictionary file of the reference if it doesn't exist
dictfile=${reffile%.*}.dict
if [ ! -e "$refdir${reffile%.*}.dict" ]; then
echo "Dictionary file of reference does not exist: creating dictionary"
java -jar $picard CreateSequenceDictionary R=$ref O=$refdir${reffile%.*}.dict
else
echo "Dictionary file found"
fi

#index the reference if it is not indexed already
idxfile=${reffile%.*}.fai
if [ ! -e "$refdir$reffile.fai" ]; then
echo "Reference is not indexed: indexing with SAMTOOLS"
samtools faidx $ref
else
echo "Index for reference file found"
fi

#call variants across entire group
if [ ! -e "${projectdir}${grpcalldir}${grpname}_samples.vcf" ]; then
echo "No VCF of called variants found for this group! Please run variantcall_group_vcf_from_bam script first"
else
echo "Using the called variants from ${grpname}_samples.vcf"
fi

###########################################################################################################################################################################


#begin processing each sample
cat $samplist | while read samp
do
#enter directory
echo "starting FASTA creation for sample $samp"

cd ${samp}_process
mkdir $gencalldir
cd $gencalldir

if [ -e "vcf_list.list" ]; then
echo "Removing an existing list of VCF coverage files"
rm vcf_list.list
fi

#get the VCF for the BEST variant calls
echo "creating a VCF of $samp from the group VCF"

java -jar $gatk \
-T SelectVariants \
-R $ref \
--variant ${projectdir}${grpcalldir}${grpname}_samples.vcf \
-sn $samp \
-env \
-o ${samp}_haploAll.vcf


#make VCF
echo "looking for the BP_RESOLUTION VCF for $samp"

if [ -e "${projectdir}${grpcalldir}${samp}.g.forbam.vcf" ]; then
echo "found BP_RESOLUTION VCF for $samp"
cp ${projectdir}${grpcalldir}${samp}.g.forbam.vcf .
BPvcf=${samp}.g.forbam.vcf
else
echo "Cannot find BP_RESOLUTION VCF for $samp! Creating a VCF of the entire reference for $samp"
java -jar $gatk \
-T HaplotypeCaller \
-R $ref \
-I ../${samp}${bamprefix} \
-ploidy 1 \
-ERC BP_RESOLUTION \
-o $samp.v.vcf
BPvcf=$samp.v.vcf
fi
#

#make a hack vcf from the ERC vcf
echo "selecting low coverage variants from the VCF for $samp"

#get the VCF
vcfheadl=$(awk '/#CHROM/{ print NR; exit }' $BPvcf)

#get coverage
# java -jar $gatk \
# -T DepthOfCoverage \
# -R $ref \
# -o ${samp}_coverage_gatk_table \
# --outputFormat table \
# -I ../../call_group_genotypes2.0/${samp}.gatk.bam

#get the vcf without the header
tail -n +$vcfheadl $BPvcf > $samp.nohead.vcf

#run the VCFs

for cov in $mincoverage
do
#run the script to create a min-cov vcf
echo "running the min.cov script for coverage $cov"
Rscript ${projectdir}filter_cover_for_vcfs.R ${projectdir}coverage/${samp}_coverage.txt $samp.nohead.vcf $cov

#
BPvcf=${samp}.g.forbam.vcf

#recreate the vcf
head -n $vcfheadl $BPvcf > $samp.new$cov.vcf
cat $samp.new$cov.vcf gens.$cov.vcf > $samp.new$cov.full.vcf

#remove the NON_REFs from the vcfs
vcffile=$samp.new$cov.full.vcf
echo "altering $vcffile"
sed -i 's/,<NON_REF>//g' $vcffile
sed -i 's/A\t<NON_REF>/N\tA/g' $vcffile
sed -i 's/C\t<NON_REF>/N\tC/g' $vcffile
sed -i 's/G\t<NON_REF>/N\tG/g' $vcffile
sed -i 's/T\t<NON_REF>/N\tT/g' $vcffile


#make fasta for each coverage
echo "creating a FASTA of bases above $cov reads for $samp"

java -Xmx2g -jar $gatk \
-T FastaAlternateReferenceMaker \
-R $ref \
--variant ${samp}_haploAll.vcf \
-snpmask $samp.new$cov.full.vcf \
-snpmaskPriority \
-L $mtchr \
-o ${samp}.min${cov}.filt.fasta
done

cd ..
cd ..
done
