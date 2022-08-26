#!/bin/bash
#SBATCH --job-name=aDNA_endogenous      # job name (shows up in the queue)
#SBATCH --account=uoo02327     # Project Account
#SBATCH --time=00:05:00         # Walltime (HH:MM:SS)
#SBATCH --cpus-per-task=1      # number of cores per task
#SBATCH --mem-per-cpu=1500      # memory/cpu (in MB)
#SBATCH --ntasks=1              # number of tasks (e.g. MPI)
#SBATCH --partition=large       # specify a partition
#SBATCH --hint=nomultithread    # don't use hyperthreading
#SBATCH --chdir=/nesi/nobackup/uoo02327/denise/AncPop_analysis   # directory where you run the job
#SBATCH --output=arraylogs/%x-%A_%a.out      # %x is job name, %A is job number and %a is array task id
#SBATCH --mail-type=ALL         # Optional: Send email notifications
#SBATCH --mail-user=marde569@student.otago.ac.nz     # Use with --mail-type option

#30-11-2018 Denise, to compile some stats on the deduplicated and rescaled alignments

# loading necessary modules and environments:
module load SAMtools

grep '^H' all_contamination_stats.txt | awk '{print $1}' > files.txt

echo "File"$'\t'"TotalReads"$'\t'"MappedKaka"$'\t'"%Endogenous"$'\t'"AfterDedup"$'\t'"%EndogAfterDedup" >> endogenous_stats.txt

cat files.txt | while read f
do
  echo -n $f$'\t' >> endogenous_stats.txt
  total=`awk -v var="$f" '$0~var{print $2+$3}' all_contamination_stats.txt`
  mapped=`grep -A23 $f all_contamination_stats.txt | awk '/Nestor/{sum+=$2} END {print sum}'`
  dedup=`samtools view -c ${f}_process/${f}_merged_rescale_rdgrps.bam`
  endog=`awk -v n="$total" -v m="$mapped" 'BEGIN{ print(m/n*100)}'`
  real_end=`awk -v n="$total" -v m="$dedup" 'BEGIN{ print(m/n*100)}'`
  echo $total$'\t'$mapped$'\t'$endog$'\t'$dedup$'\t'$real_end >> endogenous_stats.txt
done
