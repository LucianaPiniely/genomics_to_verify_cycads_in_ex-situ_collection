##Installing HybPiper 2
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

#Create a conda environment called "hybpiper" with hybpiper installed 
conda create --name hybpiper -c chrisjackson-pellicle hybpiper

conda activate hybpiper

#Prepare reference file
##GoFlag bait file
#To create this custom reference file create this script using nano command and name it extract_sequences.sh

#!/bin/bash

###Insert the taxon name for your reference of interest between the first two "/   /"

awk '/Encephalartos_barteri/ {print; getline; print}' combinedTarget.fasta > ref_seq.fasta

sed -E "s/>L([0-9]*)_.*_.*_(.*_.*)_[1-2]__REF/>\2-\1/g" ref_seq.fasta > target.fa

rm ref_seq.fasta


##we included the unpaired reads into HybPiper's assembly,
#you need to concatenate the two files together (combine them into one file).
#To concatenate the unpaired files for each sample:
for file in *R1_unpaired*
do
file2=${file//_R1_/_R2_}
file3=${file//R1_unpaired/unpaired}
cat $file $file2 > $file3
done

#after concatenate you need to remove the unpaired reads for R1 and R2
rm *R1_unpaired.fastq.gz *R2_unpaired.fastq.gz

#Now you should only have these three files per sample:
#HE1_R1_paired.fastq.gz, HE1_R2_paired.fastq.gz and HE1_unpaired.fastq.gz


#Assembly
for i in *R1_paired.fastq.gz
do f2=${i//_R1_/_R2_}
f3=${i//R1_paired.fastq.gz/unpaired.fastq.gz}
f4=${i//_R1_paired.fastq.gz/}
echo hybpiper assemble -t_dna target.fa -r $i $f2 --unpaired $f3 --prefix $f4 --bwa --run_intronerate >> assemble.sh
done 

#If you are running Hybpiper 2.1.6 or newer versions, Intronerate is run by default and you have to remove --run_intronerate from the script.

module load Hybpiper

#create a script to run the for loop
#I create a script called slurmAssemble.sh in the current directory with the for loop
#!/bin/bash
#SBATCH --account=youraccount
#SBATCH --partition=normal
#SBATCH --time=24:00:00
#SBATCH --mem=25G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --job-name=assemble
#SBATCH --output=assemble.log
#SBATCH --error=error.assemble.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=youremail

module load hybpiper
./assemble.sh

#make the shell script executable:
chmod +x assemble.sh

#run the script
sbatch slurmAssemble.sh

##SUMMARY STATISTICS
#make a namefile for the assembled data
for i in *R1_paired.fastq.gz
do
name=${i//_R1_paired.fastq.gz/}
echo $name >> namefile.txt
done

#This line generates some tsv's containing statistics about each individual's assemblies
hybpiper stats -t_dna target.fa gene namefile.txt

#This line generates a heatmap based off of the summary statistics in the previous line
hybpiper recovery_heatmap seq_lengths.tsv

hybpiper recovery_heatmap_all seq_lengths.tsv