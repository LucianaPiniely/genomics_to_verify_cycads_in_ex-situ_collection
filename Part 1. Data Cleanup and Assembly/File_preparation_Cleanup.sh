##File preparation
#Most files coming straight from the sequencing facility have a few extra characters in the names.
#Clean up the file names that have something like samplename_S79_R1_001.fastq.gz:
#to remove the S79:

for filename in *fastq.gz; 
do mv "$filename" "$(echo "$filename" | sed 's/_S[0-9]*_/_/')"; done

#to remove the _001_:

for filename in *fastq.gz;  do mv "$filename" "$(echo "$filename" | sed 's/_001.fastq.gz/.fastq.gz/')"; done

#do note that if you are trying to clean up something like .1. in a file name, in sed, "."" means the first character in the file name and that will be what is edited.

###TRIMMOMATIC
#To run trimmomatic on one sample, cd into the folder containing your raw reads and run:
java -jar /software/trimmomatic/0.39/trimmomatic-0.39.jar PE -phred33 \
input_R1.fastq.gz input_R2.fastq.gz \
output_R1_paired.fastq.gz output_R1_unpaired.fastq.gz \
output_R2_paired.fastq.gz output_R2_unpaired.fastq.gz \
ILLUMINACLIP:/software/trimmomatic/0.39/adapters/TruSeq3-PE.fa:2:30:10:2:true \
LEADING:10 TRAILING:10 SLIDINGWINDOW:4:20 MINLEN:40

#To run trimmomatic on all samples in a folder, you can use a for loop to create a shell script that will run trimmomatic on all samples.
for _R1 in *_R1*; do
R2=${_R1/_R1.fastq.gz/_R2.fastq.gz}
R1p=${_R1/\.fastq.gz/_paired.fastq.gz}
R1u=${_R1/\.fastq.gz/_unpaired.fastq.gz}
R2p=${R2/\.fastq.gz/_paired.fastq.gz}
R2u=${R2/\.fastq.gz/_unpaired.fastq.gz}
echo java -jar /software/trimmomatic/0.39/trimmomatic-0.39.jar PE -phred33 $_R1 $R2 $R1p $R1u $R2p $R2u ILLUMINACLIP:/software/trimmomatic/0.39/adapters/TruSeq3-PE.fa:2:30:10:2:true LEADING:10 TRAILING:10 SLIDINGWINDOW:4:20 MINLEN:40 >> trimmomatic.sh
done

#make the shell script executable:
chmod +x trimmomatic.sh

#make a slurm script to run trimmomatic on the cluster-called slurmtrim.sh

nano slurmtrim.sh #view and edit slurm

#A sample batch script:
#!/bin/bash
#SBATCH --account=p32835
#SBATCH --partition=normal
#SBATCH --time=15:00:00
#SBATCH --mem=5G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --job-name=trimmomatic
#SBATCH --output=trimmomatic.log
#SBATCH --error=error.trimmomatic.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lucianapiniely2026@u.northwestern.edu
#module load hybpiper
#or
#module load java
./trimmomatic.sh

#run the script
sbatch slurmtrim.sh