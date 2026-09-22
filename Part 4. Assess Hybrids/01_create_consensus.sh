cd ~
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
conda create -n myEnv
conda activate myEnv
conda install -c bioconda -c conda-forge bcftools=1.9

conda install -c bioconda -c conda-forge bwa

conda install -c bioconda bbmap

#lets generate shell script to generate consensus sequences from hybpiper output

nano generate_consensus_sequences.sh

#!/bin/bash

set +e
set -u
set -o pipefail

CONTIGTYPE="normal"
THREADS=40
SAMPLENAME=""
ALLELE_FREQ=0.15
READ_DEPTH=10
ALLELE_COUNT=4
PIPERDIR="/PATH/trimmed_reads_Encephalartos"
OUTDIR_BASE="../HybPhaser"
NAMELIST="not a file"
CLEANUP="FALSE"
BBMAP="FALSE"
MINID=0.95

while getopts 'ict:s:a:f:d:p:o:n:b:m:' OPTION; do
    case "$OPTION" in

        i) CONTIGTYPE="ISC" ;;
        c) CLEANUP="TRUE" ;;
        t) THREADS=$OPTARG ;;
        s) SAMPLENAME=$OPTARG ;;
        f) ALLELE_FREQ=$OPTARG ;;
        a) ALLELE_COUNT=$OPTARG ;;
        d) READ_DEPTH=$OPTARG ;;
        p) PIPERDIR=$OPTARG ;;
        o) OUTDIR_BASE=$OPTARG ;;
        n) NAMELIST=$OPTARG ;;
        b) BBMAP="TRUE" ;;
        m) MINID=$OPTARG ;;

        ?)
            echo "Usage: generate_consensus_sequences.sh [options]" >&2
            exit 1
            ;;
    esac
done
shift "$(($OPTIND -1))"


# SAMPLE LIST HANDLING


if [[ -f $NAMELIST ]]; then
    SAMPLES=$(<$NAMELIST)

elif [[ $SAMPLENAME != "" ]]; then
    SAMPLES=$SAMPLENAME

else
    echo "No sample name given or namelist file does not exist!"
    exit 1
fi


#main loop
for SAMPLE in $SAMPLES
do
    SECONDS=0

    echo "Collecting read and contig files for $SAMPLE"

    OUTDIR=$OUTDIR_BASE"/01_data"

    if [[ ! -f "$PIPERDIR/$SAMPLE/genes_with_seqs.txt" ]]; then
        echo "!!! Could not find HybPiper output files for $SAMPLE"
        continue
    fi

    # READS — paired FASTQ only

    READ1="$PIPERDIR/${SAMPLE}_R1_paired.fastq.gz"
    READ2="$PIPERDIR/${SAMPLE}_R2_paired.fastq.gz"

    if [[ ! -f "$READ1" || ! -f "$READ2" ]]; then
        echo "Missing paired reads for $SAMPLE"
        continue
    fi

    mkdir -p "$OUTDIR/$SAMPLE/reads"
    cp "$READ1" "$OUTDIR/$SAMPLE/reads/"
    cp "$READ2" "$OUTDIR/$SAMPLE/reads/"

    # CONTIGS 

    if [[ $CONTIGTYPE == "normal" ]]; then
        mkdir -p $OUTDIR/$SAMPLE/contigs
        cp $PIPERDIR/$SAMPLE/*/*/sequences/FNA/*.FNA $OUTDIR/$SAMPLE/contigs/ 2> /dev/null

        for i in $OUTDIR/$SAMPLE/contigs/*.FNA
        do
            FILE=${i/*\/contigs\//}
            GENE=${FILE/.FNA/}
            sed -i "s/ .*//" $i
            sed -i "s/\(>.*\)/\1-$GENE/" $i
        done

        for f in $OUTDIR/$SAMPLE/contigs/*.FNA; do mv "$f" "${f/.FNA/.fasta}"; done

        CONTIGPATH="$OUTDIR/$SAMPLE/contigs"
        mkdir -p "$OUTDIR/$SAMPLE/consensus"

    else
        mkdir -p $OUTDIR/$SAMPLE/intronerated_contigs
        cp $PIPERDIR/$SAMPLE/*/*/sequences/intron/*_supercontig.fasta $OUTDIR/$SAMPLE/intronerated_contigs/ 2> /dev/null

        for f in $OUTDIR/$SAMPLE/intronerated_contigs/*.*
        do
            sed -i "s/ .*//" $f
            awk '{if(NR==1){print $0}else{if($0~/^>/){print "\n"$0}else{printf $0}}}' "$f" > "${f/_supercontig/_intronerated}"
            rm "$f" -f
            echo "" >> "${f/_supercontig/_intronerated}"
        done

        CONTIGPATH="$OUTDIR/$SAMPLE/intronerated_contigs"
        mkdir -p "$OUTDIR/$SAMPLE/intronerated_consensus"
    fi

    mkdir -p "$OUTDIR/$SAMPLE/mapping_files"

    # FUNCTION: MAKE CONSENSUS

    makecons () {
        local CONTIG=$1
        FILE=${CONTIG/*\/*contigs\//}
        GENE=${FILE/.fasta/}

        echo -e '\e[1A\e[K'Generating consensus sequences for $SAMPLE - $GENE

        if [[ $CONTIGTYPE == "normal" ]]; then 
            CONSENSUS=$OUTDIR/$SAMPLE/consensus/$GENE".fasta"
            BAM=$OUTDIR/$SAMPLE/mapping_files/$GENE".bam"
            VCFZ=$OUTDIR/$SAMPLE/mapping_files/$GENE".vcf.gz"
        else
            GENE=${GENE/_intronerated/}
            CONSENSUS=$OUTDIR/$SAMPLE/intronerated_consensus/$GENE"_intronerated.fasta"
            BAM=$OUTDIR/$SAMPLE/mapping_files/$GENE"_intronerated.bam"
            VCFZ=$OUTDIR/$SAMPLE/mapping_files/$GENE"_intronerated.vcf.gz"
        fi

        READS1="$OUTDIR/$SAMPLE/reads/${SAMPLE}_R1_paired.fastq.gz"
        READS2="$OUTDIR/$SAMPLE/reads/${SAMPLE}_R2_paired.fastq.gz"

        bwa index $CONTIG 2> /dev/null
        bwa mem $CONTIG $READS1 $READS2 -t 1 -v 1 2> /dev/null | samtools sort > $BAM

        bcftools mpileup -I -Ov -f $CONTIG $BAM 2> /dev/null \
            | bcftools call -mv -A -Oz -o $VCFZ 2> /dev/null

        bcftools index -f --threads 1 $VCFZ 2> /dev/null 

        bcftools consensus \
            -I \
            -i "(DP4[2]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3]) >= $ALLELE_FREQ && \
                (DP4[0]+DP4[1]+DP4[2]+DP4[3]) >= $READ_DEPTH && \
                (DP4[2]+DP4[3]) >= $ALLELE_COUNT" \
            -f $CONTIG $VCFZ \
            2> /dev/null \
            | awk '{if(NR==1){print $0}else{if($0~/^>/){print "\n"$0}else{printf $0}}}' \
            > $CONSENSUS

        echo "" >> $CONSENSUS

        rm $CONTIG.* 2> /dev/null 
        rm $VCFZ".csi" 2> /dev/null 

        if [[ $CLEANUP == "TRUE" ]]; then 
            rm $VCFZ $BAM
        fi
    }

    # RUN CONSENSUS GENERATION

    max_jobs="$THREADS"
    current_jobs=0

    for CONTIG in $CONTIGPATH/*.fasta; do
        ((current_jobs >= max_jobs)) && wait -n
        makecons "$CONTIG" &
        ((++current_jobs))
    done

    wait

    DURATION_SAMPLE=$SECONDS
    echo -e '\e[1A\e[K'Generated consensus for $SAMPLE in $(($DURATION_SAMPLE / 60)) minutes and $(($DURATION_SAMPLE % 60)) seconds.

done

#activate your env
conda activate myEnv

#run script
bash generate_consensus_sequences.sh -n your.txt -p YOURreadsFOLDER -t 40 -o $HOME/HybPhaser_output

####STEP2
#2A
#SNPs count

#conda install R
#make config.txt and edit as needed per step

```
######################################################
### Configuration File for all HybPhaser R scripts ###
######################################################

# General settings
path_to_output_folder = "yourfolder"
fasta_file_with_targets = "yourtargetfile"
targets_file_format = "DNA"    # "DNA" or "AA" 
path_to_namelist = ".txt "
intronerated_contig = "no"    


###############################
### Part 1:  SNP Assessment ###
###############################

name_for_dataset_optimization_subset = "" 

# missing data
remove_samples_with_less_than_this_propotion_of_loci_recovered = .2
remove_samples_with_less_than_this_propotion_of_target_sequence_length_recovered = .1
remove_loci_with_less_than_this_propotion_of_samples_recovered = .2
remove_loci_with_less_than_this_propotion_of_target_sequence_length_recovered = .1
##lowered filters for target sequence length recovered because we find it less important that our sequences match the targets well, more important that we get good proportion of sequences recovered.
# Paralogs 
remove_loci_for_all_samples_with_more_than_this_mean_proportion_of_SNPs = "outliers"   # any number between 0 and 1, "none" or  "outliers" 
file_with_putative_paralogs_to_remove_for_all_samples = ""
remove_outlier_loci_for_each_sample = "yes"       



##################################
### Part 2:  Clade Association ###
##################################

# set variables to determine thresholds for the dataset optimization and run the script from the main script

path_to_clade_association_folder = "/home/ngavinsmyth/Noras_Sequences/HybPhaser2/04_clade_association"
csv_file_with_clade_reference_names = "/home/ngavinsmyth/Noras_Sequences/HybPhaser2/03_sequence_lists/trimmed_loci_consensus_unique_alleles/4471_ref.csv"
path_to_reference_sequences = "/home/ngavinsmyth/Noras_Sequences/HybPhaser2/03_sequence_lists/samples_consensus"
path_to_read_files_cladeassociation = "/home/ngavinsmyth/Noras_Sequences/HybPhaser2/mapped_reads"
read_type_cladeassociation = "single-end" 
ID_read_pair1 = ""
ID_read_pair2 = ""
file_with_samples_included = ""
path_to_bbmap = "/opt/apps/bbmap/" 
no_of_threads_clade_association = "auto"
run_clade_association_mapping_in_R = "no"   
java_memory_usage_clade_association = ""    # e.g. 2G (for 2GB) needed when java -Xmx error comes up. should be max. 85% of physical memory



#######################
### Part 3: Phasing ###
#######################

# set variables and direct the scripts to the right folders and files before running the single scripts (from inside this file)

path_to_phasing_folder = ""
csv_file_with_phasing_prep_info = ""
path_to_read_files_phasing = ""
read_type_4phasing = "paired-end"   
ID_read_pair1 = "_R1.fastq.gz"    # e.g. "_R1.fastq.gz"
ID_read_pair2 = "_R2.fastq.gz"    # e.g. "_R2.fastq.gz"
reference_sequence_folder = ""
folder_for_phased_reads = ""
folder_for_phasing_stats = ""
path_to_bbmap_executables = ""
no_of_threads_phasing = "auto" 
java_memory_usage_phasing = ""    # e.g. 2G (for 2GB) needed when java -Xmx error comes up. should be max. 85% of physical memory
run_bash_script_in_R = "no"      



################################
### Part 4: Merging Datasets ###
################################

# This script is to generate sequence lists that combine the phased sequences with the normal non-phased ones (but exclude the normal ones of the phased accessions)
# It is also possible to make subsets of samples or loci using lists of included or excluded samples/loci

path_to_sequence_lists_normal = ""
path_to_sequence_lists_phased = ""
path_of_sequence_lists_output = ""

path_to_namelist_normal = ""
path_to_namelist_phased = ""

# to generate subsets one can chose to define samples that are either included or excluded 
file_with_samples_included = ""   
file_with_samples_excluded = ""

# Similarly one can choose to exclude loci or use a list with only the included loci. 
file_with_loci_excluded = ""
file_with_loci_included = ""

exchange_phased_with_not_phased_samples = "yes"   # "yes" or "no" 
include_phased_seqlists_when_non_phased_locus_absent = "no"  # "yes" or "no"

`````#END OF CONFIG.TXT


#below is the R script
#Rscript named 1a_count_snps.R

```
##################
### SNPs count ###
##################

# load config
if (!(exists("config_file"))) {config_file <- "./config.txt"}
source(config_file)

# load packages
library(ape)
library(seqinr)
library(stringr)

# generate folder

output_Robjects <- file.path(path_to_output_folder,"00_R_objects", name_for_dataset_optimization_subset)
dir.create(output_Robjects, showWarnings = F, recursive = T)

#####################################################################################
### Counting polymorphic sites (masked as ambiguity codes in conseneus sequences) ###
#####################################################################################

# getting targets and sample names
if(targets_file_format == "AA"){
  targets <- read.fasta(fasta_file_with_targets, seqtype = "AA", as.string = TRUE, set.attributes = FALSE)
  } else if(targets_file_format == "DNA"){
    targets <- read.fasta(fasta_file_with_targets, seqtype = "DNA", as.string = TRUE, set.attributes = FALSE)
  } else {
    print("Warning! Target file type not set properly. Should be 'DNA' or 'AA'!")
}
targets_name <- unique(gsub(".*-","",labels(targets)))

samples <- readLines(path_to_namelist)



if(intronerated_contig=="yes"){
  intronerated_name <- "intronerated" 
  intronerated_underscore <- "_"
} else {
  intronerated_name <- ""
  intronerated_underscore <- ""
  }

# function for counting ambiguities
seq_stats <- function(file){
  fasta <- read.fasta(file, as.string=TRUE, set.attributes = FALSE )
  seq <- gsub("N|[?]|-","",fasta[[1]])
  c(round(str_length(seq),0),(str_count(seq,"Y|K|R|S|M|y|k|r|s|m|w") + str_count(seq,"W|D|H|B|V|w|d|h|b|v")*2) )
}

# generate empty tables for SNPS and sequence length
tab_snps <- data.frame(loci=targets_name)
tab_length <- data.frame(loci=targets_name)

# fill tables with information on snps and sequence length for each sample and locus 
start_time <- Sys.time()
for(sample in samples){
  #sample <- samples[1]
  print(paste(sample," ", round(difftime(Sys.time(),start_time, units="secs"),0),"s", sep="") )
  
  tab_sample <- data.frame(targets=targets_name, seq_length=NA, ambis=NA, ambi_prop=NA)
  
  
  consensus_files <- list.files(file.path(path_to_output_folder,"01_data",sample, paste(intronerated_name,intronerated_underscore,"consensus", sep="")),full.names = T)
  
  for(consensus_file in consensus_files) {
    
    gene <- gsub("(_intronerated|).fasta","",gsub(".*/","",consensus_file))
    
    file.path(path_to_output_folder,"01_data",sample, "consensus", paste(gene,intronerated_underscore,intronerated_name,".fasta",sep=""))
    
    if(file.info(consensus_file)$size >1){
      stats <- seq_stats(consensus_file)
    } else {
      stats <- c(NA,NA)
    }  
    tab_sample$seq_length[match(gene,tab_sample$targets)] <- stats[1]
    tab_sample$ambis[match(gene,tab_sample$targets)] <- stats[2]  
  
  }
  
  tab_sample$ambi_prop <- tab_sample$ambis/tab_sample$seq_length  
  tab_snps[match(sample,samples)] <- tab_sample$ambi_prop
  tab_length[match(sample,samples)] <- tab_sample$seq_length

}


colnames(tab_snps) <- samples
rownames(tab_snps) <- targets_name
colnames(tab_length) <- samples
rownames(tab_length) <- targets_name

### Generate output tables and save data in Robjects


saveRDS(tab_snps,file=file.path(output_Robjects,"Table_SNPs.Rds"))
saveRDS(tab_length,file=file.path(output_Robjects,"Table_consensus_length.Rds"))
```

#run on terminal
Rscript 1a_count_snps.R

##STEP 2 B
#1b_assess_dataset.R
#####################################################################################
### Generating files for assessment of missing data,. paralogs and heterozygosity ###
#####################################################################################

# load config
if (!(exists("config_file"))) {config_file <- "./config.txt"}
source(config_file)

# load packages
library(ape)

# generate folders
if(name_for_dataset_optimization_subset != ""){
  folder_subset_add <- paste("_",name_for_dataset_optimization_subset, sep="")
}else {
  folder_subset_add <- ""
} 

output_Robjects <- file.path(path_to_output_folder,"00_R_objects", name_for_dataset_optimization_subset)

output_assess <- file.path(path_to_output_folder,paste("02_assessment", folder_subset_add, sep=""))

dir.create(output_assess, showWarnings = F)
dir.create(output_Robjects, showWarnings = F)

# check if SNPs count for subset has been done, if not use SNPs count of first run subset with name "")
if(file.exists(file=file.path(output_Robjects,"Table_SNPs.Rds"))){
  tab_snps <- readRDS(file=file.path(output_Robjects,"Table_SNPs.Rds"))
  tab_length <- readRDS(file=file.path(output_Robjects,"Table_consensus_length.Rds"))
} else {
  tab_snps <- readRDS(file=file.path(path_to_output_folder,"00_R_objects/Table_SNPs.Rds"))
  tab_length <- readRDS(file=file.path(path_to_output_folder,"00_R_objects/Table_consensus_length.Rds"))
}


tab_snps <- as.matrix(tab_snps)
loci <- t(tab_snps)

file.copy(from = config_file, to=file.path(output_assess,"0_used_config_file.txt"))

##########################################################
### Dataset optimization step 1: Reducing missing data ### 
##########################################################


nloci <- length(colnames(loci))
nsamples <- length(rownames(loci))

failed_loci <- which(colSums(is.na(loci))==nrow(loci))
failed_samples <- which(colSums(is.na(tab_snps))==nrow(tab_snps))

# per locus

seq_per_locus <- vector()
for(i in 1:nloci){
  seq_per_locus[i] <- length(which(!(is.na(loci[,i]))))
}
names(seq_per_locus) <- colnames(loci)
seq_per_locus_prop <- seq_per_locus/nsamples


# per sample

seq_per_sample <- vector()
for(i in 1:length(colnames(tab_snps))){
  seq_per_sample[i] <- length(which(!(is.na(tab_snps[,i]))))
}
names(seq_per_sample) <- colnames(tab_snps)
seq_per_sample_prop <- seq_per_sample/nloci


# proportion of target sequence length

if(targets_file_format == "AA"){
  targets_length_all <- lengths(read.FASTA(fasta_file_with_targets, type = "AA"))*3
} else if(targets_file_format == "DNA"){
  targets_length_all <- lengths(read.FASTA(fasta_file_with_targets))
} else { 
  print("Warning! Target file type not set properly. Should be 'DNA' or 'AA'!")
}


gene_names <- unique(gsub(".*-","",gsub(" .*","",names(targets_length_all))))
max_target_length <- vector()

for(i in 1:length(gene_names)){
  max_target_length[i] <- max(targets_length_all[grep(paste("\\b",gene_names[i],"\\b",sep=""),names(targets_length_all))])
}

names(max_target_length) <- gene_names
comb_target_length <- sum(max_target_length)
comb_seq_length_samples <- colSums(tab_length, na.rm = T)
prop_target_length_per_sample <- comb_seq_length_samples/comb_target_length
mean_seq_length_loci <- rowMeans(tab_length, na.rm = T)
names(mean_seq_length_loci) <- gene_names
prop_target_length_per_locus <- mean_seq_length_loci/max_target_length



# application of thresholds

outsamples_missing_loci <- seq_per_sample_prop[which(seq_per_sample_prop < remove_samples_with_less_than_this_propotion_of_loci_recovered)]
outsamples_missing_target <- prop_target_length_per_sample[which(prop_target_length_per_sample < remove_samples_with_less_than_this_propotion_of_target_sequence_length_recovered)]
outsamples_missing <- unique(names(c(outsamples_missing_loci,outsamples_missing_target)))

outloci_missing_samples <- seq_per_locus_prop[which(seq_per_locus_prop < remove_loci_with_less_than_this_propotion_of_samples_recovered)]
outloci_missing_target <- prop_target_length_per_locus[which(prop_target_length_per_locus < remove_loci_with_less_than_this_propotion_of_target_sequence_length_recovered)]
outloci_missing <- unique(names(c(outloci_missing_samples,outloci_missing_target)))


# removing bad loci and samples from the table 

tab_snps_cl1 <- tab_snps

if(length(outsamples_missing) != 0){
  tab_snps_cl1 <- tab_snps_cl1[,-which(colnames(tab_snps) %in% outsamples_missing)]
}

if(length(outloci_missing) != 0){
  tab_snps_cl1 <- tab_snps_cl1[-which(rownames(tab_snps) %in% outloci_missing),]
}

loci_cl1 <- t(tab_snps_cl1)



# output

# graphics

for(i in 1:2){
  
  if(i==1){
    pdf(file=file.path(output_assess,"1_Data_recovered_overview.pdf"), width = 11, height=7)
  } else {
    png(file=file.path(output_assess,"1_Data_recovered_overview.png"), width = 1400, height=1000)
    par(cex.axis=2, cex.lab=2, cex.main=2)
  }
  par(mfrow=c(2,3))
  
  boxplot(seq_per_sample_prop, main=paste("Samples: prop. of",nloci,"loci recovered"), xlab=paste("mean:",round(mean(seq_per_sample_prop, na.rm = TRUE),2), " | median:", round(median(seq_per_sample_prop, na.rm = TRUE),2)," | threshold:",remove_samples_with_less_than_this_propotion_of_loci_recovered," (",length(outsamples_missing_loci)," out)",sep=""))
  abline(h=remove_samples_with_less_than_this_propotion_of_loci_recovered, lty=2, col="red")
  
  boxplot(prop_target_length_per_sample, main=paste("Samples: prop. of target sequence length recovered"), xlab=paste("mean:",round(mean(prop_target_length_per_sample, na.rm = TRUE),2)," | median:", round(median(prop_target_length_per_sample, na.rm = TRUE),2)," | threshold:",remove_samples_with_less_than_this_propotion_of_target_sequence_length_recovered," (",length(outsamples_missing_target)," out)",sep=""))
  abline(h=remove_samples_with_less_than_this_propotion_of_target_sequence_length_recovered, lty=2, col="red")
  
  plot(prop_target_length_per_sample,seq_per_sample_prop, main = "Prop. of loci vs\n prop. of target length", xlab = "Prop. of target length", ylab= "Prop. of loci" )
  abline(h=remove_samples_with_less_than_this_propotion_of_loci_recovered, lty=2, col="red")
  abline(v=remove_samples_with_less_than_this_propotion_of_target_sequence_length_recovered, lty=2, col="red")
  
  
  boxplot(seq_per_locus_prop, main=paste("Loci: prop. of",nsamples,"samples recovered"), xlab=paste("mean:",round(mean(seq_per_locus_prop, na.rm = TRUE),2), " | median:", round(median(seq_per_locus_prop, na.rm = TRUE),2)," | threshold:",remove_loci_with_less_than_this_propotion_of_samples_recovered," (",length(outloci_missing_samples)," out)",sep=""))
  abline(h=remove_loci_with_less_than_this_propotion_of_samples_recovered, lty=2, col="red")
  
  boxplot(prop_target_length_per_locus, main=paste("Loci: prop. of target sequence length recovered"), xlab=paste("mean:",round(mean(prop_target_length_per_locus, na.rm = TRUE),2)," | median:", round(median(prop_target_length_per_locus, na.rm = TRUE),2)," | threshold:",remove_loci_with_less_than_this_propotion_of_target_sequence_length_recovered," (",length(outloci_missing_target)," out)",sep=""))
  abline(h=remove_loci_with_less_than_this_propotion_of_target_sequence_length_recovered, lty=2, col="red")
  
  plot(prop_target_length_per_locus,seq_per_locus_prop, main = "Prop. of samples vs\n prop. of target length", xlab = "Prop. of target length", ylab= "Prop. of samples" )
  abline(h=remove_loci_with_less_than_this_propotion_of_samples_recovered, lty=2, col="red")
  abline(v=remove_loci_with_less_than_this_propotion_of_target_sequence_length_recovered, lty=2, col="red")
  
  dev.off()
}


# tables

tab_seq_per_sample <- cbind(seq_per_sample,round(seq_per_sample_prop,3), round(prop_target_length_per_sample,3))
colnames(tab_seq_per_sample) <- c("No. loci", "Prop. of loci", "Prop. of target length")
write.csv(tab_seq_per_sample, file.path(output_assess, "1_Data_recovered_per_sample.csv"))

tab_seq_per_locus <- cbind(seq_per_locus,round(seq_per_locus_prop,3), round(prop_target_length_per_locus,3))
colnames(tab_seq_per_locus) <- c("No. samples", "Prop.of samples", "Prop. of target length")
write.csv(tab_seq_per_locus, file.path(output_assess, "1_Data_recovered_per_locus.csv"))


# summary text file
summary_file=file.path(output_assess,"1_Summary_missing_data.txt")
cat(file=summary_file, append = FALSE, "Dataset optimisation: Samples and loci removed to reduce missing data\n")

cat(file=summary_file, append = T, "\n", length(failed_samples)," samples failed completely:\n", paste(names(failed_samples)),"\n", sep="")
cat(file=summary_file, append = T, "\n", length(outsamples_missing_loci)," samples are below the threshold (",remove_samples_with_less_than_this_propotion_of_loci_recovered,") for proportion of recovered loci:\n", paste(names(outsamples_missing_loci),"\t",round(outsamples_missing_loci,3),"\n"), sep="")
cat(file=summary_file, append = T, "\n", length(outsamples_missing_target)," samples are below the threshold (",remove_samples_with_less_than_this_propotion_of_target_sequence_length_recovered,") for recovered target sequence length\n", paste(names(outsamples_missing_target),"\t",round(outsamples_missing_target,3),"\n"), sep="")
cat(file=summary_file, append = T, "\nIn total ", length(outsamples_missing), " samples were removed:\n", paste(outsamples_missing,"\n"), sep="")

cat(file=summary_file, append = T, "\n", length(failed_loci)," loci failed completely:\n", paste(names(failed_loci)),"\n", sep="")
cat(file=summary_file, append = T, "\n", length(outloci_missing_samples)," loci are below the threshold (", remove_loci_with_less_than_this_propotion_of_samples_recovered,") for proportion of recovered samples:\n", paste(names(outloci_missing_samples),"\t",round(outloci_missing_samples,3),"\n"), sep="")
cat(file=summary_file, append = T, "\n", length(outloci_missing_target)," loci are below the threshold (",remove_loci_with_less_than_this_propotion_of_target_sequence_length_recovered,") for proportion of recovered target sequence length:\n", paste(names(outloci_missing_target),"\t",round(outloci_missing_target,3),"\n"), sep="")
cat(file=summary_file, append = T, "\nIn total ", length(outloci_missing), " loci were removed:\n", paste(outloci_missing,"\n"), sep="")



############################################################################################
### Dataset optimization step 2, removing paralogs for a) all samples and b) each sample ### 
############################################################################################


### 2a) Paralogs across multiple samples (removing loci with unusually high proportions of SNPs across all samples)
###################################################################################################################

loci_cl1_colmeans <- colMeans(as.matrix(loci_cl1), na.rm = T)
nloci_cl1 <- length(colnames(loci_cl1))
nsamples_cl1 <- length(colnames(tab_snps_cl1))

loci_cl1_colmeans_mean <- round(mean(loci_cl1_colmeans),4)
loci_cl1_colmeans_median <- round(median(loci_cl1_colmeans),4)


# applying chosen threshold
if (length(remove_loci_for_all_samples_with_more_than_this_mean_proportion_of_SNPs) != 1 || remove_loci_for_all_samples_with_more_than_this_mean_proportion_of_SNPs == "none"){
  outloci_para_all <- vector()
  threshold_value <- 1
} else if (remove_loci_for_all_samples_with_more_than_this_mean_proportion_of_SNPs == "outliers"){
  threshold_value <- 1.5*IQR(loci_cl1_colmeans, na.rm = TRUE )+quantile(loci_cl1_colmeans, na.rm = TRUE )[4]
  outloci_para_all_values <- loci_cl1_colmeans[which(loci_cl1_colmeans > threshold_value)]
  outloci_para_all <- names(outloci_para_all_values)
} else if (remove_loci_for_all_samples_with_more_than_this_mean_proportion_of_SNPs == "file"){
  if(file.exists(file_with_putative_paralogs_to_remove_for_all_samples) == FALSE){
    print("File with list of paralogs to remove for all samples does not exist.")
  } else {
    threshold_value <- 1
    outloci_para_all <- readLines(file_with_putative_paralogs_to_remove_for_all_samples)
    outloci_para_all_values <-loci_cl1_colmeans[which(names(loci_cl1_colmeans) %in% outloci_para_all)]
  }
} else {
  threshold_value <- remove_loci_for_all_samples_with_more_than_this_mean_proportion_of_SNPs
  outloci_para_all_values <- loci_cl1_colmeans[which(loci_cl1_colmeans > threshold_value)]
  outloci_para_all <- names(outloci_para_all_values)
}


# color outliers red
colour_outparaall <- rep("black",nloci_cl1)
colour_outparaall[which(colnames(loci_cl1[,order(loci_cl1_colmeans)]) %in% outloci_para_all)] <- "red"
loci_cl1_order_means <- loci_cl1[,order(loci_cl1_colmeans)]


# generate bar graph
for(i in 1:2){
  
  if(i==1){
    pdf(file=file.path(output_assess,"2a_Paralogs_for_all_samples.pdf"), width = 11, height=7)
  } else {
    png(file=file.path(output_assess,"2a_Paralogs_for_all_samples.png"), width = 1400, height=1000)
    par(cex.axis=2, cex.lab=2, cex.main=2)
  }
  
  layout(matrix(c(1,2),2,2, byrow=TRUE), widths=c(5,1))
  
  barplot(sort(loci_cl1_colmeans), col=colour_outparaall, border = NA, las=2,
          main=paste("Mean % SNPs across samples (n=",nsamples_cl1,") for each locus (n=", nloci_cl1,")", sep=""))
  if(length(threshold_value)>0 && remove_loci_for_all_samples_with_more_than_this_mean_proportion_of_SNPs != "file"){abline(h=threshold_value, col="red", lty=2)}
  boxplot(loci_cl1_colmeans, las=2)
  if(length(threshold_value)>0 && remove_loci_for_all_samples_with_more_than_this_mean_proportion_of_SNPs != "file"){abline(h=threshold_value, col="red", lty=2)}
  
  dev.off()
}


# removing marked loci from table
if(length(outloci_para_all)==0) {tab_snps_cl2a <- tab_snps_cl1
} else { tab_snps_cl2a <- tab_snps_cl1[-which(rownames(tab_snps_cl1) %in%  outloci_para_all),]}




### 2b) Paralogs for each sample (removing outlier loci for each sample)
##########################################################################


tab_snps_cl2b <- tab_snps_cl2a

if(!exists("remove_outlier_loci_for_each_sample")) {remove_outlier_loci_for_each_sample <- "no"}

# generate tables without zeros to count only loci with SNPs
tab_snps_cl2a_nozero <- tab_snps_cl2a
tab_snps_cl2a_nozero[which(tab_snps_cl2a_nozero==0)] <- NA
tab_snps_cl2b_nozero <- tab_snps_cl2a_nozero

if(remove_outlier_loci_for_each_sample == "yes" ){

  outloci_para_each <- list()
  outloci_para_each <- sapply(colnames(tab_snps_cl2a),function(x) NULL)
  threshold_para_each <- outloci_para_each
  for(i in 1:length(colnames(tab_snps_cl2a))){
    threshold_i <- 1.5*IQR(tab_snps_cl2a_nozero[,i], na.rm = TRUE ) + quantile(tab_snps_cl2a_nozero[,i] , na.rm = TRUE)[[4]]
    outlier_loci_i <- tab_snps_cl2a_nozero[which(tab_snps_cl2a_nozero[,i] > threshold_i),i]
    outloci_para_each[i] <- list(outlier_loci_i)
    threshold_para_each[i] <- threshold_i
    
    tab_snps_cl2b[which(rownames(tab_snps_cl2a) %in% outlier_loci_i),i] <- NA
    tab_snps_cl2b_nozero[which(rownames(tab_snps_cl2b_nozero) %in% names(outlier_loci_i)),i] <- NA
    
    outliers_color <- "red"
  }
  
} else {
  outloci_para_each <- list()
  outloci_para_each <- sapply(colnames(tab_snps_cl2a),function(x) NULL)
  tab_snps_cl2b <- tab_snps_cl2a
  outliers_color <- "black"
}

tab_length_cl2b <- tab_length[which(rownames(tab_length) %in% rownames(tab_snps_cl2b)),which(colnames(tab_length) %in% colnames(tab_snps_cl2b))]


### output 

# generate graphic
for (i in 1:2){
  if(i==1){
    pdf(file=file.path(output_assess,"2b_Paralogs_for_each_sample.pdf"), width = 10, h=14)
  } else {
    png(file=file.path(output_assess,"2b_Paralogs_for_each_sample.png"), width = 1000, h=1400)
  }
  par(mfrow=c(1,1))
  boxplot(as.data.frame(tab_snps_cl2a_nozero[,order(colMeans(as.matrix(tab_snps_cl2a_nozero), na.rm = T))]), 
          horizontal=T, las=1, yaxt='n',
          main="Proportions of SNPs for all loci per sample\n(only loci with any SNPs)",
          ylab="Samples",
          xlab="Proportion of SNPs",
          col="grey", pars=list(outcol=outliers_color),
          outpch=20
          )
  dev.off()  
}



# write summary text file 
cl2b_file <- file.path(output_assess,"2_Summary_Paralogs.txt")
cat(file=cl2b_file,"Removal of putative paralog loci.")
cat(file=cl2b_file,"Paralogs removed for all samples:\n", append = T)
cat(file=cl2b_file, paste("Variable 'remove_loci_for_all_samples_with_more_than_this_mean_proportion_of_SNPs' set to: ", remove_loci_for_all_samples_with_more_than_this_mean_proportion_of_SNPs,"\n", sep=""), append=T)
if(remove_loci_for_all_samples_with_more_than_this_mean_proportion_of_SNPs=="file"){
  cat(file=cl2b_file, paste("Loci listed in this file were removed: '", file_with_putative_paralogs_to_remove_for_all_samples,"'\n"), append=T)
} else if(remove_loci_for_all_samples_with_more_than_this_mean_proportion_of_SNPs=="none"){
  cat(file=cl2b_file,"None!\n", append = T)
} else {
  cat(file=cl2b_file, paste("Resulting threshold value (mean proportion of SNPs):", round(threshold_value,5),"\n"), append=T)
}
if(length(outloci_para_all)>0){
  cat(file=cl2b_file, paste(length(outloci_para_all)," loci were removed:\n",sep=""), append=T)
  cat(file=cl2b_file, "locus\tmean_prop_SNPs\n", append=T)
  cat(file=cl2b_file, paste(paste(outloci_para_all, round(outloci_para_all_values,4), sep="\t"), collapse="\n"), append=T)
  cat(file=cl2b_file,  "\n\n", append = T)
}

cat(file= file.path(output_assess,"2a_List_of_paralogs_removed_for_all_samples.txt"), paste(c(outloci_para_all,"\n"),collapse = "\n"))



if(remove_outlier_loci_for_each_sample=="yes"){
  cat(file=cl2b_file, "Paralogs removed for each sample:\n", append=T)
  cat(file=cl2b_file, "Sample\tthreshold\t#removed\tnames\n", append=T)
  for(i in 1:length(names(outloci_para_each))){
    cat(file=cl2b_file, names(outloci_para_each)[i],"\t", append=T)  
    cat(file=cl2b_file, round(threshold_para_each[[i]],5), length(outloci_para_each[[i]]), paste(names(outloci_para_each[[i]]),collapse=", "), sep="\t", append=T)  
    cat(file=cl2b_file, "\n", append=T)  
  }
} else{
  cat(file=cl2b_file, "The step for removing paralogs for each samples was skipped.\n", append=T) 
}


# tables

write.csv(tab_snps_cl2b, file = file.path(output_assess,"0_Table_SNPs.csv"))
write.csv(tab_length_cl2b, file = file.path(output_assess,"0_Table_consensus_length.csv"))

# txt file with included samples
write(rownames(loci)[which(!(rownames(loci) %in% outsamples_missing))], file=file.path(output_assess,"0_namelist_included_samples.txt"))

#write(paste(rownames(tab_snps_cl2a),colMeans(t(tab_snps_cl2a), na.rm = T)), file=file.path(output_assess,"Mean_SNPs_loci.txt"))


### save Data as R objects

saveRDS(tab_snps_cl2b,file=file.path(output_Robjects,"Table_SNPs_cleaned.Rds"))
saveRDS(tab_length_cl2b,file=file.path(output_Robjects,"Table_consensus_length_cleaned.Rds"))

saveRDS(outloci_missing,file=file.path(output_Robjects,"outloci_missing.Rds"))
saveRDS(outsamples_missing,file=file.path(output_Robjects,"outsamples_missing.Rds"))
saveRDS(outloci_para_all,file=file.path(output_Robjects,"outloci_para_all.Rds"))
saveRDS(outloci_para_each,file=file.path(output_Robjects,"outloci_para_each.Rds"))



#############################################################################################################
### Generating summary table and graphs for assessment of Locus heterozygosity and allele divergence of samples ###
#############################################################################################################


tab_length <- as.matrix(tab_length)

tab_length_cl2b <- tab_length[which(rownames(tab_length) %in% rownames(tab_snps_cl2b)),which(colnames(tab_length) %in% colnames(tab_snps_cl2b))]

targets_length_cl2b <- sum(max_target_length[which(gsub(".*-","",names(max_target_length)) %in% rownames(tab_snps_cl2b))])

########## generating summary table 

nloci_cl2 <- length(tab_snps_cl2b[,1])
tab_het_ad <- data.frame("sample"=colnames(tab_snps_cl2b))

for(i in 1:length(colnames(tab_snps_cl2b))){
  tab_het_ad$bp[i] <- sum(tab_length_cl2b[,i], na.rm = T)
  tab_het_ad$bpoftarget[i] <- round(sum(tab_length_cl2b[,i], na.rm = T)/targets_length_cl2b,3)*100
  tab_het_ad$paralogs_all[i] <- length(outloci_para_all)
  tab_het_ad$paralogs_each[i] <- length(outloci_para_each[[i]])
  tab_het_ad$nloci[i] <- nloci_cl1-length(outloci_para_all)-length(which(is.na(tab_snps_cl2b[,i])))-length(outloci_para_each[[i]])
  tab_het_ad$allele_divergence[i] <- 100*round(sum(tab_length_cl2b[,i] * tab_snps_cl2b[,i], na.rm = T) / sum(tab_length_cl2b[,i], na.rm = T),5)
  tab_het_ad$locus_heterozygosity[i] <- 100*round(1 - length(which(tab_snps_cl2b[,i]==0))/ (nloci_cl2-length(which(is.na(tab_snps_cl2b[,i])))),4)
  tab_het_ad$'loci with >0.5% SNPs'[i] <- 100*round(1 - length(which(tab_snps_cl2b[,i]<0.005))/ (nloci_cl2-length(which(is.na(tab_snps_cl2b[,i])))),4)
  tab_het_ad$'loci with >1% SNPs'[i] <- 100*round(1 - length(which(tab_snps_cl2b[,i]<0.01))/ (nloci_cl2-length(which(is.na(tab_snps_cl2b[,i])))),4)
  tab_het_ad$'loci with >2% SNPs'[i] <- 100*round(1 - length(which(tab_snps_cl2b[,i]<0.02))/ (nloci_cl2-length(which(is.na(tab_snps_cl2b[,i])))),4)
}

# output as csv file
write.csv(tab_het_ad, file = file.path(output_assess, "4_Summary_table.csv"))

#output as R-object
saveRDS(tab_het_ad, file = file.path(output_Robjects, "Summary_table.Rds"))



### Generating graphs


text_size_mod <- 1
nrows <- length(tab_het_ad[,1])
text_size <- (15+200/nrows)*text_size_mod


for(i in 1:2){
  
  if(i==1){
    pdf(file.path(output_assess,"3_LH_vs_AD.pdf"), h=10,w=10)
  } else {
    png(file.path(output_assess,"3_LH_vs_AD.png"), h=1000,w=1000)
  }
  
  plot(tab_het_ad$allele_divergence,tab_het_ad$locus_heterozygosity,
       xlab="Allele divergence [%]", ylab="Locus heterozygosity [%]", main="Locus heterozygosity vs allele divergence", las=1)
  dev.off()
}


for(i in 1:2){
  
  if(i==1){
    pdf(file.path(output_assess,"3_varLH_vs_AD.pdf"), h=10,w=10)
  } else {
    png(file.path(output_assess,"3_varLH_vs_AD.png"), h=1000,w=1000)
  }
  
  par(mfrow=c(2,2))
  plot(tab_het_ad$allele_divergence,tab_het_ad$locus_heterozygosity,
       xlab="Allele divergence [%]", ylab="Locus heterozygosity (0% SNPs) [%]", main="Locus heterozygosity (0% SNPs) vs allele divergence",las=1
  )
  plot(tab_het_ad$allele_divergence,tab_het_ad$`loci with >0.5% SNPs`,
       xlab="Allele divergence [%]", ylab="Locus heterozygosity (>0.5% SNPs) [%]", main="Locus heterozygosity (>0.5% SNPs) vs allele divergence",las=1
  )
  plot(tab_het_ad$allele_divergence,tab_het_ad$`loci with >1% SNPs`,
       xlab="Allele divergence [%]", ylab="Locus heterozygosity (>1% SNPs) [%]", main="Locus heterozygosity (>1% SNPs) vs allele divergence",las=1
  )
  plot(tab_het_ad$allele_divergence,tab_het_ad$`loci with >2% SNPs`,
       xlab="Allele divergence [%]", ylab="Locus heterozygosity (>2% SNPs) [%]", main="Locus heterozygosity (>2% SNPs) vs allele divergence",las=1
  )
  par(mfrow=c(1,1))
  dev.off()
}
```#END

##2 C
#Generate sequence lists for the cleaned dataset
#r script 1c_generate_sequence_lists.R
####################################
### Generation of sequence lists ###
####################################

# load config
if (!(exists("config_file"))) {config_file <- "./config.txt"}
source(config_file)

# load packages
library(ape)
library(seqinr)
library(stringr)

if(name_for_dataset_optimization_subset != ""){
  folder_subset_add <- paste("_",name_for_dataset_optimization_subset, sep="")
} else {
  folder_subset_add <- ""
} 

output_Robjects <- file.path(path_to_output_folder,"00_R_objects", name_for_dataset_optimization_subset)
output_sequences <- file.path(path_to_output_folder,paste("03_sequence_lists", folder_subset_add, sep=""))

if(intronerated_contig=="yes"){
  intronerated_name <- "intronerated" 
  intronerated_underscore <- "_"
} else {
  intronerated_name <- ""
  intronerated_underscore <- ""
}

targets <- read.fasta(fasta_file_with_targets, as.string=TRUE, set.attributes = FALSE)
targets_name <- unique(gsub(".*-","",labels(targets)))
samples <- readLines(path_to_namelist)


outsamples_missing <- readRDS(file=file.path(output_Robjects,"outsamples_missing.Rds"))
outloci_missing <- readRDS(file=file.path(output_Robjects,"outloci_missing.Rds"))
outloci_para_all <- readRDS(file=file.path(output_Robjects,"outloci_para_all.Rds"))
outloci_para_each <- readRDS(file=file.path(output_Robjects,"outloci_para_each.Rds"))
tab_snps_cl2b <- readRDS(file=file.path(output_Robjects,"Table_SNPs_cleaned.Rds"))

tab_snps <- as.matrix(tab_snps_cl2b)
loci <- t(tab_snps)
failed_loci <- which(colSums(is.na(loci))==nrow(loci))
failed_samples <- which(colSums(is.na(tab_snps))==nrow(tab_snps))


#############################################

folder4seq_consensus_loci <- file.path(output_sequences,"loci_consensus")
folder4seq_consensus_samples <- file.path(output_sequences,"samples_consensus")
folder4seq_contig_loci <- file.path(output_sequences,"loci_contigs")
folder4seq_contig_samples <- file.path(output_sequences,"samples_contigs")


unlink(c(folder4seq_consensus_loci, folder4seq_consensus_samples,  folder4seq_contig_loci, folder4seq_contig_samples),recursive = T) # delete directory, if it existed in order to prevent errors

dir.create(output_sequences, showWarnings = F)
dir.create(folder4seq_consensus_loci, showWarnings = F)
dir.create(folder4seq_consensus_samples, showWarnings = F)
dir.create(folder4seq_contig_loci, showWarnings = F)
dir.create(folder4seq_contig_samples, showWarnings = F)




#########################################################################
### concatenate consensus files across all samples to lists per locus ###
#########################################################################
# this extracts sequences from all subfolders in the HybPiper folder and collates them into one file per locus


# if you are on Linux, then the bash commands cat and sed are used. If not, then R file operations are used, which are slower

if(Sys.info()['sysname']=="Linux"){
  if(intronerated_contig=="yes"){
    for(locus in targets_name){
      command_cat_consensus <- paste("cat",file.path(path_to_output_folder,"01_data/*/intronerated_consensus/",paste(locus,"_intronerated.fasta",sep="")),">",file.path(folder4seq_consensus_loci,paste(locus,"_intronerated_consensus.fasta",sep="")))
      system(command_cat_consensus, ignore.stderr = TRUE)
      command_cat_contig <- paste("cat",file.path(path_to_output_folder,"01_data/*/intronerated_contigs/",paste(locus,"_intronerated.fasta",sep="")),">",file.path(folder4seq_contig_loci,paste(locus,"_intronerated_contig.fasta",sep="")))
      system(command_cat_contig, ignore.stderr = TRUE)
      command_remove_locus_in_seqnames_consensus <- (paste("sed -i 's/-",locus,"//g' ", file.path(folder4seq_consensus_loci,paste(locus,"_intronerated_consensus.fasta",sep="")), sep=""))
      system(command_remove_locus_in_seqnames_consensus) 
      command_remove_locus_in_seqnames_contig <- (paste("sed -i 's/-",locus,"//g' ", file.path(folder4seq_contig_loci,paste(locus,"_intronerated_contig.fasta",sep="")), sep=""))
      system(command_remove_locus_in_seqnames_contig)
    }
  } else {
    for(locus in targets_name){
      command_cat_consensus <- paste("cat",file.path(path_to_output_folder,"01_data/*/consensus/",paste(locus,".fasta",sep="")),">",file.path(folder4seq_consensus_loci,paste(locus,"_consensus.fasta",sep="")))
      command_cat_contig    <- paste("cat",file.path(path_to_output_folder,"01_data/*/contigs/",paste(locus,".fasta",sep="")),">",file.path(folder4seq_contig_loci,paste(locus,"_contig.fasta",sep="")))
      system(command_cat_consensus, ignore.stderr = TRUE)
      system(command_cat_contig, ignore.stderr = TRUE)
      command_remove_locus_in_seqnames_consensus <- (paste("sed -i 's/-",locus,"//g' ", file.path(folder4seq_consensus_loci,paste(locus,"_consensus.fasta",sep="")), sep=""))
      command_remove_locus_in_seqnames_contig <- (paste("sed -i 's/-",locus,"//g' ", file.path(folder4seq_contig_loci,paste(locus,"_contig.fasta",sep="")), sep=""))
      system(command_remove_locus_in_seqnames_consensus)
      system(command_remove_locus_in_seqnames_contig)
    } 
  }  
} else {
  if(intronerated_contig=="yes"){
    for(locus in targets_name){
      #list all fasta files from that locus for all samples
      fasta_files <- list.files(path=file.path(path_to_output_folder,"01_data/"), pattern=paste(locus,"_intronerated.fasta",sep=""), recursive=TRUE, full.names = TRUE)
      #select consensus/contig files
      consensus_files <- grep("consensus",fasta_files,value = TRUE)
      contigs_files <- grep("contigs",fasta_files,value = TRUE)
      #define output files
      output_file_consensus <- file.path(folder4seq_consensus_loci,paste(locus,"_intronerated_consensus.fasta",sep=""))
      output_file_contigs <- file.path(folder4seq_contig_loci,paste(locus,"_intronerated_contigs.fasta",sep=""))
      #generate output files
      file.create(output_file_consensus, overwrite=TRUE)
      file.create(output_file_contigs, overwrite=TRUE)
      #append sample fastas to the empty output file
      file.append(output_file_consensus,consensus_files)
      file.append(output_file_contigs,contigs_files)
      #read lines of each file and remove the "-locus" of the sequence names
      lines_consensus <- gsub(paste("-",locus,sep=""),"",readLines(output_file_consensus))
      lines_contigs <- gsub(paste("-",locus,sep=""),"",readLines(output_file_contigs))
      #write lines into files
      write(lines_consensus, file=output_file_consensus)
      write(lines_contigs, file=output_file_contigs)
    }
  } else {
    for(locus in targets_name){
      #list all fasta files from that locus for all samples
      fasta_files <- list.files(path=file.path(path_to_output_folder,"01_data/"), pattern=paste(locus,".fasta",sep=""), recursive=TRUE, full.names = TRUE)
      #select consensus/contig files
      consensus_files <- grep("consensus",fasta_files,value = TRUE)
      contigs_files <- grep("contigs",fasta_files,value = TRUE)
      #define output files
      output_file_consensus <- file.path(folder4seq_consensus_loci,paste(locus,"_consensus.fasta",sep=""))
      output_file_contigs <- file.path(folder4seq_contig_loci,paste(locus,"_contigs.fasta",sep=""))
      #generate output files
      file.create(output_file_consensus, overwrite=TRUE)
      file.create(output_file_contigs, overwrite=TRUE)
      #append sample fastas to the empty output file
      file.append(output_file_consensus,consensus_files)
      file.append(output_file_contigs,contigs_files)
      #read lines of each file and remove the "-locus" of the sequence names
      lines_consensus <- gsub(paste("-",locus,sep=""),"",readLines(output_file_consensus))
      lines_contigs <- gsub(paste("-",locus,sep=""),"",readLines(output_file_contigs))
      #write lines into files
      write(lines_consensus, file=output_file_consensus)
      write(lines_contigs, file=output_file_contigs)
    } 
  }  
}


# changing interleaved HybPiper files to non-interleaved fasta files
for(file in list.files(folder4seq_contig_loci, full.names = T)){
  file <- list.files(folder4seq_contig_loci, full.names = T)[3]
  lines <- readLines(file)
  conx <- file(file)
  writeLines(str_split(paste(gsub("(>.*)",":\\1:",lines),collapse =""), pattern = ":")[[1]][-1], conx)
  close(conx)
}


# check whether accessions are in the hybpiper folder but not in the samples list.
# if the sample list is smaller some sequences have to be removed from the loci lists

hybpiper_result_dirs <- list.dirs(file.path(path_to_output_folder,"01_data"), full.names = FALSE, recursive = FALSE)
dirs_not_in_sample_list <- hybpiper_result_dirs[which(!(hybpiper_result_dirs %in% samples))]

if(length(dirs_not_in_sample_list) !=0 ){
  # consensus
  for(raw_consensus_file in list.files(folder4seq_consensus_loci, full.names = TRUE)){
    locus_consensus <- readLines(raw_consensus_file)
    lines_with_samplename <- which(gsub(">","",locus_consensus) %in% dirs_not_in_sample_list)
    if(length(lines_with_samplename) !=0){
      lines_to_remove <- c(lines_with_samplename,lines_with_samplename+1)
      locus_file_red <- locus_consensus[-lines_to_remove]
      conn <- file(raw_consensus_file)
      writeLines(locus_file_red, conn)
      close(conn)
    }
  }
  # contig
  for(raw_contig_file in list.files(folder4seq_contig_loci, full.names = TRUE)){
    locus_contig <- readLines(raw_contig_file)
    lines_with_samplename <- which(gsub(">","",locus_contig) %in% dirs_not_in_sample_list)
    if(length(lines_with_samplename) !=0){
      lines_to_remove <- c(lines_with_samplename,lines_with_samplename+1)
      locus_file_hp_red <- locus_contig[-lines_to_remove]
      conn <- file(raw_contig_file)
      writeLines(locus_file_hp_red, conn)
      close(conn)
    }
  }
}


### remove loci from dataset optimization (missing data and paralogs)
######################################################################

## remove loci (failed/missing data/paralogs for all) for all samples

loci_files_consensus <- list.files(path = folder4seq_consensus_loci, full.names = T )
loci_files_contig <- list.files(path = folder4seq_contig_loci, full.names = T )

loci_to_remove <- c(names(failed_loci), outloci_missing, outloci_para_all)

if(length(loci_to_remove)!=0){
	if(intronerated_contig=="no"){
		loci_files_to_remove_consensus <- loci_files_consensus[which(gsub(".*/(.*)_consensus.fasta","\\1",loci_files_consensus) %in% loci_to_remove)]
		loci_files_to_remove_contig <- loci_files_contig[which(gsub(".*/(.*)_contig.fasta","\\1",loci_files_contig) %in% loci_to_remove)]
	} else {
		loci_files_to_remove_consensus <- loci_files_consensus[which(gsub(".*/(.*)_intronerated_consensus.fasta","\\1",loci_files_consensus) %in% loci_to_remove)]
		loci_files_to_remove_contig <- loci_files_contig[which(gsub(".*/(.*)_intronerated_contig.fasta","\\1",loci_files_contig) %in% loci_to_remove)]
	}
	file.remove(loci_files_to_remove_consensus)
	file.remove(loci_files_to_remove_contig)
}


# get vector of all samples that should be removed from every locus

samples_to_remove_4all <- vector()

if(length(failed_samples) > 0 ){
  samples_to_remove_4all <- names(failed_samples)
} 

if(length(outsamples_missing) > 0 ){
  samples_to_remove_4all <-  unique(c(samples_to_remove_4all, outsamples_missing))
} 


## remove samples to be removed from all and sequences in each locus file from paralogs for each sample



for(locus in rownames(tab_snps_cl2b)){
  if(!(locus %in% names(failed_loci))){
  
    if(length(grep(paste("\\b",locus,"\\b",sep=""),outloci_para_each)) >0 ){
      samples_to_remove <- c(samples_to_remove_4all, names(outloci_para_each[grep(paste("\\b",locus,"\\b",sep=""),outloci_para_each)]))
    } else {
      samples_to_remove <- samples_to_remove_4all
    }
    
    
    if(length(samples_to_remove) !=0 ){
      
      # consensus
      locus_consensus <- readLines(file.path(folder4seq_consensus_loci,paste(locus,"_",intronerated_name, intronerated_underscore,"consensus.fasta",sep="")))
      lines_with_samplename <- which(gsub(">","",locus_consensus) %in% samples_to_remove)
      if(length(lines_with_samplename) !=0){
        lines_to_remove <- c(lines_with_samplename,lines_with_samplename+1)
        locus_file_red <- locus_consensus[-lines_to_remove]
        conn <- file(file.path(folder4seq_consensus_loci,paste(locus,"_",intronerated_name, intronerated_underscore,"consensus.fasta",sep="")))
        writeLines(locus_file_red, conn)
        close(conn)
      }
      
      # contig
      locus_contig <- readLines(file.path(folder4seq_contig_loci,paste(locus,"_",intronerated_name, intronerated_underscore,"contig.fasta",sep="")))
      lines_with_samplename <- which(gsub(">","",locus_contig) %in% samples_to_remove)
      if(length(lines_with_samplename) !=0){
        lines_to_remove <- c(lines_with_samplename,lines_with_samplename+1)
        locus_file_hp_red <- locus_contig[-lines_to_remove]
        conn <- file(file.path(folder4seq_contig_loci,paste(locus,"_",intronerated_name, intronerated_underscore,"contig.fasta",sep="")))
        writeLines(locus_file_hp_red, conn)
        close(conn)
      }
    }  
  }
} 


#########################################################################
### concatenate consensus files across all loci to lists per sample   ###
#########################################################################

samples <- readLines(path_to_namelist)

# remove failed samples from list
if(length(failed_samples) != 0){
  samples <- samples[-which(samples %in% names(failed_samples))]
}

# collect all sequences
# if you are on Linux, then the bash commands cat and sed are used. If not, then R file operations are used, which are slower
if(Sys.info()['sysname']=="Linux"){
  if(intronerated_contig=="yes"){
    for(sample in samples){
      command_cat_loci_consensus <- paste("cat",file.path(path_to_output_folder,"01_data/",sample,"/intronerated_consensus/*.fasta"),">",file.path(folder4seq_consensus_samples,paste(sample,"_intronerated_consensus.fasta",sep="")))
      system(command_cat_loci_consensus)
      command_cat_loci_contig <- paste("cat",file.path(path_to_output_folder,"01_data/",sample,"/intronerated_contigs/*.fasta"),">",file.path(folder4seq_contig_samples,paste(sample,"_intronerated_contig.fasta",sep="")))
      system(command_cat_loci_contig)
    } 
  } else {
    for(sample in samples){
      command_cat_loci_consensus <- paste("cat",file.path(path_to_output_folder,"01_data/",sample,"/consensus/*.fasta"),">",file.path(folder4seq_consensus_samples,paste(sample,"_consensus.fasta",sep="")))
      system(command_cat_loci_consensus)
      command_cat_loci_contig <- paste("cat",file.path(path_to_output_folder,"01_data/",sample,"/contigs/*.fasta"),">",file.path(folder4seq_contig_samples,paste(sample,"_contig.fasta",sep="")))
      system(command_cat_loci_contig)
    }
  }  
} else {
  if(intronerated_contig=="yes"){
    for(sample in samples){
      #list all fasta files from the sample for all loci
      consensus_files_samples <- list.files(path=file.path(path_to_output_folder,"01_data/",sample,"/intronerated_consensus/"),pattern="*.fasta", full.names = TRUE)
      contigs_files_samples <- list.files(path=file.path(path_to_output_folder,"01_data/",sample,"/intronerated_contigs/"),pattern="*.fasta", full.names = TRUE)
      #define output files
      output_file_consensus_samples <- file.path(folder4seq_consensus_samples,paste(sample,"_intronerated_consensus.fasta",sep=""))          
      output_file_contigs_samples <- file.path(folder4seq_contig_samples,paste(sample,"_intronerated_contigs.fasta",sep=""))          
      #create output files
      file.create(output_file_consensus_samples, overwrite=TRUE)
      file.create(output_file_contigs_samples, overwrite=TRUE)
      #append fasta files to empty output file
      file.append(output_file_consensus_samples,consensus_files_samples)                
      file.append(output_file_contigs_samples,contigs_files_samples)        
    } 
  } else {
    for(sample in samples){
      
      #list all fasta files from the sample for all loci
      consensus_files_samples <- list.files(path=file.path(path_to_output_folder,"01_data/",sample,"/consensus/"),pattern="*.fasta", full.names = TRUE)
      contigs_files_samples <- list.files(path=file.path(path_to_output_folder,"01_data/",sample,"/contigs/"),pattern="*.fasta", full.names = TRUE)
      #define output files
      output_file_consensus_samples <- file.path(folder4seq_consensus_samples,paste(sample,"_consensus.fasta",sep=""))          
      output_file_contigs_samples <- file.path(folder4seq_contig_samples,paste(sample,"_contigs.fasta",sep=""))          
      #create output files
      file.create(output_file_consensus_samples, overwrite=TRUE)
      file.create(output_file_contigs_samples, overwrite=TRUE)
      #append fasta files to empty output file
      file.append(output_file_consensus_samples,consensus_files_samples)                
      file.append(output_file_contigs_samples,contigs_files_samples)        
    }
  }
}

# changing interleaved HybPiper files to non-interleaved fasta files
for(file in list.files(folder4seq_contig_samples, full.names = T)){
  lines <- readLines(file)
  conx <- file(file)
  writeLines(str_split(paste(gsub("(>.*)",":\\1:",lines),collapse =""), pattern = ":")[[1]][-1], conx)
  close(conx)
}




# remove samples (missing data, paralogs for all)
###################################################

## copy all sequence lists to new folder before removing parts of it


## remove all outlier loci (empty and high allele divergence) (cleaning step 1 and 2)
sample_files_consensus <- list.files(path = folder4seq_consensus_samples, full.names = T )
sample_files_contig <- list.files(path = folder4seq_contig_samples, full.names = T )


if(length(outsamples_missing)==0){
  samples_files_to_remove_consensus=""
  samples_files_to_remove_contig=""
} else {
  samples_files_to_remove_consensus <- sample_files_consensus[which(gsub(".*/(.*)_consensus.fasta","\\1",sample_files_consensus) %in% outsamples_missing)]
  samples_files_to_remove_contig <- sample_files_contig[which(gsub(".*/(.*)_contigs.fasta","\\1",sample_files_contig) %in% outsamples_missing)]

  file.remove(samples_files_to_remove_consensus)
  file.remove(samples_files_to_remove_contig)
}


samples_in <- samples
if(length(outsamples_missing) != 0){
  samples_in <- samples_in[-which(samples %in% outsamples_missing)]
}


loci_to_remove_4all <- vector()

if(length(failed_loci) > 0){
  loci_to_remove_4all <- names(failed_loci)
}

if(length(outloci_missing) > 0){
  loci_to_remove_4all <- c(loci_to_remove_4all, outloci_missing)
}

if(length(outloci_para_all) > 0){
  loci_to_remove_4all <- c(loci_to_remove_4all, outloci_para_all)
}

# remove outlier loci per sample in sample lists
for(sample in samples_in){
  
  
  if(length(which(names(outloci_para_each)%in% sample)) > 0 ){
    loci_to_remove <- c(loci_to_remove_4all, names(outloci_para_each[[which(names(outloci_para_each)%in% sample)]]))
  } else {
    loci_to_remove <- loci_to_remove_4all
  }
  
  
  
  if(length(loci_to_remove)!=0){
    
    
    #consensus
    consensus_file2clean <- file.path(folder4seq_consensus_samples, paste(sample, intronerated_underscore, intronerated_name,"_consensus.fasta",sep=""))
    samples_consensus <- readLines(consensus_file2clean)
    lines_with_lociname <- which(gsub(">.*-","",samples_consensus) %in% loci_to_remove)
    if(length(lines_with_lociname) !=0){
      lines_to_remove <- c(lines_with_lociname,lines_with_lociname+1)
      sample_file_consensus_red <- samples_consensus[-lines_to_remove]
      conn <- file(consensus_file2clean)
      writeLines(sample_file_consensus_red, conn)
      close(conn)
    }
    
    #contig
    contig_file2clean <- file.path(folder4seq_contig_samples, paste(sample,intronerated_underscore, intronerated_name,"_contig.fasta",sep=""))
    samples_contig <- readLines(contig_file2clean)
    lines_with_lociname <- which(gsub(">.*-","",samples_contig) %in% loci_to_remove)
    if(length(lines_with_lociname) !=0){
      lines_to_remove <- c(lines_with_lociname,lines_with_lociname+1)
      sample_file_contig_red <- samples_contig[-lines_to_remove]
      conn <- file(contig_file2clean)
      writeLines(sample_file_contig_red, conn)
      close(conn)
    }
  }
} 
```#END

#run
Rscript 1c_generate_sequence_lists.R