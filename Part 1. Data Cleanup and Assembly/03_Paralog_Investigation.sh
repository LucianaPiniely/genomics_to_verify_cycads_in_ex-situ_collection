# Investigate paralogs
#retrieve paralogs
hybpiper paralog_retriever -t_dna target.fa  namefile.txt  --fasta_dir_all ./paralogs

#we had so much paralogs highest was 500, so we filetered those with average of 3 and run to paragone to determine orthologs


#paragone pipeline
#we followed the instructions on install and running the pipeline 
https://github.com/chrisjackson-pellicle/ParaGone

nextflow run paragone.nf -c paragone.config -profile standard_singularity --gene_fasta_directory paralogs_all --internal_outgroups 20030281xC,20080218xD,20110155xH --mo --mi --rt
# for my outgroup we used Zamia individuals, 20030281xC,20080218xD,20110155xH


#mi was the best fit for our case so we took its alignments and do more filetering
#we kept only loci with > 50% individuals 
#the remain were the cleaned alignments we used for downstream analyses
