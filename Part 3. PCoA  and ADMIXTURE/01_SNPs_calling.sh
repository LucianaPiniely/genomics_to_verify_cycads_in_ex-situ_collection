#We need to call SNPs from our clean alignments

#change data from fasta to vcf using snp-sites
https://github.com/sanger-pathogens/snp-sites

#install snp-sites using bioconda
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda
conda install snp-sites

#create a directory to store the vcf files
mkdir -p vcf_files

#convert fasta files to vcf files using snp-sites
for f in *.fasta; do
    base=$(basename "$f" .fasta)
    snp-sites -v -o vcf_files/${base}.vcf "$f"
    echo "Done: $f"
done

# compress and index all VCF files first (required for bcftools merge)
for f in vcf_files/*.vcf; do
    bgzip -f "$f"
    tabix -p vcf "${f}.gz"
    echo "Indexed: $f"
done

#merging all VCF into one 
 #using bcftools could have been easier but the issue is that bcftools merge duplicated sample names across loci instead of merging them properly.
#In our dataset the problem is each VCF file has the same sample names but different loci — bcftools merge treats them as different samples.
#So we used a python script for merging (vcf_merging.py)


