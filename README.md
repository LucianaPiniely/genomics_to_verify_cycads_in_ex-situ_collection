# Using Genomics to verify cycads in Ex situ collection Case study: Tanzanian cycads

Ex situ collections are critical for cycad conservation, but their value depends on correct species identification. 
Older accessions, confiscated material, and hobbyist donations often have incomplete records and outdated nomenclature, and
cycads are hard to identify morphologically because they take years to produce cones and vegetative characters overlap among close relatives.
In this project we used Tanzanian Encephalartos as a model to test whether target-enrichment data can verify species labels and detect mislabeled accessions.
We used Illumina and hybridize our libraries using GoFlag 408 bait set, baits specific for gymnosperms.
This is the workflow of the analyses 

## Part 1. Data Cleanup and Assembly
We removed extra characters in the names and also Illumina adapters. We used trimmomatic to remove the adapter.
For Assembly we used Hybpiper 2
We used GoFlag bait file and created a custom reference file to include Encephalartos only. 
Next, we screened the assemblies for potential paralogs using HybPiper's `paralog_retriever`. This step helped identify loci with multiple competing copies that could represent duplicated genes rather than true orthologs. After reviewing the results, we retained loci with fewer than three paralog copies per gene on average.
To resolve these putative paralogs and infer orthologous relationships, we used ParaGone. This pipeline combines multiple tools for sequence alignment, trimming, error filtering, gene tree construction, and paralog resolution. We ran all three ParaGone orthology-resolution approaches (MI, MO, and RT) and compared the resulting phylogenies to determine which method best reflected the relationships among our samples. 
Once the best-performing approach was identified, alignments were filtered again to retain only loci present in more than 50% of samples. These final alignments were used in all downstream analyses.
Scripts for this part are under Part 1 folder.

## Part 2. Phylogenetic relationship
Using the selected alignments, we reconstructed a phhylogeny. 
We made two phylogenies , one with individuals with locality information and second with all individuals from the gardens. These phylogenies gave us an idea whether some individuals are showing signs of mislabeling
Scripts on Part 2 folder

## Part 3. PCoA and Admixture
We then used the clean alignments to estimate PCoA and admixture
The main thing we noted here is that after SNPs calling, each locus was represented by its own VCF file. The next step was to merge these files into a single VCF. Instead of using bcftools merge which is commonly used , in our case its wasn't suitable for our dataset. Each locus-specific VCF contained the same set of sample names, but different genomic loci. So we customized the Python script (vcf_merging.py). The script maintained a single copy of each sample, extracted variants from every locus-specific VCF, and merged them into a one file while preserving the correct sample identities across loci. 
 So we used this merged VCF provided and do the filtering using then estimates PCoA and ADMIXTURE.
 All scripts for this part are under Part 3 folder
