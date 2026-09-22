# Using Genomics to verify cycads in Ex situ collection Case study: Tanzanian cycads

Ex situ collections are critical for cycad conservation, but their value depends on correct species identification. 
Older accessions, confiscated material, and hobbyist donations often have incomplete records and outdated nomenclature, and
cycads are hard to identify morphologically because they take years to produce cones and vegetative characters overlap among close relatives.
In this project we used Tanzanian Encephalartos as a model to test whether target-enrichment data can verify species labels and detect mislabeled accessions.
We used Illumina and hybridize our libraries using GoFlag 408 bait set, baits specific for gymnosperms.
This is the workflow of the analyses 

## Part 1. Data Cleanup and Assembly
We removed extra characters in the names and also Illumina adapters. We used trimmomatic to remove the adapter
For Assembly we used Hybpiper 2
We used GoFlag bait file and created a custom reference file to include Encephalartos only. 
Because target-enrichment data can contain duplicated gene copies, we investigated paralogs using HybPiper and
then used ParaGone to separate orthologous from paralogous sequences. We compared the orthology-resolution approaches
and chose the best-performing dataset 

## Part 2. Phylogenetic relationship
Using the selected alignments, we reconstructed a phhylogeny. 
We made two phylogenies , one with individuals with locality information and second with all individuals from the gardens

## Assessment of collections identities using PCoA and Admixture
We then used the clean alignments to estimate PCoA and admixture
More details are explained on the scripts
