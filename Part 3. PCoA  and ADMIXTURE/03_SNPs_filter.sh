#we used plink

#Convert VCF to PLINK binary format

plink --vcf mergedgardenonly.vcf --make-bed --out gardenonly --allow-extra-chr --const-fid

#Assign unique SNP IDs
plink --bfile gardenonly --set-missing-var-ids @:# --make-bed --out gardenonly_uniqueid --allow-extra-chr

# Sort the bim file by locus and position 
sort -k1,1 -k4,4n gardenonly_uniqueid.bim > tmp.bim
mv tmp.bim gardenonly_uniqueid.bim

#Rebuild binary files after sorting
plink --bfile gardenonly_uniqueid --make-bed --out gardenonly_full --allow-extra-chr

#LD pruning
plink --bfile gardenonly_full --indep-pairwise 50 10 0.1 --out pruned --allow-extra-chr

#Extract the LD-pruned dataset
plink --bfile gardenonly_full --extract pruned.prune.in --make-bed --out gardenonly_pruned --allow-extra-chr

wc -l < gardenonly_pruned.bim


#I also used created another dataset (dataset 2)that include  samples that  didn't show clear distance on PCoA 
#samples to exclude (samples_exclude.txt) by using --remove on plink

#Finally we ended up with these files gardenonly_pruned.bim, "gardenonly_pruned.bed","gardenonly_pruned.fam" 
#For dataset 2, gardenonly_pruned_subset.bed, gardenonly_pruned_subset.bim, gardenonly_pruned_subset.fam.
#we gona use for PCoA and ADMIXTURE

