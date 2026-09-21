# we used those filtered alignments from Paragone

#Packages to install
#IQTREE
#Astral 
#MAFFT
#TrimAl

#we retrimmed using trimal with strictplus parameter to remove 
#trimming (TrimAl)
for i in *.fasta; do
trimal -in $i -out ${i//.fasta/.trimmed.fasta} -strictplus
done    

#realign (MAFFT)
for i in *.trimmed.fasta; do
mafft --auto ${i} > ${i%.*}_realigned.fasta;
done 
    

#Create gene tree (IQTREE)
for i in *.realigned.fasta; do
iqtree -s $i -nt AUTO -m TEST -B 1000 -nm 5000 -wbtl;
done        

#After we got our tree files we concatenated treefiles
cat *treefile > all.tre

#create a phylogeny (Astral)
java -jar /yourAstraldirectory/.astral/ASTRAL/Astral/astral.5.7.8.jar -i all.tre -t 2 -a YOURtipnames.txt -o Phylogeny.tre 2>astral.log

#Here we created two phylogeny
#One with samples with locality and one with all samples from garden
#Check individuals.csv