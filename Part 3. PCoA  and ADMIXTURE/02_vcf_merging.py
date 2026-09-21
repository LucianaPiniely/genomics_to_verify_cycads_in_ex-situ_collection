#same chromosome name so position conflicted between loci
#so we used python script
python3 << 'EOF'
import glob
import gzip

# Read  sample list from namefile.txt
samples = [s.strip() for s in open('namefile.txt') if s.strip()]

# Open output merged VCF file and write header
# The header defines the VCF format version and lists all samples
with open('merged_gardenonly.vcf', 'w') as out:

# write VCF format version
out.write('##fileformat=VCFv4.1\n')

# write column header with all samples
# FORMAT column contains GT (genotype) information
out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' +
             '\t'.join(samples) + '\n')


for vcf in vcf_files:# Loop through each locus VCF file, extract locus name from filename to use as chromosome name
 locus = vcf.replace('vcf_files/', '').replace('.vcf.gz', '') #extract locus name from filename to use as chromosome name
with gzip.open(vcf, 'rt') as f: #open each gzipped VCF file
 header_samples = [] #store sample names from this locus header, each locus may have different samples present

for line in f: 

if line.startswith('#CHROM'): #Read sample order from this locus header, The #CHROM line lists samples present in this specific locus
# Not all samples may be present in every locus
# samples start at column 9 (after CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT)
 header_samples = line.strip().split('\t')[9:]
    continue

# skip all other header lines
if line.startswith('#'):
        continue

#Process each variant record
# Split the line into VCF fields
    fields = line.strip().split('\t')
    chrom     = fields[0]   # chromosome (will be replaced by locus name)
    pos       = fields[1]   # position of the SNP in the alignment
    id_       = fields[2]   # variant ID (usually .)
    ref       = fields[3]   # reference allele
    alt       = fields[4]   # alternative allele
    qual      = fields[5]   # quality score
    filter_   = fields[6]   # filter status
    info      = fields[7]   # info field
    fmt       = fields[8]   # format field (GT)
    genotypes = fields[9:]  # genotypes for samples in this locus

# Map genotypes to sample list
# initialize all samples as missing
gt_map = {s: '.' for s in samples}

 # fill in actual genotypes for samples present in this locus
 for i, s in enumerate(header_samples):
                    if s in gt_map:
                        gt_map[s] = genotypes[i]

# Write the variant record to merged VCF
# get genotypes in master sample order
new_gts = '\t'.join(gt_map[s] for s in samples)

# write record: locus as CHROM, original POS, standard fields, all genotypes
out.write(f'{locus}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t{new_gts}\n')

EOF

#The next step is to convert merged VCF to PLINK format (plink_format.sh)
