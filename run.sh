FOUNDERS=~/br/data/genotype/founders.vcf.gz
GEN_MAP_DIR=~/br/data/genotype/genetic_map

# Convert genotypes to VCF and subset SNPs
## Uncomment this and include the --keep argument to test on a smaller sample:
# cut -d' ' -f1,2 data/P50_round2_3473_unpruned_sorted.fam | head -200 > samples.txt
plink2 --bfile data/P50_round2_3473_unpruned_sorted \
    --extract data/observed_SNPs.txt \
    --export vcf bgz id-paste=iid \
    --out data/geno
    # --keep samples.txt \
tabix data/geno.vcf.gz

# Remove extra INFO and FORMAT fields, add IDs, subset SNPs, and fix strain order in founder VCF
bcftools annotate $FOUNDERS \
    --set-id +'%CHROM:%POS' \
    --remove INFO,FORMAT |
    bcftools view \
    --samples ACI-N,BN-N,BUF-N,F344-N,M520-N,MR-N,WKY-N,WN-N \
    --include ID=@data/observed_SNPs.txt \
    -O z -o data/founders.vcf.gz
tabix data/founders.vcf.gz

# Prepare genotype files for qtl2
python3 qtl2_geno.py \
    --vcf data/geno.vcf.gz \
    --founder-vcf data/founders.vcf.gz \
    --snps data/observed_SNPs.txt \
    --gmap-dir $GEN_MAP_DIR \
    --geno-out qtl2/geno.csv \
    --founder-out qtl2/founder_geno.csv \
    --pmap-out qtl2/pmap.csv \
    --gmap-out qtl2/gmap.csv

# Prepare phenotypes for qtl2
python3 qtl2_pheno.py

Rscript run_qtl2.R
