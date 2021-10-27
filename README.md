# haplotype-mapping
### Estimate founder haplotypes and map phenotypes to them

Run `run.sh` or single commands within it as needed.

See the [R/qtl2 user guide](https://kbroman.org/qtl2/assets/vignettes/user_guide.html) for much more in-depth instructions.

## Additional required inputs

- `founders.vcf.gz`: VCF with founder genotypes. Specify its location in `run.sh`. `run.sh` processes and subsets it and puts the new version in `data/`.

- `genetic_map/MAP4chr{1..20}.txt`: The genetic maps used by qtl2 to convert physical coordinates (bp) to genetic coordinates (centimorgans). Specify its location in `run.sh`.

- `data/obesity_residuals_2020.csv`

- `data/P50_round2_3473_unpruned_sorted.{bed,bim,fam}`

## Notes

- `qtl2_geno.py` checks the REF allele for founders and the cohort and outputs only SNPs where the REFs match. However, there are a lot of mismatches, perhaps because REF vs. ALT was decided arbitrarily. So you probably want to harmonize those somehow so that all valid SNPs are included.

- SNPs are subsetted to those observed in genotype sequencing. The imputed SNPs don't really add information useful for haplotype estimation since they were imputed using estimated haplotypes. And even using only the observed SNPs takes a lot of time and memory. FYI I obtained the observed subset using:

`cat /projects/ps-palmer/apurva/phased_genotypes/chr{1..20}/chr*_hap | grep -v -- '---' | cut -d' ' -f2 | gzip -c > observed_SNPs.txt.gz`

- You can adjust the density of the pseudomarker grid in `run_qtl2.R` to speed up mapping.

- Instead of the basic `calc_genoprob` function, `run_qtl2.R` uses `calc_genoprob_fst` from the `qtl2fst` library, which computes probabilities one chromosome at a time, saving them to disk and returning an index pointing to the data.

- Many cryptic errors while running qtl2 can occur from running out of memory due to large data size.
