## Calculate and save 3D array of R/qtl2 genotype probabilities
## Based on https://kbroman.org/qtl2/assets/vignettes/user_guide.html

library(qtl2)
library(qtl2fst) # For memory-intensive calc_genoprob with many SNPs

CORES <- 6
GRID_SPACING <- 0.1 # Distance between pseudomarkers used for mapping, in centimorgans.

cross <- read_cross2("qtl2/control.yaml")

cat("Calculating haplotype probabilities...\n")
gmap <- insert_pseudomarkers(cross$gmap, step = GRID_SPACING)
pr <- calc_genoprob_fst(cross, "pr", "qtl2/pr", map = gmap, error_prob = 0.01, quiet = FALSE, cores = CORES)
## Uncomment this and comment out the above line to load a previously calculated haplotype probability array:
# pr <- readRDS("qtl2/pr/pr_fstindex.rds")

cat("Thinning out dense regions of SNPs...\n")
grid <- calc_grid(cross$gmap, step = GRID_SPACING)
pr <- probs_to_grid(pr, grid)  # Subset genotype probabilities to grid
saveRDS(pr, "qtl2/pr/pr_grid.rds")
## Uncomment this and comment out the above lines to reload previously computed probabilities:
# pr <- readRDS("qtl2/pr/pr_grid.rds")

cat("Calculating LOCO kinship matrices...\n")
kinship <- calc_kinship(pr, "loco", quiet = FALSE, cores = CORES)

cat("Running genome scan...\n")
out <- scan1(pr, cross$pheno, kinship, quiet = FALSE, cores = CORES)
saveRDS(out, "qtl2/scan1_out.rds")
## To load previously computed scan1 output:
# out <- readRDS("qtl2/scan1_out.rds")

## These are the coefficients per founder per phenotype. You may not need them.
cat("Calculating genome scan coefficients...\n")
apr <- genoprob_to_alleleprob(pr, quiet = FALSE, cores = CORES)
coef <- list()
for (chrom in names(apr)) {
    print(chrom)
    coef[[chrom]] <- list()
    for (phen in colnames(cross$pheno)) {
        print(phen)
        co <- scan1coef(apr[, chrom], cross$pheno[, phen], kinship[[chrom]], quiet = FALSE, cores = CORES)
        coef[[chrom]][[phen]] <- co
    }
}

## Save results with marker positions (in Mb) for plotting etc.
pmap <- interp_map(gmap, cross$gmap, cross$pmap)
scan <- list(pmap = pmap, out = out, coef = coef)
saveRDS(scan, "qtl2/scan.rds")
