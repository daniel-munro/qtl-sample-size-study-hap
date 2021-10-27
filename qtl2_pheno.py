import pandas as pd

pheno = pd.read_csv("data/obesity_residuals_2020.csv", index_col=0)
pheno = pheno.drop(columns='idx')
pheno.to_csv('qtl2/pheno.csv', na_rep='NA')

covar = pheno[[]].copy()
covar['generations'] = 90
covar.to_csv('qtl2/covar.csv')
