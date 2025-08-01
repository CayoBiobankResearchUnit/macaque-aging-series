library(tidyverse)
library(dplyr)
library(stringr)
library(ashr)
library(mashr)

#for cell type
args <- commandArgs(trailingOnly = TRUE)
#ct = args[1]
#load(sprintf('/scratch/wyang114/outputs/nebula/by_region_NBGMM_mash_inputs/cpc_p05/%s_mash_input.RData',ct))
ct='olig-only'
load('/scratch/klwatki5/R_analysis/bbi/mashr/input/sex/autosomes/oligodendrocytes_mash_sex_auto_DEG_input_no_OPC.Rdata')
coef_mat = auto_coef_mat
se_mat = auto_se_mat

print(ct)
dim(coef_mat)
dim(se_mat)
coef_mat[1:3,1:3]
se_mat[1:3,1:3]

#Read in data
Bhat = coef_mat%>%as.matrix()
Shat = se_mat%>%as.matrix()
data   = mash_set_data(Bhat, Shat)

#Select strong signals
m.1by1 = mash_1by1(data)
strong = get_significant_results(m.1by1,0.05) 

#Obtain initial data-driven covariance matrices
U.pca = cov_pca(data,5,subset=strong)
print(names(U.pca))

#Apply extreme deconvolution
start = Sys.time()
U.ed = cov_ed(data, U.pca, subset=strong)
end = Sys.time()
print(end-start)

#Run mash-fit the model
m.ed = mash(data, U.ed)


#get posterior summary
lfsr = get_lfsr(m.ed)
pm = get_pm(m.ed)
psd = get_psd(m.ed)

print(paste0('strong set size:',length(strong)))
print(paste0('log-likelihood',get_loglik(m.ed),digits = 10))

#save(lfsr,pm,psd,file=sprintf('/scratch/wyang114/outputs/nebula/by_region_NBGMM_mash_outputs/cpc_p05/%s_mash_results.RData',ct))
save(lfsr,pm,psd,file='/home/wyang114/output/nebula/%s_mash_sex_results.RData',ct)
