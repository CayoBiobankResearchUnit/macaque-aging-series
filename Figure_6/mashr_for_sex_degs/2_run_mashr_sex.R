#!/usr/bin/env Rscript

library(tidyverse)
library(dplyr)
library(stringr)
library(ashr)
library(mashr)

#for cell type
args <- commandArgs(trailingOnly = TRUE)
ct = args[1]

# Autosomes 
load(paste0('/scratch/klwatki5/R_analysis/bbi/mashr/mash_new/input/',ct,'_mash_sex_auto_DEG_input.Rdata'))

# Read in data
Bhat = auto_coef_mat%>%as.matrix()
Shat = auto_se_mat%>%as.matrix()
data   = mash_set_data(Bhat, Shat)

#Select strong signals
m.1by1 = mash_1by1(data)
strong = get_significant_results(m.1by1,0.05)
length(strong)

#Obtain initial data-driven covariance matrices
U.pca = cov_pca(data,5,subset=strong)
print(names(U.pca))

#Apply extreme deconvolution
U.ed = cov_ed(data, U.pca, subset=strong)

#Run mash-fit the model
U.c = cov_canonical(data)  
m.ed = mash(data, U.ed)

# Save intermediate files in case of crash 
save(m.ed,m.1by1,U.pca,U.ed, file=paste0('/scratch/klwatki5/R_analysis/bbi/mashr/mash_new/results/canon_dd/full_',ct,'_mash_sex_auto_0.05.RData'))

#get posterior summary
lfsr = get_lfsr(m.ed)
pm = get_pm(m.ed)
psd = get_psd(m.ed)

save(lfsr,pm,psd,file=sprintf('/scratch/klwatki5/R_analysis/bbi/mashr/mash_new/results/canon_dd/%s_mash_sex_auto_DEG_results_0.05.RData',ct))

################################################################################
# XY code 
################################################################################
# Sex chromosomes
load(paste0('/scratch/klwatki5/R_analysis/bbi/mashr/mash_new/input/',ct,'_mash_sex_xy_DEG_input.RData'))
print(ct)
dim(xy_coef_mat)
dim(xy_se_mat)
xy_coef_mat[1:3,1:3]
xy_se_mat[1:3,1:3]

#Read in data
Bhat = xy_coef_mat%>%as.matrix()
Shat = xy_se_mat%>%as.matrix()
data   = mash_set_data(Bhat, Shat) 

#Select strong signals
m.1by1 = mash_1by1(data)
strong = get_significant_results(m.1by1,0.05)
length(strong)

#Obtain initial data-driven covariance matrices
U.pca = cov_pca(data,5,subset=strong)
print(names(U.pca))

#Apply extreme deconvolution
U.ed = cov_ed(data, U.pca, subset=strong)

#Run mash-fit the model
m.ed = mash(data, U.ed)
save(m.ed,m.1by1,U.pca,U.ed, file=paste0('/scratch/klwatki5/R_analysis/bbi/mashr/mash_new/results/strong_0.05/full_',ct,'_mash_sex_xy.RData'))

#get posterior summary
lfsr = get_lfsr(m.ed)
pm = get_pm(m.ed)
psd = get_psd(m.ed)

save(lfsr,pm,psd,file=sprintf('/scratch/klwatki5/R_analysis/bbi/mashr/mash_new/results/strong_0.05/%s_mash_sex_xy_DEG_results.RData',ct))


