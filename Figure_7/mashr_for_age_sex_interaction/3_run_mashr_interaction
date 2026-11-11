library(tidyverse)
library(dplyr)
library(stringr)
library(ashr)
library(mashr)

#for cell type
args <- commandArgs(trailingOnly = TRUE)
ct = args[1]

################################################################################
# AGE
################################################################################
load(paste0('/scratch/klwatki5/R_analysis/bbi/mashr/mash_new/interaction/input/',ct,'_mash_age_DEG_input.RData'))
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
strong = get_significant_results(m.1by1,0.1)

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
save(data,m.1by1,U.pca,U.ed,m.ed, file=sprintf('/scratch/klwatki5/R_analysis/bbi/mashr/mash_new/interaction/results/full_mashr_output/%s_mash_age_objects.RData',ct))

#get posterior summary
lfsr = get_lfsr(m.ed)
pm = get_pm(m.ed)
psd = get_psd(m.ed)

print(paste0('strong set size:',length(strong)))
print(paste0('log-likelihood',get_loglik(m.ed),digits = 10))

save(lfsr,pm,psd,file=sprintf('/scratch/klwatki5/R_analysis/bbi/mashr/mash_new/interaction/results/%s_mash_age_DEG_results.RData',ct))

################################################################################
# SEX - AUTOSOMES
################################################################################
load(paste0('/scratch/klwatki5/R_analysis/bbi/mashr/mash_new/interaction/input/',ct,'_mash_sex_auto_DEG_input.RData'))

# Read in data
Bhat = auto_coef_mat%>%as.matrix()
Shat = auto_se_mat%>%as.matrix()
data   = mash_set_data(Bhat, Shat)

#Select strong signals
m.1by1 = mash_1by1(data)
strong = get_significant_results(m.1by1,0.1)
length(strong)

#Obtain initial data-driven covariance matrices
U.pca = cov_pca(data,5,subset=strong)
print(names(U.pca))

#Apply extreme deconvolution
U.ed = cov_ed(data, U.pca, subset=strong)

#Run mash-fit the model
m.ed = mash(data, U.ed)
save(data,m.1by1,U.pca,U.ed,m.ed, file=sprintf('/scratch/klwatki5/R_analysis/bbi/mashr/mash_new/interaction/results/full_mashr_output/%s_mash_sex_auto_objects.RData',ct))

#get posterior summary
lfsr = get_lfsr(m.ed)
pm = get_pm(m.ed)
psd = get_psd(m.ed)

print(paste0('strong set size:',length(strong)))
print(paste0('log-likelihood',get_loglik(m.ed),digits = 10))

save(lfsr,pm,psd,file=sprintf('/scratch/klwatki5/R_analysis/bbi/mashr/mash_new/interaction/results/%s_mash_sex_auto_DEG_results.RData',ct))

################################################################################
# SEX - XY
################################################################################
load(paste0('/scratch/klwatki5/R_analysis/bbi/mashr/mash_new/interaction/input/',ct,'_mash_sex_xy_DEG_input.RData'))

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
strong = get_significant_results(m.1by1,0.1)
length(strong)

#Obtain initial data-driven covariance matrices
U.pca = cov_pca(data,5,subset=strong)
print(names(U.pca))

#Apply extreme deconvolution
U.ed = cov_ed(data, U.pca, subset=strong)

#Run mash-fit the model
m.ed = mash(data, U.ed)
save(data,m.1by1,U.pca,U.ed,m.ed, file=sprintf('/scratch/klwatki5/R_analysis/bbi/mashr/mash_new/interaction/results/full_mashr_output/%s_mash_sex_xy_objects.RData',ct))

#get posterior summary
lfsr = get_lfsr(m.ed)
pm = get_pm(m.ed)
psd = get_psd(m.ed)

print(paste0('strong set size:',length(strong)))
print(paste0('log-likelihood',get_loglik(m.ed),digits = 10))

save(lfsr,pm,psd,file=sprintf('/scratch/klwatki5/R_analysis/bbi/mashr/mash_new/interaction/results/%s_mash_sex_xy_DEG_results.RData',ct))

################################################################################
# INTERACTION TERM
################################################################################
load(paste0('/scratch/klwatki5/R_analysis/bbi/mashr/mash_new/interaction/input/',ct,'_mash_interaction_DEG_input.RData'))
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
strong = get_significant_results(m.1by1,0.1)

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
save(data,m.1by1,U.pca,U.ed,m.ed, file=sprintf('/scratch/klwatki5/R_analysis/bbi/mashr/mash_new/interaction/results/full_mashr_output/%s_mash_interaction_objects.RData',ct)) 

#get posterior summary
lfsr = get_lfsr(m.ed)
pm = get_pm(m.ed) 
psd = get_psd(m.ed)

print(paste0('strong set size:',length(strong)))
print(paste0('log-likelihood',get_loglik(m.ed),digits = 10))

save(lfsr,pm,psd,file=sprintf('/scratch/klwatki5/R_analysis/bbi/mashr/mash_new/interaction/results/%s_mash_interaction_DEG_results.RData',ct))

################################################################################
# INTERCEPT
################################################################################
load(paste0('/scratch/klwatki5/R_analysis/bbi/mashr/mash_new/interaction/input/',ct,'_mash_intercept_DEG_input.RData'))
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
strong = get_significant_results(m.1by1,0.1)

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
save(data,m.1by1,U.pca,U.ed,m.ed, file=sprintf('/scratch/klwatki5/R_analysis/bbi/mashr/mash_new/interaction/results/full_mashr_output/%s_mash_intercept_objects.RData',ct))

#get posterior summary
lfsr = get_lfsr(m.ed)
pm = get_pm(m.ed)
psd = get_psd(m.ed)

print(paste0('strong set size:',length(strong)))
print(paste0('log-likelihood',get_loglik(m.ed),digits = 10))

save(lfsr,pm,psd,file=sprintf('/scratch/klwatki5/R_analysis/bbi/mashr/mash_new/interaction/results/%s_mash_intercept_DEG_results.RData',ct))

