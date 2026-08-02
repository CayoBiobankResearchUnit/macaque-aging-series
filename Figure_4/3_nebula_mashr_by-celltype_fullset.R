#' @file 3_nebula_mashr_by-celltype_fullset.R
#' @title Empirical Bayes Correction of nebula output by cell class
#' @author Wei Yang
#'
#' @description
#' Applies Multivariate Adaptive Shrinkage (\code{mashr}) to perform Empirical Bayes 
#' correction on effect sizes (logFC) and standard errors obtained from \code{nebula}.
#'
#' @section Inputs:
#' \itemize{
#'   \item \code{<cell class>_mash_input.RData}: matrices of logFC and se prepared for mashR 
#'   input from 2_nebula_mashr_prepare_files.R
#' }
#'
#' @section Outputs:
#' \itemize{
#'   \item \code{<cell class>_mash_age_results.RData}: matrices of LFSR, pm and psd from mashR. 
#' }
#'
#' @usage
#' #SBATCH -N 1
#' #SBATCH -c 1  
#' #SBATCH -t 0-1:00:00   
#' #SBATCH --mem=50G
#' module load gsl-2.7.1-gcc-11.2.0
#' module load r-4.2.2-gcc-11.2.0
#' module load gcc-11.2.0-gcc-11.2.0
#' Rscript nebula_mashr_by-celltype_fullset.R $cell_type
NULL

#for cell type
args <- commandArgs(trailingOnly = TRUE)
ct = args[1]

library(tidyverse)
library(dplyr)
library(stringr)
library(ashr)
library(mashr)

set.seed(2025)


load(sprintf('/scratch/wyang114/outputs/nebula/by_region_NBGMM_mash_inputs/cpc_p05/%s_mash_input.RData',ct))
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

save(lfsr,pm,psd,file='/scratch/wyang114/outputs/nebula/mash_outputs/%s_mash_age_results.RData',ct)
