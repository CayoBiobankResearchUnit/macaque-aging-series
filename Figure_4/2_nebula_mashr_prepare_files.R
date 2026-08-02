#' @file 2_nebula_mashr_prepare_files.R
#' @title Prepare matrices for Empirical Bayes Correction 
#' @author Wei Yang
#'
#' @description
#' Applies Multivariate Adaptive Shrinkage (\code{mashr}) to perform Empirical Bayes 
#' correction on effect sizes (logFC) and standard errors obtained from \code{nebula}.
#'
#' @section Inputs:
#' \itemize{
#'   \item \code{<cell class>_<subcluster>.txt}: data frame of nebula output with logFC, se and p-value, 
#'   adjusted p-value after FDR correction for each subcluster from 1_nebula_by_subtype-region.R.
#' }
#'
#' @section Outputs:
#' \itemize{
#'   \item \code{<cell class>_mash_input.RData}: matrices of logFC and se prepared for mashR input.
#' }
#'
#' @usage
#' Rscript 2_nebula_mashr_prepare_files.R
NULL

library(tidyverse)
library(dplyr)
library(stringr)

####organize nebula outputs from all subclusters into a dataframe
data_dir = '/scratch/wyang114/outputs/nebula/by_region_NBLMM/'
file_list <- list.files(path = data_dir, pattern = "\\.txt$", recursive = TRUE, full.names = TRUE)
df_list <- vector("list", length(file_list))
for (i in seq_along(file_list)) {
  df_list[[i]] <- read.table(file_list[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
}
combined_df <- do.call(rbind, df_list)
write.csv(combined_df, file = file.path(data_dir, "by_region_all_NBGMM_cpc.csv"), row.names = FALSE)


####Use 30 donors and 300 cells as cutoff. (filter all basket cells in MB)
df = read.csv(file.path(data_dir,'by_region_all_NBGMM_cpc.csv'))
cell_count = read.csv(file.path(data_dir,'combined_cell_type_region_counts.csv'))
cell_count$identity = paste0(cell_count$cell_type,'_louvain',cell_count$louvain,'_',cell_count$region)
all_identity = cell_count$identity%>%unique()
e.filter.others = cell_count%>%filter(cell_type!='basket')%>%
  filter(unique_animal_count<30 | cell_count < 300)%>%pull(identity) #removing 572 region x subtypes
e.filter.basket = cell_count%>%filter(cell_type=='basket')%>%
  filter(region=='MB'| unique_animal_count<30 | cell_count < 300)%>%pull(identity) #removing 13 region x subtypes
e.filter = c(e.filter.basket,e.filter.others)

###############split into celltype-based datasets
#filtering
df_filter = df%>%filter(cpc>=0.05)%>%  #filter out genes with lower cpc
  filter(!id %in% e.filter)%>%         # filter out region x subtype lower than threshold
  filter(!is.na(p_age))                # filter out genes with NA in p_age or p_sex


ct_list = df_filter$cell_type%>%unique
for(ct in ct_list){
  df_subset = df_filter%>%filter(cell_type==ct)
  print(min(df_subset$cpc))
  #Use logFC_age or logFC_sex to make coefficient matrix
  coef_mat = df_subset%>%pivot_wider(id_cols = gene,names_from = id, values_from = logFC_age)%>%
    column_to_rownames('gene')
  #Set coefficient of missing value to 0
  coef_mat[is.na(coef_mat)]=0
  sum(is.na(coef_mat))
  
  #Use se_age or se_sex to make coefficient matrix
  se_mat = df_subset%>%pivot_wider(id_cols = gene,names_from = id, values_from = se_age)%>%
    column_to_rownames('gene')
  #Set SE of missing value to 1000
  se_mat[is.na(se_mat)]=1000
  sum(is.na(se_mat))
  
  print(ct)
  print(dim(coef_mat))
  print(dim(se_mat))
  save(coef_mat,se_mat,file=sprintf('/scratch/wyang114/outputs/nebula/by_region_NBGMM_mash_inputs/cpc_p05/%s_mash_input.RData',ct))
}