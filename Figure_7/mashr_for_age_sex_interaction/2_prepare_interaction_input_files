library(tidyverse)
library(dplyr)
library(stringr)

df = read_csv("/scratch/klwatki5/R_analysis/bbi/mashr/mash_new/interaction/interaction_nebula_results_filtered_for_input.csv") # 597359 rows after filtering genes for main effect sig
var = read.csv('/scratch/klwatki5/R_analysis/bbi/all_var_meta.txt.gz')[-1] 
sex_chr_genes = var %>% filter(chr %in% c('X','Y'))

df_filter = df
rm(df)
ct_list = df_filter$cell_type%>%unique
# SEX
for(ct in ct_list){
  df_subset = df_filter%>%filter(cell_type==ct)
  #print(min(df_subset$cpc))
  #Use logFC_age or logFC_age to make coefficient matrix
  coef_mat = df_subset%>%pivot_wider(id_cols = gene,names_from = id, values_from = logFC_sexM)%>%
    column_to_rownames('gene')
  #Set coefficient of missing value to 0
  coef_mat[is.na(coef_mat)]=0
  sum(is.na(coef_mat))
  
  #Use se_age or se_age to make coefficient matrix
  se_mat = df_subset%>%pivot_wider(id_cols = gene,names_from = id, values_from = se_sexM)%>%
    column_to_rownames('gene')
  #Set SE of missing value to 1000
  se_mat[is.na(se_mat)]=1000
  sum(is.na(se_mat))
  
  print(ct)
  print(dim(coef_mat))
  print(dim(se_mat))
  auto_coef_mat = coef_mat[!rownames(coef_mat) %in% sex_chr_genes$ensembl_gene_id, ]
  auto_se_mat = se_mat[!rownames(se_mat) %in% sex_chr_genes$ensembl_gene_id, ]
  
  save(auto_coef_mat,auto_se_mat,file=(paste0('/scratch/klwatki5/R_analysis/bbi/mashr/mash_new/interaction/input/',ct,'_mash_sex_auto_DEG_input.RData')))
  
  # Make sex chr-only input 
  xy_coef_mat = coef_mat[rownames(coef_mat) %in% sex_chr_genes$ensembl_gene_id, ]
  xy_se_mat = se_mat[rownames(se_mat) %in% sex_chr_genes$ensembl_gene_id, ]
  
  save(xy_coef_mat,xy_se_mat,file=(paste0('/scratch/klwatki5/R_analysis/bbi/mashr/mash_new/interaction/input/',ct,'_mash_sex_xy_DEG_input.RData')))
  
}

# AGE
for(ct in ct_list){
  df_subset = df_filter%>%filter(cell_type==ct)
  print(min(df_subset$cpc))
  #Use logFC_age or logFC_age to make coefficient matrix
  coef_mat = df_subset%>%pivot_wider(id_cols = gene,names_from = id, values_from = logFC_age_centered)%>%
    column_to_rownames('gene')
  #Set coefficient of missing value to 0
  coef_mat[is.na(coef_mat)]=0
  sum(is.na(coef_mat))
  
  #Use se_age or se_age to make coefficient matrix
  se_mat = df_subset%>%pivot_wider(id_cols = gene,names_from = id, values_from = se_age_centered)%>%
    column_to_rownames('gene')
  #Set SE of missing value to 1000
  se_mat[is.na(se_mat)]=1000
  sum(is.na(se_mat))
  
  print(ct)
  print(dim(coef_mat))
  print(dim(se_mat))
  save(coef_mat,se_mat,file=sprintf('/scratch/klwatki5/R_analysis/bbi/mashr/mash_new/interaction/input/%s_mash_age_DEG_input.RData',ct))
}
################################################################################
# Interaction
################################################################################
# INTERACTION
for(ct in ct_list){
  df_subset = df_filter%>%filter(cell_type==ct)
  print(min(df_subset$cpc))
  #Use logFC_age or logFC_age to make coefficient matrix
  coef_mat = df_subset%>%pivot_wider(id_cols = gene,names_from = id, values_from = logFC_interaction)%>%
    column_to_rownames('gene')
  #Set coefficient of missing value to 0
  coef_mat[is.na(coef_mat)]=0
  sum(is.na(coef_mat))
  
  #Use standard error to make SE matrix
  se_mat = df_subset%>%pivot_wider(id_cols = gene,names_from = id, values_from = se_interaction)%>%
    column_to_rownames('gene')
  #Set SE of missing value to 1000
  se_mat[is.na(se_mat)]=1000
  sum(is.na(se_mat))
  
  print(ct)
  print(dim(coef_mat))
  print(dim(se_mat))
  
  save(coef_mat,se_mat,file=(paste0('/scratch/klwatki5/R_analysis/bbi/mashr/mash_new/interaction/input/',ct,'_mash_interaction_DEG_input.RData')))
}

################################################################################
# Intercept
################################################################################
# INTERCEPT
for(ct in ct_list){
  df_subset = df_filter%>%filter(cell_type==ct)
  print(min(df_subset$cpc))
  #Use logFC_age or logFC_age to make coefficient matrix
  coef_mat = df_subset%>%pivot_wider(id_cols = gene,names_from = id, values_from = logFC_intercept)%>%
    column_to_rownames('gene')
  #Set coefficient of missing value to 0
  coef_mat[is.na(coef_mat)]=0
  sum(is.na(coef_mat))
  
  #Use standard error to make SE matrix
  se_mat = df_subset%>%pivot_wider(id_cols = gene,names_from = id, values_from = se_intercept)%>%
    column_to_rownames('gene')
  #Set SE of missing value to 1000
  se_mat[is.na(se_mat)]=1000
  sum(is.na(se_mat))
  
  print(ct)
  print(dim(coef_mat))
  print(dim(se_mat))
  
  save(coef_mat,se_mat,file=(paste0('/scratch/klwatki5/R_analysis/bbi/mashr/mash_new/interaction/input/',ct,'_mash_intercept_DEG_input.RData')))
}
