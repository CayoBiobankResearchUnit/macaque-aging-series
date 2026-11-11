library(tidyverse)
library(dplyr)
library(stringr)

df = read.csv('/scratch/klwatki5/R_analysis/bbi/nebula/infants_removed/ir_add_nebula_size_filt.csv') # 7,686,734
df = df %>% mutate(cell_type = case_when(
  (cell_type == "oligodendrocytes" & louvain %in% c(12,13)) ~ "OPC",
  TRUE ~ cell_type)) %>% 
  mutate(id = paste0(cell_type,"_louvain",louvain,"_", region)) %>%
  mutate(gene.id = paste0(gene,"_",id))

# rm(list = setdiff(ls(), "df"))

################################################################################
# Re-do size filtering with OPC and oligs separately 
################################################################################
####Use 30 donors and 300 cells as cutoff. (filter all basket cells in MB)
cell_count = read_csv("/scratch/klwatki5/R_analysis/bbi/nebula/infants_removed/cell_counts_updated.csv")
cell_count = cell_count %>% mutate(cell_type = case_when(
  (cell_type == "oligodendrocytes" & louvain %in% c(12,13)) ~ "OPC",
  TRUE ~ cell_type)) %>% 
  mutate(identity = paste0(cell_type,"_louvain",louvain,"_", region)) 

all_identity = cell_count$identity%>%unique()

e.filter.others = cell_count%>% filter(unique_animal_count<30 | cell_count < 300)%>%pull(identity)

e.filter.basket = cell_count%>%filter(cell_type=='basket cells')%>%
  filter(region=='MB'| unique_animal_count<30 | cell_count < 300)%>%pull(identity) 
e.filter = c(e.filter.basket,e.filter.others)

###############
#filtering
df_filter = df %>% 
  filter(!id %in% e.filter) %>%         # filter out region x subtype lower than threshold
  filter(!is.na(p_sexM))                # filter out genes with NA in p_age or p_sex
df_filter = df_filter %>% filter(!(cell_type == "medium spiny neurons" & region == "HIP")) # 5,146,072

unique(df_filter$cell_type) # 12
length(unique(df_filter$gene)) # 15117
length(unique(df_filter$id)) # 844

ct_list = df_filter$cell_type%>%unique

################################################################################
################################################################################
################################################################################
### UP TO THIS POINT, AGE AND SEX DEG ANALYSIS IS IDENTICAL. SEE "Figure_4" DIRECTORY FOR PREVIOUS STEPS (NEBULA)
################################################################################
################################################################################
################################################################################
# Read in gene metadata 
var = read.csv('/scratch/klwatki5/R_analysis/bbi/all_var_meta.txt.gz')[-1] 

# Identify genes on sex chromosomes for filtering 
sex_chr_genes = var %>% filter(chr %in% c('X','Y'))

# SEX
for(ct in ct_list){
  df_subset = df_filter%>%filter(cell_type==ct)
  coef_mat = df_subset%>%pivot_wider(id_cols = gene,names_from = id, values_from = logFC_sexM)%>%
    column_to_rownames('gene')
  #Set coefficient of missing value to 0
  coef_mat[is.na(coef_mat)]=0
  sum(is.na(coef_mat))
  
  se_mat = df_subset%>%pivot_wider(id_cols = gene,names_from = id, values_from = se_sexM)%>%
    column_to_rownames('gene')
  #Set SE of missing value to 1000
  se_mat[is.na(se_mat)]=1000
  sum(is.na(se_mat))
  
  print(ct)
  print(dim(coef_mat))
  print(dim(se_mat))

  # AUTOSOME FILES 
  auto_coef_mat = coef_mat[!rownames(coef_mat) %in% sex_chr_genes$ensembl_gene_id, ]
  auto_se_mat = se_mat[!rownames(se_mat) %in% sex_chr_genes$ensembl_gene_id, ]
  
  save(auto_coef_mat,auto_se_mat,file=(paste0('/scratch/klwatki5/R_analysis/bbi/mashr/mash_new/input/',ct,'_mash_sex_auto_DEG_input.Rdata')))  
  
  # XY CHROMOSOME FILES 
  xy_coef_mat = coef_mat[rownames(coef_mat) %in% sex_chr_genes$ensembl_gene_id, ]
  xy_se_mat = se_mat[rownames(se_mat) %in% sex_chr_genes$ensembl_gene_id, ]
  
  save(xy_coef_mat,xy_se_mat,file=(paste0('/scratch/klwatki5/R_analysis/bbi/mashr/mash_new/input/',ct,'_mash_sex_xy_DEG_input.RData')))

}
