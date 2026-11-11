################################################################################
# Raw Results!  
################################################################################
var_info <- read.csv('/scratch/klwatki5/R_analysis/bbi/all_var_meta.txt.gz')[-1]

# Define the cell type
cell_type = 'vascular cells'
################################################################################
# AGE DEG RESLTS
################################################################################
# Define the base directory
base_dir = '/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/1_raw_matrices'

#Load the new age DEG results
load(file.path(base_dir, paste0(cell_type, '_mash_age_DEG_results.RData')))

# Convert the age matrices to data frames 
age_lfsr_df <- as.data.frame(lfsr)
age_pm_df <- as.data.frame(pm)
age_psd_df <- as.data.frame(psd)

# Add the gene names as a column
age_lfsr_df$gene <- rownames(age_lfsr_df)
age_pm_df$gene <- rownames(age_pm_df) 
age_psd_df$gene <- rownames(age_psd_df)

# Gather the age data into a long format
age_df <- age_lfsr_df %>%
  dplyr::select(gene, everything()) %>%
  pivot_longer(
    cols = -gene,
    names_to = c("cell_type", "louvain", "region"), 
    names_pattern = "([^_]+)_([^_]+)_([^_]+)",
    values_to = "age_lfsr"
  ) %>%
  inner_join(age_pm_df %>% 
               dplyr::select(gene, everything()) %>%
               pivot_longer(
                 cols = -gene,
                 names_to = c("cell_type", "louvain", "region"),
                 names_pattern = "([^_]+)_([^_]+)_([^_]+)", 
                 values_to = "age_pm"
               ), 
             by = c("gene", "cell_type", "louvain", "region")) %>%
  inner_join(age_psd_df %>%
               dplyr::select(gene, everything()) %>%
               pivot_longer(
                 cols = -gene,
                 names_to = c("cell_type", "louvain", "region"),
                 names_pattern = "([^_]+)_([^_]+)_([^_]+)",
                 values_to = "age_psd"
               ),
             by = c("gene", "cell_type", "louvain", "region")) %>%
  mutate(louvain = str_extract(louvain, "\\d+"),
         gene_id = paste(gene, cell_type, louvain, region, sep = "_"))

df = age_df %>% rename("ensembl_gene_id" = gene)

# Read in gene metadata and add external gene names
var_info <- read.csv('/scratch/klwatki5/R_analysis/bbi/all_var_meta.txt.gz')[-1]

df <- df %>%
  left_join(var_info[, c("ensembl_gene_id", "external_gene_name","chr")],
            by = "ensembl_gene_id")

# Check dimensions 
nrow(unique(df[, c("region", "cell_type","louvain")]))
length(unique(df$gene_id)) # gene_id = gene + CTRC = unique tests in this dataset, should = nrow 
length(unique(df$ensembl_gene_id))
summary(df$age_lfsr)
summary(df$age_pm) 

# Save the combined dataframe
write.csv(df, file=paste0('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/2_raw_celltype_dataframes/', cell_type, '_age_DEG_mashr_df.csv'))

################################################################################
# SEX DEG RESLTS - AUTOSOMES
################################################################################

# Define the base directory
base_dir = '/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/1_raw_matrices/sex_auto/'

#Load the new sex DEG results
#load(file.path(base_dir, paste0(cell_type,'_mash_sex_auto_raw_results.RData')))
load(file.path(base_dir, paste0(cell_type,'_mash_sex_auto_DEG_results.RData')))
# Convert the sex matrices to data frames
sexM_lfsr_df <- as.data.frame(lfsr)
sexM_pm_df <- as.data.frame(pm)
sexM_psd_df <- as.data.frame(psd)

# Add the gene names as a column
sexM_lfsr_df$gene <- rownames(sexM_lfsr_df)
sexM_pm_df$gene <- rownames(sexM_pm_df)
sexM_psd_df$gene <- rownames(sexM_psd_df)

# Gather the sex data into a long format
sex_df <- sexM_lfsr_df %>%
  dplyr::select(gene, everything()) %>%
  pivot_longer(
    cols = -gene,
    names_to = c("cell_type", "louvain", "region"),
    names_pattern = "([^_]+)_([^_]+)_([^_]+)",
    values_to = "sexM_lfsr"
  ) %>%
  inner_join(sexM_pm_df %>%
               dplyr::select(gene, everything()) %>%
               pivot_longer(
                 cols = -gene,
                 names_to = c("cell_type", "louvain", "region"),
                 names_pattern = "([^_]+)_([^_]+)_([^_]+)",
                 values_to = "sexM_pm"
               ),
             by = c("gene", "cell_type", "louvain", "region")) %>%
  inner_join(sexM_psd_df %>%
               dplyr::select(gene, everything()) %>%
               pivot_longer(
                 cols = -gene,
                 names_to = c("cell_type", "louvain", "region"),
                 names_pattern = "([^_]+)_([^_]+)_([^_]+)",
                 values_to = "sexM_psd"
               ),
             by = c("gene", "cell_type", "louvain", "region")) %>%
  mutate(louvain = str_extract(louvain, "\\d+"),
         gene_id = paste(gene, cell_type, louvain, region, sep = "_"))

df = sex_df %>% rename("ensembl_gene_id" = gene)

# Read in gene metadata and add external gene names

df <- df %>%
  left_join(var_info[, c("ensembl_gene_id", "external_gene_name","chr")],
            by = "ensembl_gene_id")
# Check dimensions
nrow(unique(df[, c("region", "cell_type","louvain")]))
length(unique(df$gene_id)) # gene_id = gene + CTRC = unique tests in this dataset, should = nrow
length(unique(df$ensembl_gene_id))
summary(df$sexM_lfsr)
summary(df$sexM_pm)

# Save the combined autosome object to merge with XY dataframe
sex_auto_raw_df  = df

################################################################################
# SEX DEG RESLTS - XY CHROMOSOMES
################################################################################

# Define the base directory
base_dir = '/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/1_raw_matrices/sex_xy/'

#Load the new sex DEG results
load(file.path(base_dir, paste0(cell_type,'_mash_sex_xy_DEG_results.RData')))

# Convert the sex matrices to data frames
sexM_lfsr_df <- as.data.frame(lfsr)
sexM_pm_df <- as.data.frame(pm)
sexM_psd_df <- as.data.frame(psd)

# Add the gene names as a column
sexM_lfsr_df$gene <- rownames(sexM_lfsr_df)
sexM_pm_df$gene <- rownames(sexM_pm_df)
sexM_psd_df$gene <- rownames(sexM_psd_df)

# Gather the sex data into a long format
sex_df <- sexM_lfsr_df %>%
  dplyr::select(gene, everything()) %>%
  pivot_longer(
    cols = -gene,
    names_to = c("cell_type", "louvain", "region"),
    names_pattern = "([^_]+)_([^_]+)_([^_]+)",
    values_to = "sexM_lfsr"
  ) %>%
  inner_join(sexM_pm_df %>%
               dplyr::select(gene, everything()) %>%
               pivot_longer(
                 cols = -gene,
                 names_to = c("cell_type", "louvain", "region"),
                 names_pattern = "([^_]+)_([^_]+)_([^_]+)",
                 values_to = "sexM_pm"
               ),
             by = c("gene", "cell_type", "louvain", "region")) %>%
  inner_join(sexM_psd_df %>%
               dplyr::select(gene, everything()) %>%
               pivot_longer(
                 cols = -gene,
                 names_to = c("cell_type", "louvain", "region"),
                 names_pattern = "([^_]+)_([^_]+)_([^_]+)",
                 values_to = "sexM_psd"
               ),
             by = c("gene", "cell_type", "louvain", "region")) %>%
  mutate(louvain = str_extract(louvain, "\\d+"),
         gene_id = paste(gene, cell_type, louvain, region, sep = "_"))

df = sex_df %>% rename("ensembl_gene_id" = gene)

df <- df %>%
  left_join(var_info[, c("ensembl_gene_id", "external_gene_name","chr")],
            by = "ensembl_gene_id")

# Check dimensions
nrow(unique(df[, c("region", "cell_type","louvain")]))
length(unique(df$gene_id)) # gene_id = gene + CTRC = unique tests in this dataset, should = nrow
length(unique(df$ensembl_gene_id))
summary(df$sexM_lfsr)
summary(df$sexM_pm)
length(unique(df$cell_type, df$louvain, df$region))

# Bind with autosome raw df
df_raw_merged = rbind(df, sex_auto_raw_df)
# Save the combined dataframe
write.csv(df_raw_merged, file=paste0('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/2_raw_celltype_dataframes/', cell_type, '_sex_DEG_mashr_df.csv'))

################################################################################
# INTERACTION TERM DEGS
################################################################################
# Define the base directory
base_dir = '/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/1_raw_matrices'

#Load the new age DEG results
load(file.path(base_dir, paste0(cell_type, '_mash_interaction_DEG_results.RData')))
load('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/1_raw_matrices/vascular cells_mash_interaction_DEG_results.RData')

# Convert the age matrices to data frames 
interaction_lfsr_df <- as.data.frame(lfsr)
interaction_pm_df <- as.data.frame(pm)
interaction_psd_df <- as.data.frame(psd)

# Add the gene names as a column
interaction_lfsr_df$gene <- rownames(interaction_lfsr_df)
interaction_pm_df$gene <- rownames(interaction_pm_df) 
interaction_psd_df$gene <- rownames(interaction_psd_df)

# Gather the age data into a long format
interaction_df <- interaction_lfsr_df %>%
  dplyr::select(gene, everything()) %>%
  pivot_longer(
    cols = -gene,
    names_to = c("cell_type", "louvain", "region"), 
    names_pattern = "([^_]+)_([^_]+)_([^_]+)",
    values_to = "interaction_lfsr"
  ) %>%
  inner_join(interaction_pm_df %>% 
               dplyr::select(gene, everything()) %>%
               pivot_longer(
                 cols = -gene,
                 names_to = c("cell_type", "louvain", "region"),
                 names_pattern = "([^_]+)_([^_]+)_([^_]+)", 
                 values_to = "interaction_pm"
               ), 
             by = c("gene", "cell_type", "louvain", "region")) %>%
  inner_join(interaction_psd_df %>%
               dplyr::select(gene, everything()) %>%
               pivot_longer(
                 cols = -gene,
                 names_to = c("cell_type", "louvain", "region"),
                 names_pattern = "([^_]+)_([^_]+)_([^_]+)",
                 values_to = "interaction_psd"
               ),
             by = c("gene", "cell_type", "louvain", "region")) %>%
  mutate(louvain = str_extract(louvain, "\\d+"),
         gene_id = paste(gene, cell_type, louvain, region, sep = "_"))

df = interaction_df %>% rename("ensembl_gene_id" = gene)

df <- df %>%
  left_join(var_info[, c("ensembl_gene_id", "external_gene_name","chr")],
            by = "ensembl_gene_id")

# Check dimensions 
nrow(unique(df[, c("region", "cell_type","louvain")]))
length(unique(df$gene_id)) # gene_id = gene + CTRC = unique tests in this dataset, should = nrow 
length(unique(df$ensembl_gene_id))
summary(df$interaction_lfsr)
summary(df$interaction_pm) 
any(is.na(df))
# Save the combined dataframe
write.csv(df, file=paste0('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/2_raw_celltype_dataframes/', cell_type, '_interaction_DEG_mashr_df.csv'))

################################################################################
# INTERCEPT TERM
################################################################################
# # Define the base directory
# base_dir = '/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/1_raw_matrices/intercept/'
# 
# #Load the new age DEG results
# load(file.path(base_dir, paste0(cell_type, '_mash_intercept_DEG_results.RData')))
# 
# 
# # Convert the age matrices to data frames 
# intercept_lfsr_df <- as.data.frame(lfsr)
# intercept_pm_df <- as.data.frame(pm)
# intercept_psd_df <- as.data.frame(psd)
# 
# # Add the gene names as a column
# intercept_lfsr_df$gene <- rownames(intercept_lfsr_df)
# intercept_pm_df$gene <- rownames(intercept_pm_df) 
# intercept_psd_df$gene <- rownames(intercept_psd_df)
# 
# # Gather the age data into a long format
# intercept_df <- intercept_lfsr_df %>%
#   dplyr::select(gene, everything()) %>%
#   pivot_longer(
#     cols = -gene,
#     names_to = c("cell_type", "louvain", "region"), 
#     names_pattern = "([^_]+)_([^_]+)_([^_]+)",
#     values_to = "intercept_lfsr"
#   ) %>%
#   inner_join(intercept_pm_df %>% 
#                dplyr::select(gene, everything()) %>%
#                pivot_longer(
#                  cols = -gene,
#                  names_to = c("cell_type", "louvain", "region"),
#                  names_pattern = "([^_]+)_([^_]+)_([^_]+)", 
#                  values_to = "intercept_pm"
#                ), 
#              by = c("gene", "cell_type", "louvain", "region")) %>%
#   inner_join(intercept_psd_df %>%
#                dplyr::select(gene, everything()) %>%
#                pivot_longer(
#                  cols = -gene,
#                  names_to = c("cell_type", "louvain", "region"),
#                  names_pattern = "([^_]+)_([^_]+)_([^_]+)",
#                  values_to = "intercept_psd"
#                ),
#              by = c("gene", "cell_type", "louvain", "region")) %>%
#   mutate(louvain = str_extract(louvain, "\\d+"),
#          gene_id = paste(gene, cell_type, louvain, region, sep = "_"))
# 
# df = intercept_df %>% rename("ensembl_gene_id" = gene)
# 
# df <- df %>%
#   left_join(var_info[, c("ensembl_gene_id", "external_gene_name","chr")],
#             by = "ensembl_gene_id")
# 
# # Check dimensions 
# nrow(unique(df[, c("region", "cell_type","louvain")]))
# length(unique(df$gene_id)) # gene_id = gene + CTRC = unique tests in this dataset, should = nrow 
# length(unique(df$ensembl_gene_id))
# summary(df$intercept_lfsr)
# summary(df$intercept_pm) 
# 
# # Save the combined dataframe
# write.csv(df, file=paste0('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/2_raw_celltype_dataframes/', cell_type, '_intercept_DEG_mashr_df.csv'))

###############################################################################
# NEXT STEP: FILTER MISSING AND CLEAN COLUMNS 
###############################################################################
df_nebula_filter = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/interaction_nebula_results_filtered_for_input.csv')
rm(list = setdiff(ls(), c("df_nebula_filter")))

cell_type = "vascular cells"
################################################################################
# AGE
################################################################################
# Add the column needed for combining Nebula and mashr data frames
df = read.csv(paste0('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/2_raw_celltype_dataframes/', cell_type, '_age_DEG_mashr_df.csv'))[-1]
df_raw <- df %>%  mutate(id = paste0(cell_type, "_louvain", louvain, "_", region), gene.id = paste0(ensembl_gene_id, "_", id)) 
df_age_nebula_filt = left_join(df_raw, df_nebula_filter, by = c('gene.id', 'id')) %>% filter(!is.na(logFC_age_centered))

# Clean columns
df_age_nebula_filt <- df_age_nebula_filt %>%
  dplyr::select( -gene.id, -gene, -external_gene_name.y, -chr.y, -cell_type.y, -region.y, -louvain.y) %>%
  rename("cell_type" = cell_type.x) %>%
  rename("external_gene_name" = external_gene_name.x) %>%
  rename("chr" = chr.x)  %>%
  rename("region" = region.x) %>% 
  rename("louvain" = louvain.x)

# Check dimensions 
nrow(df)
nrow(unique(df_age_nebula_filt[, c("region", "cell_type","louvain")]))
nrow(df_age_nebula_filt)
length(unique(df_age_nebula_filt$ensembl_gene_id))
summary(df_age_nebula_filt$age_lfsr)
summary(df_age_nebula_filt$age_pm) 

colnames(df_age_nebula_filt)

write_csv(df_age_nebula_filt, file = paste0('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/',cell_type,'_age_nebula_filt.csv'))

################################################################################
# SEX
################################################################################
# Add the column needed for combining Nebula and mashr data frames
df = read.csv(paste0('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/2_raw_celltype_dataframes/', cell_type, '_sex_DEG_mashr_df.csv'))[-1]
df_raw <- df %>%  mutate(id = paste0(cell_type, "_louvain", louvain, "_", region), gene.id = paste0(ensembl_gene_id, "_", id))
df_sex_nebula_filt = left_join(df_raw, df_nebula_filter, by = c('gene.id', 'id')) %>% filter(!is.na(logFC_sexM))

# Clean columns
df_sex_nebula_filt <- df_sex_nebula_filt %>%
  dplyr::select( -gene.id, -gene, -external_gene_name.y, -chr.y, -cell_type.y, -region.y, -louvain.y) %>%
  rename("cell_type" = cell_type.x) %>%
  rename("external_gene_name" = external_gene_name.x) %>%
  rename("chr" = chr.x)  %>%
  rename("region" = region.x) %>%
  rename("louvain" = louvain.x)

# Check dimensions
nrow(unique(df_sex_nebula_filt[, c("region", "cell_type","louvain")]))
length(unique(df_sex_nebula_filt$id))
length(unique(df_sex_nebula_filt$ensembl_gene_id))
summary(df_sex_nebula_filt$sexM_lfsr)
summary(df_sex_nebula_filt$sexM_pm)

write_csv(df_sex_nebula_filt, file = paste0('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/',cell_type,'_sex_nebula_filt.csv'))

################################################################################
# INTERACTION TERM
################################################################################
# Add the column needed for combining Nebula and mashr data frames
df = read.csv(paste0('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/2_raw_celltype_dataframes/', cell_type, '_interaction_DEG_mashr_df.csv'))[-1]
df_raw <- df %>%  mutate(id = paste0(cell_type, "_louvain", louvain, "_", region), gene.id = paste0(ensembl_gene_id, "_", id)) 
df_interaction_nebula_filt = left_join(df_raw, df_nebula_filter, by = c('gene.id', 'id')) %>% filter(!is.na(logFC_interaction))

# Clean columns
df_interaction_nebula_filt <- df_interaction_nebula_filt %>%
  dplyr::select( -gene.id, -gene, -external_gene_name.y, -chr.y, -cell_type.y, -region.y, -louvain.y) %>%
  rename("cell_type" = cell_type.x) %>%
  rename("external_gene_name" = external_gene_name.x) %>%
  rename("chr" = chr.x)  %>%
  rename("region" = region.x) %>% 
  rename("louvain" = louvain.x)

# Check dimensions 
nrow(unique(df_interaction_nebula_filt[, c("region", "cell_type","louvain")]))
length(unique(df_interaction_nebula_filt$id))
length(unique(df_interaction_nebula_filt$ensembl_gene_id))
summary(df_interaction_nebula_filt$interaction_lfsr)
summary(df_interaction_nebula_filt$interaction_pm) 
any(is.na(df_interaction_nebula_filt))

write_csv(df_interaction_nebula_filt, file = paste0('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/',cell_type,'_interaction_nebula_filt.csv'))

################################################################################
# INTERCEPT TERM
################################################################################
# Add the column needed for combining Nebula and mashr data frames
df = read.csv(paste0('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/2_raw_celltype_dataframes/', cell_type, '_intercept_DEG_mashr_df.csv'))[-1]
df_raw <- df %>%  mutate(id = paste0(cell_type, "_louvain", louvain, "_", region), gene.id = paste0(ensembl_gene_id, "_", id))
df_intercept_nebula_filt = left_join(df_raw, df_nebula_filter, by = c('gene.id', 'id')) %>% filter(!is.na(logFC_intercept))

# Clean columns
df_intercept_nebula_filt <- df_intercept_nebula_filt %>%
  dplyr::select( -gene.id, -gene, -external_gene_name.y, -chr.y, -cell_type.y, -region.y, -louvain.y) %>%
  rename("cell_type" = cell_type.x) %>%
  rename("external_gene_name" = external_gene_name.x) %>%
  rename("chr" = chr.x)  %>%
  rename("region" = region.x) %>%
  rename("louvain" = louvain.x)

# Check dimensions
nrow(unique(df_intercept_nebula_filt[, c("region", "cell_type","louvain")]))
length(unique(df_intercept_nebula_filt$gene_id))
length(unique(df_intercept_nebula_filt$ensembl_gene_id))
summary(df_intercept_nebula_filt$intercept_lfsr)
summary(df_intercept_nebula_filt$intercept_pm)
any(is.na(df_intercept_nebula_filt))

write_csv(df_intercept_nebula_filt, file = paste0('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/',cell_type,'_intercept_nebula_filt.csv'))

################################################################################
# Read in all Nebula-filtered cell type files, combine into one age result file. 
################################################################################
# AGE
astrocytes_age = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/astrocytes_age_nebula_filt.csv')
basket_age = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/basket cells_age_nebula_filt.csv')
cerebellar_age = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/cerebellar neurons_age_nebula_filt.csv')
ependymal_age = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/ependymal cells_age_nebula_filt.csv')
gaba_age = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/GABAergic neurons_age_nebula_filt.csv')
glut_age = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/glutamatergic neurons_age_nebula_filt.csv')
msn_age = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/medium spiny neurons_age_nebula_filt.csv')
microglia_age = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/microglia_age_nebula_filt.csv')
midbrain_age = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/midbrain neurons_age_nebula_filt.csv')
oligodendrocytes_age = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/oligodendrocytes_age_nebula_filt.csv')
OPC_age = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/OPC_age_nebula_filt.csv')
vascular_age=read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/vascular cells_age_nebula_filt.csv')


age_merged = rbind(astrocytes_age, basket_age, cerebellar_age, ependymal_age, gaba_age, glut_age, microglia_age, midbrain_age, msn_age, oligodendrocytes_age, OPC_age, vascular_age)
any(is.na(age_merged))
write_csv(age_merged, file= '/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/4_merged_dataframes/age_merged.csv')

age_sig = age_merged %>% filter(age_lfsr<0.05)
length(unique(age_sig$gene_id))
length(unique(age_sig$ensembl_gene_id))

################################################################################
# SEX
################################################################################
astrocytes_sex = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/astrocytes_sex_nebula_filt.csv')
basket_sex = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/basket cells_sex_nebula_filt.csv')
cerebellar_sex = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/cerebellar neurons_sex_nebula_filt.csv')
ependymal_sex = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/ependymal cells_sex_nebula_filt.csv')
gaba_sex = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/GABAergic neurons_sex_nebula_filt.csv')
glut_sex = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/glutamatergic neurons_sex_nebula_filt.csv')
msn_sex = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/medium spiny neurons_sex_nebula_filt.csv')
microglia_sex = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/microglia_sex_nebula_filt.csv')
midbrain_sex = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/midbrain neurons_sex_nebula_filt.csv')
oligodendrocytes_sex = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/oligodendrocytes_sex_nebula_filt.csv')
OPC_sex = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/OPC_sex_nebula_filt.csv')
vascular_sex=read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/vascular cells_sex_nebula_filt.csv')

sex_merged = rbind(astrocytes_sex, basket_sex, cerebellar_sex, ependymal_sex, gaba_sex, glut_sex, microglia_sex, midbrain_sex, msn_sex, oligodendrocytes_sex, OPC_sex, vascular_sex)
any(is.na(sex_merged))
write_csv(sex_merged, file= '/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/4_merged_dataframes/sex_merged.csv')

sex_sig= sex_merged %>% filter(sexM_lfsr<0.05)
length(unique(sex_sig$gene_id))
length(unique(sex_sig$ensembl_gene_id))

################################################################################
# INTERACTION
################################################################################
astrocytes_interaction = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/astrocytes_interaction_nebula_filt.csv')
basket_interaction = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/basket cells_interaction_nebula_filt.csv')
cerebellar_interaction = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/cerebellar neurons_interaction_nebula_filt.csv')
ependymal_interaction = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/ependymal cells_interaction_nebula_filt.csv')
gaba_interaction = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/GABAergic neurons_interaction_nebula_filt.csv')
glut_interaction = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/glutamatergic neurons_interaction_nebula_filt.csv')
msn_interaction = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/medium spiny neurons_interaction_nebula_filt.csv')
microglia_interaction = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/microglia_interaction_nebula_filt.csv')
midbrain_interaction = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/midbrain neurons_interaction_nebula_filt.csv')
oligodendrocytes_interaction = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/oligodendrocytes_interaction_nebula_filt.csv')
OPC_interaction = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/OPC_interaction_nebula_filt.csv')
vascular_interaction=read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/vascular cells_interaction_nebula_filt.csv')


interaction_merged = rbind(astrocytes_interaction, basket_interaction, cerebellar_interaction, ependymal_interaction, gaba_interaction, glut_interaction, microglia_interaction, midbrain_interaction, msn_interaction, oligodendrocytes_interaction, OPC_interaction, vascular_interaction)
any(is.na(interaction_merged))
write_csv(interaction_merged, file= '/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/4_merged_dataframes/interaction_merged.csv')

interaction_sig = interaction_merged %>% filter(interaction_lfsr<0.05)
length(unique(interaction_sig$gene_id))
length(unique(interaction_sig$ensembl_gene_id))


################################################################################
# Intercept
################################################################################
astrocytes_intercept = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/astrocytes_intercept_nebula_filt.csv')
basket_intercept = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/basket cells_intercept_nebula_filt.csv')
cerebellar_intercept = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/cerebellar neurons_intercept_nebula_filt.csv')
ependymal_intercept = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/ependymal cells_intercept_nebula_filt.csv')
gaba_intercept = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/GABAergic neurons_intercept_nebula_filt.csv')
glut_intercept = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/glutamatergic neurons_intercept_nebula_filt.csv')
msn_intercept = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/medium spiny neurons_intercept_nebula_filt.csv')
microglia_intercept = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/microglia_intercept_nebula_filt.csv')
midbrain_intercept = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/midbrain neurons_intercept_nebula_filt.csv')
oligodendrocytes_intercept = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/oligodendrocytes_intercept_nebula_filt.csv')
OPC_intercept = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/OPC_intercept_nebula_filt.csv')
vascular_intercept=read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/3_nebula_filtered_celltype_dataframes/vascular cells_intercept_nebula_filt.csv')


intercept_merged = rbind(astrocytes_intercept, basket_intercept, cerebellar_intercept, ependymal_intercept, gaba_intercept, glut_intercept, microglia_intercept, midbrain_intercept, msn_intercept, oligodendrocytes_intercept, OPC_intercept, vascular_intercept)
any(is.na(intercept_merged))
write_csv(intercept_merged, file= '/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/4_merged_dataframes/intercept_merged.csv')

################################################################################
# Combine all results 
################################################################################
age_merged = read_csv("/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/4_merged_dataframes/age_merged.csv")
sex_merged = read_csv("/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/4_merged_dataframes/sex_merged.csv")
int_merged = read_csv("/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/4_merged_dataframes/interaction_merged.csv")
icp_merged = read_csv("/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/4_merged_dataframes/intercept_merged.csv")

# 1. First identify the unique columns from each dataset
age_specific_cols <- c("age_lfsr", "age_pm", "age_psd")
sex_specific_cols <- c("sexM_lfsr", "sexM_pm", "sexM_psd")
int_specific_cols = c("interaction_lfsr", "interaction_pm", "interaction_psd")
icp_specific_cols = c("intercept_lfsr", "intercept_pm", "intercept_psd")


# 2. Get the common columns (removing the .x/.y suffix)
common_cols <- gsub("\\.x$", "", colnames(age_merged))
common_cols <- common_cols[!common_cols %in% c(age_specific_cols, sex_specific_cols, int_specific_cols)] #, icp_specific_cols)]

# 3. Create the clean merge
temp_merged <- merge(
  age_merged[, c(common_cols, age_specific_cols)],
  sex_merged[, c("gene_id", sex_specific_cols)],
  by = "gene_id",
  all = TRUE
)

# Then merge the result with int_merged
master_merged <- merge(
  temp_merged,
  int_merged[, c("gene_id", int_specific_cols)],
  by = "gene_id",
  all = TRUE
)

id_counts = master_merged %>% distinct(cell_type, id )
id_counts %>% 
  group_by(cell_type) %>% 
  summarise(unique_ids = n_distinct(id))

master_merged_2 <- merge(
  master_merged,
  icp_merged[, c("gene_id", icp_specific_cols)],
  by = "gene_id",
  all = TRUE
)

any(is.na(master_merged_2))

write_csv(master_merged_2, file='/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/4_merged_dataframes/master_merged_interaction.csv')

################################################################################
# Define effect direction quadrants 
################################################################################
df = read_csv('/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/4_merged_dataframes/master_merged_interaction.csv') # 597,359

df = df %>%
mutate(
  age_effect_female = age_pm,
  age_effect_male = age_pm + interaction_pm,
  
  # Which sex has the stronger (larger magnitude) age effect?
  stronger_in = case_when(
    abs(age_effect_male) > abs(age_effect_female) ~ "male",
    abs(age_effect_male) < abs(age_effect_female) ~ "female",
    TRUE ~ "equal"
  ),
  
  # Direction of the stronger effect
  stronger_direction = case_when(
    stronger_in == "male"   ~ sign(age_effect_male),
    stronger_in == "female" ~ sign(age_effect_female),
    TRUE ~ 0  # equal or zero slope
  )
)

df = df %>% mutate(quadrant = case_when(
  (stronger_in == "male" & age_effect_male > 0) ~ 1, 
  (stronger_in == "male" & age_effect_male < 0) ~ 2, 
  (stronger_in == "female" & age_effect_female < 0) ~ 3,
  (stronger_in == "female" & age_effect_female > 0) ~ 4,
))
df$quadrant = as.factor(df$quadrant)

interaction_summary = df %>%
  group_by(id) %>%  
  summarize(region = dplyr::first(region),           # Keep one value of region
            louvain = dplyr::first(louvain),        # Keep one value of louvain
            cell_type = dplyr::first(cell_type),    # Keep one value of cell_type
            total_genes = length(ensembl_gene_id),
            
            age_sig = sum(age_lfsr < 0.05),
            age_up = sum(age_lfsr < 0.05 & age_pm > 0),    # Only count if significant AND positive
            age_down = sum(age_lfsr < 0.05 & age_pm < 0),
            nebula_age_sig = sum(p_adj_age_centered < 0.05),
            
            sex_sig = sum(sexM_lfsr < 0.05),
            sex_M = sum(sexM_lfsr < 0.05 & sexM_pm > 0),    
            sex_F = sum(sexM_lfsr < 0.05 & sexM_pm < 0),
            sex_auto = sum(sexM_lfsr < 0.05 & !(chr %in% c("X", "Y"))),
            sex_X = sum(sexM_lfsr < 0.05 & chr == "X"),
            sex_Y = sum(sexM_lfsr < 0.05 & chr == "Y"),
            nebula_sex_sig = sum(p_adj_sexM < 0.05),
            
            interaction_sig = sum(interaction_lfsr < 0.05),
            interaction_q1 = sum(interaction_lfsr < 0.05 & quadrant == 1),
            interaction_q2 = sum(interaction_lfsr < 0.05 & quadrant == 2),
            interaction_q3 = sum(interaction_lfsr < 0.05 & quadrant == 3),
            interaction_q4 = sum(interaction_lfsr < 0.05 & quadrant == 4),
            nebula_interaction_sig = sum(p_adj_interaction < 0.05)
  ) %>%
  mutate(prop_aDEGs = age_sig/total_genes) %>% 
  mutate(prop_sDEGs = sex_sig/total_genes) %>%
  mutate(prop_interaction = interaction_sig/total_genes) %>%
  mutate(prop_q1 = interaction_q1/interaction_sig) %>% 
  mutate(prop_q2 = interaction_q2/interaction_sig) %>% 
  mutate(prop_q3 = interaction_q3/interaction_sig) %>% 
  mutate(prop_q4 = interaction_q4/interaction_sig) %>%
  mutate(prop_neb_int = nebula_interaction_sig/total_genes)

interaction_summary[is.na(interaction_summary)] = 0

write_csv(interaction_summary, '/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/4_merged_dataframes/interaction_summary_new_quadrants.csv')

write_csv(df, '/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/4_merged_dataframes/master_merged_new_quadrants.csv')
