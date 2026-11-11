library(tidyverse)

################################################################################
# Results!  
################################################################################
df_nebula_filter = read.csv('/scratch/klwatki5/R_analysis/bbi/nebula/infants_removed/ir_add_nebula_size_filt.csv')
###############################################################################
rm(list = setdiff(ls(), "df_nebula_filter"))

# Define the cell type
cell_type = 'oligodendrocytes'

################################################################################
# SEX DEG RESLTS - AUTOSOMES
################################################################################

# Define the base directory
base_dir = '/scratch/klwatki5/R_analysis/bbi/mashr/results/1_raw_matrices/autosomes/'

#Load the new sex DEG results 
load(file.path(base_dir, paste0(cell_type,'_mash_sex_auto_DEG_results.RData')))

# Load OPC/olig results with different path; had to re-run later
# load('/scratch/klwatki5/R_analysis/bbi/mashr/results/1_raw_matrices/autosomes/oligodendrocytes_mash_sex_auto_DEG_results_no_OPC.RData')
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
var_info <- read.csv('/scratch/klwatki5/R_analysis/bbi/all_var_meta.txt.gz')[-1]

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
# Filtering 
################################################################################
# Add the column needed for combining Nebula and mashr data frames
df_raw <- df %>%  mutate(id = paste0(cell_type, "_louvain", louvain, "_", region), gene.id = paste0(ensembl_gene_id, "-", id)) 
df_sex_nebula_filt = left_join(df_raw, df_nebula_filter, by = c('gene.id', 'id')) %>% filter(!is.na(logFC_sexM))

# Clean columns
df_sex_nebula_filt <- df_sex_nebula_filt %>%
  dplyr::select(-gene.id, -gene, -external_gene_name.y, -chr.y, -gene.id.1, -cell_type.y, -louvain.y, -region.y) %>%
  rename("cell_type" = cell_type.x) %>%
  rename("external_gene_name" = external_gene_name.x) %>%
  rename("chr" = chr.x) %>% 
  rename("louvain" = louvain.x) %>% 
  rename("region" = region.x) 

# Check dimensions 
nrow(unique(df_sex_nebula_filt[, c("region", "cell_type","louvain")]))
length(unique(df_sex_nebula_filt$gene_id))
length(unique(df_sex_nebula_filt$ensembl_gene_id))
summary(df_sex_nebula_filt$sexM_lfsr)
summary(df_sex_nebula_filt$sexM_pm) 

# Save the combined autosome object to merge with XY dataframe
sex_auto_nebula_df = df_sex_nebula_filt

################################################################################
# SEX DEG RESLTS - XY CHROMOSOMES
################################################################################

# Define the base directory
base_dir = '/scratch/klwatki5/R_analysis/bbi/mashr/results/1_raw_matrices/sex_chr/'

#Load the new sex DEG results
load(file.path(base_dir, paste0(cell_type,'_mash_sex_xy_DEG_results.RData')))

#Special case for oligs and OPCs 
# load('/scratch/klwatki5/R_analysis/bbi/mashr/results/1_raw_matrices/sex_chr/oligodendrocytes_mash_sex_xy_DEG_results.RData')

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
var_info <- read.csv('/scratch/klwatki5/R_analysis/bbi/all_var_meta.txt.gz')[-1]

df <- df %>%
  left_join(var_info[, c("ensembl_gene_id", "external_gene_name","chr")],
            by = "ensembl_gene_id")

# Check dimensions 
nrow(unique(df[, c("region", "cell_type","louvain")]))
length(unique(df$gene_id)) # gene_id = gene + CTRC = unique tests in this dataset, should = nrow 
length(unique(df$ensembl_gene_id))
summary(df$sexM_lfsr)
summary(df$sexM_pm) 

# Bind with autosome raw df 
df_raw_merged = rbind(df, sex_auto_raw_df)
# Save the combined dataframe
write.csv(df_raw_merged, file=paste0('/scratch/klwatki5/R_analysis/bbi/mashr/results/2_raw_celltype_dataframes/ir_', cell_type, '_sex_DEG_mashr_df.csv'))

################################################################################
# Filtering 
################################################################################
# Add the column needed for combining Nebula and mashr data frames
df_raw <- df %>%  mutate(id = paste0(cell_type, "_louvain", louvain, "_", region), gene.id = paste0(ensembl_gene_id, "-", id)) 
df_sex_nebula_filt = left_join(df_raw, df_nebula_filter, by = c('gene.id', 'id')) %>% filter(!is.na(logFC_sexM))

# Clean columns
df_sex_nebula_filt <- df_sex_nebula_filt %>%
  dplyr::select(-gene.id, -gene, -external_gene_name.y, -chr.y, -gene.id.1, -cell_type.y, -louvain.y, -region.y) %>%
  rename("cell_type" = cell_type.x) %>%
  rename("external_gene_name" = external_gene_name.x) %>%
  rename("chr" = chr.x) %>% 
  rename("louvain" = louvain.x) %>% 
  rename("region" = region.x) 

# Check dimensions 
nrow(unique(df_sex_nebula_filt[, c("region", "cell_type","louvain")]))
length(unique(df_sex_nebula_filt$gene_id))
length(unique(df_sex_nebula_filt$ensembl_gene_id))
summary(df_sex_nebula_filt$sexM_lfsr)
summary(df_sex_nebula_filt$sexM_pm) 
length(unique(df_sex_nebula_filt$id))
df_filt_merged = rbind(df_sex_nebula_filt, sex_auto_nebula_df)
write_csv(df_filt_merged, file = paste0('/scratch/klwatki5/R_analysis/bbi/mashr/results/3_nebula_filtered_celltype_dataframes/',cell_type,'_sex_nebula_filt.csv'))
################################################################################
# MERGED DATAFRAMES 
################################################################################
# Read in all Nebula-filtered cell type files, combine into one sex result file.
################################################################################
astrocytes_sex = read_csv('/scratch/klwatki5/R_analysis/bbi/mashr/results/3_nebula_filtered_celltype_dataframes/astrocytes_sex_nebula_filt.csv')
astrocytes_sig = astrocytes_sex %>% filter(sexM_lfsr < 0.05)
length(unique(astrocytes_sex$ensembl_gene_id))
length(unique(astrocytes_sig$ensembl_gene_id))

basket_sex = read_csv('/scratch/klwatki5/R_analysis/bbi/mashr/results/3_nebula_filtered_celltype_dataframes/basket_sex_nebula_filt.csv')
basket_sig = basket_sex %>% filter(sexM_lfsr < 0.05)
length(unique(basket_sex$ensembl_gene_id))
length(unique(basket_sig$ensembl_gene_id))

cerebellar_sex = read_csv('/scratch/klwatki5/R_analysis/bbi/mashr/results/3_nebula_filtered_celltype_dataframes/cerebellar_sex_nebula_filt.csv')
cerebellar_sig = cerebellar_sex %>% filter(sexM_lfsr < 0.05)
length(unique(cerebellar_sex$ensembl_gene_id))
length(unique(cerebellar_sig$ensembl_gene_id))

ependymal_sex = read_csv('/scratch/klwatki5/R_analysis/bbi/mashr/results/3_nebula_filtered_celltype_dataframes/ependymal_sex_nebula_filt.csv')
ependymal_sig = ependymal_sex %>% filter(sexM_lfsr < 0.05)
length(unique(ependymal_sex$ensembl_gene_id))
length(unique(ependymal_sig$ensembl_gene_id))

gaba_sex = read_csv('/scratch/klwatki5/R_analysis/bbi/mashr/results/3_nebula_filtered_celltype_dataframes/gaba_sex_nebula_filt.csv')
gaba_sig = gaba_sex %>% filter(sexM_lfsr < 0.05)
length(unique(gaba_sex$ensembl_gene_id))
length(unique(gaba_sig$ensembl_gene_id))

glut_sex = read_csv('/scratch/klwatki5/R_analysis/bbi/mashr/results/3_nebula_filtered_celltype_dataframes/glut_sex_nebula_filt.csv')
glut_sig = glut_sex %>% filter(sexM_lfsr < 0.05)
length(unique(glut_sex$ensembl_gene_id))
length(unique(glut_sig$ensembl_gene_id))

msn_sex = read_csv('/scratch/klwatki5/R_analysis/bbi/mashr/results/3_nebula_filtered_celltype_dataframes/msn_sex_nebula_filt.csv')
msn_sig = msn_sex %>% filter(sexM_lfsr < 0.05)
length(unique(msn_sex$ensembl_gene_id))
length(unique(msn_sig$ensembl_gene_id))

microglia_sex = read_csv('/scratch/klwatki5/R_analysis/bbi/mashr/results/3_nebula_filtered_celltype_dataframes/microglia_sex_nebula_filt.csv')
microglia_sig = microglia_sex %>% filter(sexM_lfsr < 0.05)
length(unique(microglia_sex$ensembl_gene_id))
length(unique(microglia_sig$ensembl_gene_id))

midbrain_sex = read_csv('/scratch/klwatki5/R_analysis/bbi/mashr/results/3_nebula_filtered_celltype_dataframes/midbrain_sex_nebula_filt.csv')
midbrain_sig = midbrain_sex %>% filter(sexM_lfsr < 0.05)
length(unique(midbrain_sex$ensembl_gene_id))
length(unique(midbrain_sig$ensembl_gene_id))

oligodendrocytes_sex = read_csv('/scratch/klwatki5/R_analysis/bbi/mashr/results/3_nebula_filtered_celltype_dataframes/oligodendrocytes_sex_nebula_filt.csv')
oligodendrocytes_sig = oligodendrocytes_sex %>% filter(sexM_lfsr < 0.05)
length(unique(oligodendrocytes_sex$ensembl_gene_id))
length(unique(oligodendrocytes_sig$ensembl_gene_id))

OPC_sex = read_csv('/scratch/klwatki5/R_analysis/bbi/mashr/results/3_nebula_filtered_celltype_dataframes/OPCs_sex_nebula_filt.csv')
OPC_sig = OPC_sex %>% filter(sexM_lfsr < 0.05)
length(unique(OPC_sex$ensembl_gene_id))
length(unique(OPC_sig$ensembl_gene_id))

vascular_sex = read_csv('/scratch/klwatki5/R_analysis/bbi/mashr/results/3_nebula_filtered_celltype_dataframes/vascular_sex_nebula_filt.csv')
vascular_sig = vascular_sex %>% filter(sexM_lfsr < 0.05)
length(unique(vascular_sex$ensembl_gene_id))
length(unique(vascular_sig$ensembl_gene_id))

sex_merged = rbind(astrocytes_sex, basket_sex, cerebellar_sex, ependymal_sex, gaba_sex, glut_sex, microglia_sex, midbrain_sex, msn_sex, oligodendrocytes_sex, OPC_sex, vascular_sex)

write_csv(sex_merged, file= '/scratch/klwatki5/R_analysis/bbi/mashr/results/4_merged_dataframes/sex_merged.csv')

################################################################################
# Merge with age results to create master DEG results file 
################################################################################

rm(list = setdiff(ls(), "sex_merged"))

age_merged = read_csv("/scratch/klwatki5/R_analysis/bbi/mashr/results/4_merged_dataframes/age_merged.csv") # Identical to Figure_4 file 

# 1. Identify the unique columns from each dataset
age_specific_cols <- c("age_lfsr", "age_pm", "age_psd")
sex_specific_cols <- c("sexM_lfsr", "sexM_pm", "sexM_psd")

# 2. Get the common columns (removing the .x/.y suffix)
common_cols <- gsub("\\.x$", "", colnames(age_merged))
common_cols <- common_cols[!common_cols %in% c(age_specific_cols, sex_specific_cols)]

# 3. Create the clean merge
master_merged <- merge(
  # From age_merged, select only common cols and age-specific cols
  age_merged[, c(common_cols, age_specific_cols)],
  # From sex_merged, select only the key and sex-specific cols
  sex_merged[, c("gene_id", sex_specific_cols)],
  by = "gene_id",
  all = TRUE
)

master_merged = master_merged %>% mutate(
  quadrant = case_when(
    (age_pm > 0 & sexM_pm > 0) ~ 1,
    (age_pm < 0 & sexM_pm > 0) ~ 2,
    (age_pm < 0 & sexM_pm < 0) ~ 3,
    (age_pm > 0 & sexM_pm < 0) ~ 4)) # Assign quadrants of effects to all genes for a+sDEG analysis 

write_csv(master_merged, file='/scratch/klwatki5/R_analysis/bbi/mashr/results/4_merged_dataframes/master_merged.csv') # 5,146,072 rows 
write_csv(master_merged, file='/data/CEM/smacklab/data/bican/u01_results/mashr_additive_master_merged.csv')
