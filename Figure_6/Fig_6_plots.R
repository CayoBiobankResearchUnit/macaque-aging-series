library(tidyverse)
library(ggplot2)
library(broom)
library(pheatmap)
library(reshape2)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(ggbeeswarm)
library(scales)

################################################################################ 
# AGE AND SEX DEG PLOTS 
################################################################################
cell_type_colors = c(
  "medium spiny neurons" = "#5e3c99",
  "cerebellar neurons" = "#44839f",
  "basket cells" = "#4e379a",
  "astrocytes" = "#d53a50",
  "oligodendrocytes" = "#fb8d64",
  "OPC" = "#fbae64",
  "vascular cells" = "#d1538f",
  "microglia" = "#bb471c",
  "ependymal cells" = "#fc6238",
  "glutamatergic neurons" = "#3d9995",
  "GABAergic neurons" = "#4d6aa8",
  "midbrain neurons" = "#65bfb0"
)

region_colors = c("NAc"= "#148e34","ACC" = "#5555b0","IPP" = "#6c9cc0","CN" = "#2caa4d",
                  "HIP" = "#84db55","lCb" = "#e6cd69", "mdTN" = "#117867", "dlPFC" = "#345296",
                  "M1" = "#467aa1", "EC" = "#532b70", "MB" = "#fc7462")

# Post-MASH DEG results of the additive model 
df = read_csv("/data/CEM/smacklab/data/bican/u01_results/mashr_additive_master_merged.csv")

################################################################################
# Fig 6A & 6B: sDEG LOESS model residual plots
################################################################################
sex_df = df %>% filter(sexM_lfsr < 0.05)
cell_obs = read.csv('/data/CEM/smacklab/data/bican/u01_results/u01_all_cells_metadata.csv') #dataframe with info for cluster size and UMI

# Adjust data structure to separate oligodendrocytes and OPCs 
cell_obs = cell_obs %>% 
  mutate(cell_class_annotation = case_when(
    # Change oligodendrocytes or OPCs in Louvain 12,13 to "OPC"
    (cell_class_annotation %in% c("oligodendrocytes", "oligodendrocyte precursor cells") & 
       louvain %in% c(12, 13)) ~ "OPC",
    # Change oligodendrocytes or OPCs in Louvain 0-11 to "oligodendrocytes"  
    (cell_class_annotation %in% c("oligodendrocytes", "oligodendrocyte precursor cells") & 
       louvain %in% 0:11) ~ "oligodendrocytes",
    TRUE ~ cell_class_annotation
  )) %>%
  mutate(id = paste0(cell_class_annotation, "_louvain", louvain, "_", region))

cell_obs = cell_obs%>%filter(age>3) # Remove infants 
cell_obs$celltype = cell_obs$cell_class_annotation

df_nebula = read.csv('/data/CEM/smacklab/data/bican/u01_results/ir_add_nebula_size_filt.csv') # Pre-MASH Nebula results, infants removed

# Adjust data structure to separate oligodendrocytes and OPCs 
df_nebula = df_nebula %>% mutate(cell_type = case_when(
  (cell_type == "oligodendrocytes" & louvain %in% c(12,13)) ~ "OPC",
  TRUE ~ cell_type)) %>% 
  mutate(id = paste0(cell_type,"_louvain",louvain,"_", region)) %>%
  mutate(gene.id = paste0(gene,"_",id))

summary = read_csv('/data/CEM/smacklab/data/bican/u01_results/mashr_deg_summary.csv') # Summary with n genes tested & sig DEG counts by subcluster-region 

df1 = df # create a backup df object to re-define df for downstream code
df = sex_df

###DATAFRAMES for info about number of tests, number of DEGs, UMI counts, and subcluster size
df_test = df_nebula%>%group_by(id)%>%summarise(test_count=n()) #number of test performed in NEBULA

df_deg =  df%>%group_by(cell_type,louvain,region)%>%summarise(n_DEG=n()) #raw number of DEGs

df_umi = cell_obs %>%
  mutate(celltype = cell_class_annotation) %>%
  group_by(celltype, louvain, region) %>%
  summarise(mean_umi = mean(n_umi))

df_size = cell_obs %>%
  mutate(celltype = cell_class_annotation) %>%
  group_by(celltype, louvain, region) %>%
  summarise(cluster_size = n())

# Edit column names to harmonize formatting
df_deg$celltype = df_deg$cell_type
df_test <- df_test %>%
  separate(id, 
           into = c("celltype", "louvain", "region"),
           sep = "_", remove=F) %>%
  mutate(louvain = as.numeric(gsub("louvain", "", louvain)))

df_combined <- Reduce(function(x, y) merge(x, y, by = c("celltype","louvain","region")), 
                      list(df_deg, df_umi, df_size,df_test))

df_combined$n_DEG_norm = df_combined$n_DEG/df_combined$test_count

df_combined = df_combined%>%mutate(log_mean_umi = log(mean_umi), log_cluster_size = log(cluster_size),
                                   log_nDEG_norm = log(n_DEG_norm*10000)) #log-transform


# Fit the LOESS model 
loess_model <- loess(log_nDEG_norm ~log_cluster_size+log_mean_umi, data = df_combined)
df_combined$residuals<- residuals(loess_model)

#write_csv(df_combined, "/data/CEM/smacklab/data/bican/u01_results/sDEG_loess_df_w_residuals.csv")

# Create abbreviations for each cell type 
df_combined = df_combined %>% mutate(ct_abbr = case_when(
  cell_type == "oligodendrocytes" ~ "OLIG",
  cell_type == "glutamatergic neurons" ~ "GLUT",
  cell_type == "microglia" ~ "MGL",
  cell_type == "astrocytes" ~ "AST",
  cell_type == "medium spiny neurons" ~ "MSN",
  cell_type == "GABAergic neurons" ~ "GABA",
  cell_type == "midbrain neurons" ~ "MBN",
  cell_type == "cerebellar neurons" ~ "CER",
  cell_type == "vascular cells" ~ "VASC",
  cell_type == "ependymal cells" ~ "EPEN",
  cell_type == "basket cells" ~ "BC",
  cell_type == "OPC" ~ "OPC"))

# Create the plots 
ggplot(df_combined, aes(x = reorder(ct_abbr, residuals, FUN = median, decreasing = F), 
                        y = residuals, fill = celltype)) +
  geom_boxplot(color="black") +
  scale_fill_manual(values = cell_type_colors) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylab("sDEG log fold enrichment") +
  theme(
    axis.text.y = element_text(size=12, color="black"),
   # axis.text.y = element_blank(),
    axis.text.x = element_text(size = 12, color="black"),
    axis.title.x = element_text(size=14, color="black"),
   axis.title.y = element_blank(),
    legend.position = "none") +
  coord_flip()

ggplot(df_combined, aes(x = reorder(region, residuals, FUN = median, decreasing = F), 
                        y = residuals, fill = region)) +
  geom_boxplot(color="black") +
  scale_fill_manual(values = region_colors) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylab("sDEG log fold enrichment") +
  theme(
    axis.text.y = element_text(size=12, color="black"),
    # axis.text.y = element_blank(),
    axis.text.x = element_text(size = 12, color="black"),
    axis.title.x = element_text(size=14, color="black"),
    axis.title.y = element_blank(),
    legend.position = "none") +
  coord_flip()
 # 6 h x 5 w
################################################################################
# Fig. 6C: sDEG directional binom result beeswarm plot by cell class 
################################################################################
df = read_csv("/data/CEM/smacklab/data/bican/u01_results/mashr_additive_master_merged.csv")
# Perform the binomal tests on sDEGs 
direction_summary <- df %>%
  mutate(sex_dir = sign(sexM_pm)) %>% # Male-biased = positive, female-biased = negative 
  group_by(id, louvain, region, cell_type) %>%
  filter(sexM_lfsr<0.05) %>% # significance threshold 
  summarize(male_fraction = sum(sex_dir>0)/n(),  
            binom_pval=binom.test(x=sum(sex_dir>0), # Test whether the male fraction is significantly different from 1:1, and in which direction
                                  n=n(), 
                                  alternative="two.sided")$p.value,
            conf_lower =binom.test(x=sum(sex_dir>0), 
                                   n=n(), 
                                   alternative="two.sided")$conf.int[1],
            conf_upper =binom.test(x=sum(sex_dir>0), 
                                   n=n(), 
                                   alternative="two.sided")$conf.int[2]) %>%
  mutate(p_adj = p.adjust(binom_pval, method="fdr"))

direction_summary = direction_summary %>%
  mutate(direction = case_when(
    male_fraction < 0.5 ~ "female biased",
    TRUE ~ "male biased"))

direction_summary = direction_summary %>% ungroup() %>% as_tibble()

direction_summary = direction_summary %>% mutate(ct_abbr = case_when(
  cell_type == "oligodendrocytes" ~ "OLIG",
  cell_type == "glutamatergic neurons" ~ "GLUT",
  cell_type == "microglia" ~ "MGL",
  cell_type == "astrocytes" ~ "AST",
  cell_type == "medium spiny neurons" ~ "MSN",
  cell_type == "GABAergic neurons" ~ "GABA",
  cell_type == "midbrain neurons" ~ "MBN",
  cell_type == "cerebellar neurons" ~ "CER",
  cell_type == "vascular cells" ~ "VASC",
  cell_type == "ependymal cells" ~ "EPEN",
  cell_type == "basket cells" ~ "BC",
  cell_type == "OPC" ~ "OPC"))


ggplot(direction_summary, aes(x=male_fraction,y=reorder(ct_abbr, male_fraction, FUN = median, decreasing = T), color=direction, alpha=p_adj<0.05))+
  geom_beeswarm()+
  geom_vline(xintercept = 0.5, lty="dashed")+
  scale_color_manual(values = c("female biased" = "green4", "male biased" = "purple3")) +
  scale_alpha_manual(values=c(0.30,1))+ #, name="FDR < 0.05")+
  theme_classic()+
  xlab("Proportion male-biased")+
  theme(
    axis.text.y= element_text(size=12, color="black"),
    axis.text.x = element_text(size=12, color="black"),
    axis.ticks.x = element_blank(),
    axis.title.y= element_blank(),
    axis.title.x = element_text(size=14, color="black"),
    legend.position = "none") # 8 w x 6 h
################################################################################
# Fig. 6D: a+sDEG hypergeometric test result box plot by cell class 
################################################################################
summary = read_csv('/data/CEM/smacklab/data/bican/u01_results/mashr_deg_summary.csv') # Summary with n genes tested & sig DEG counts by subcluster-region 
summary_10 = summary %>% filter(age_sig >= 10 | sex_sig >= 10)

test_both_DEG_enrichment <- function(data, 
                                     background_rate = NULL) {
  background_rate <- sum(data$both_sig) / sum(data$total_genes) # background rate for both_sig DEGs
  
  # Initialize results data frame
  results <- data.frame(
    id = data$id,
    observed_DEGs = data$both_sig,
    fold_enrichment = NA,
    pvalue = NA
  )
  
  # Perform hypergeometric test for each row
  for(i in 1:nrow(data)) {
    # Get counts for this subcluster
    x <- data$both_sig[i]        # Number of both_sig DEGs in this subcluster
    m <- data$age_sig[i]        # N aDEGs in this subcluster-region
    n <- data$total_genes[i] - m            # N genes that are not aDEGs in this subcluster-region
    k <- data$sex_sig[i]                # N sDEGs in this subcluster
    
    if(data$sex_sig[i] > 0 && data$age_sig[i] > 0) {
      p1 <- x/k
      p2 <- m/(m+n)
      
      results$fold_enrichment[i] <- p1/p2  # Fixed: added [i]
    } else {
      results$fold_enrichment[i] <- NA  # Handle division by zero
    }
    
    # Hypergeometric test (testing for over-representation)
    if(x > 0) {
      results$pvalue[i] <- 1 - phyper(x - 1, m, n, k)
    } else {
      results$pvalue[i] <- 1
    }
  }
  
  # Adjust p-values for multiple testing
  results$padj <- p.adjust(results$pvalue, method = "BH")
  
  return(results)
}
# Run the analysis
both_DEG_results<- test_both_DEG_enrichment(summary_10)
summary(both_DEG_results$fold_enrichment)


both_DEG_results <- transform(both_DEG_results,
                              cell_type = sub("(.*)_louvain\\d+_.*", "\\1", id),
                              louvain = as.numeric(sub(".*_louvain(\\d+)_.*", "\\1", id)),
                              region = sub(".*_louvain\\d+_(.*)", "\\1", id)) %>% na.omit()

both_DEG_results = both_DEG_results %>% mutate(log2FC = log2(fold_enrichment))
both_results_filt = both_DEG_results %>% filter(log2FC > -50) # Remove any -INF log2FC values for plotting

both_results_filt = both_results_filt %>% mutate(ct_abbr = case_when(
  cell_type == "oligodendrocytes" ~ "OLIG",
  cell_type == "glutamatergic neurons" ~ "GLUT",
  cell_type == "microglia" ~ "MGL",
  cell_type == "astrocytes" ~ "AST",
  cell_type == "medium spiny neurons" ~ "MSN",
  cell_type == "GABAergic neurons" ~ "GABA",
  cell_type == "midbrain neurons" ~ "MBN",
  cell_type == "cerebellar neurons" ~ "CER",
  cell_type == "vascular cells" ~ "VASC",
  cell_type == "ependymal cells" ~ "EPEN",
  cell_type == "basket cells" ~ "BC",
  cell_type == "OPC" ~ "OPC"))

ggplot(both_results_filt, aes(x=reorder(ct_abbr, log2(fold_enrichment), FUN=median),
                              y=log2(fold_enrichment), fill=cell_type)) + 
  geom_hline(yintercept = 0, lty="dashed") + 
  geom_boxplot(outlier.shape="") +
  #geom_point(aes(alpha=padj<0.05), shape=21, color="black", size=2, stroke=1) +
  coord_flip() + 
  scale_alpha_manual(values=c(0.2,1)) + 
  scale_fill_manual(values=cell_type_colors)+
  ylab("log2(fold enrichment) of overlap of aDEGs and sDEGs") + 
  xlab("") + 
  theme_classic() +
  theme(legend.position ="none",
        axis.text.x = element_text(size=12, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_text(size=14))
# 8 w x 6 h
################################################################################
# Fig. 6E: a+sDEG directional binom result beeswarm plot by subcluster-region 
################################################################################
df = read_csv("/data/CEM/smacklab/data/bican/u01_results/mashr_additive_master_merged.csv")

summary = read_csv("/data/CEM/smacklab/data/bican/u01_results/mashr_deg_summary.csv")
summary_10 = summary %>% filter(age_sig >= 10 | sex_sig >=10)
purple_10 = summary_10 %>% pull(id)
df_10 = df %>% filter(id %in% purple_10)

direction_summary = df_10 %>% 
  mutate(age_sex_dir = age_pm * sexM_pm) %>%
  group_by(id, louvain, region, cell_type) %>%
  filter(age_lfsr<0.05 & sexM_lfsr<0.05) %>% ##you can vary this
  summarize(male_older_fraction = sum(age_sex_dir>0)/n(),
            binom_pval=binom.test(x=sum(age_sex_dir>0), 
                                  n=n(), 
                                  alternative="two.sided")$p.value,
            conf_lower =binom.test(x=sum(age_sex_dir>0), 
                                   n=n(), 
                                   alternative="two.sided")$conf.int[1],
            conf_upper =binom.test(x=sum(age_sex_dir>0), 
                                   n=n(), 
                                   alternative="two.sided")$conf.int[2]) %>%
  mutate(p_adj = p.adjust(binom_pval, method="fdr"))

direction_summary = direction_summary %>%
  mutate(direction = case_when(
    male_older_fraction < 0.5 ~ "female older",
    TRUE ~ "male older"))

sig_result = direction_summary %>% filter(p_adj < 0.05)

direction_summary = direction_summary %>% filter(male_older_fraction >0 & male_older_fraction <1 & male_older_fraction !=0.5) %>% filter(cell_type %in% sig_result$cell_type) %>% filter(!cell_type == "ependymal cells")

direction_summary = direction_summary %>% mutate(ct_abbr = case_when(
  cell_type == "oligodendrocytes" ~ "OLIG",
  cell_type == "glutamatergic neurons" ~ "GLUT",
  cell_type == "microglia" ~ "MGL",
  cell_type == "astrocytes" ~ "AST",
  cell_type == "medium spiny neurons" ~ "MSN",
  cell_type == "GABAergic neurons" ~ "GABA",
  cell_type == "midbrain neurons" ~ "MBN",
  cell_type == "cerebellar neurons" ~ "CER",
  cell_type == "vascular cells" ~ "VASC",
  cell_type == "ependymal cells" ~ "EPEN",
  cell_type == "basket cells" ~ "BC",
  cell_type == "OPC" ~ "OPC"))

ggplot(direction_summary, aes(x= male_older_fraction, y=reorder(ct_abbr, male_older_fraction, FUN = median, decreasing = T), color=direction, alpha=p_adj<0.05))+
  geom_beeswarm(cex=1.25)+
  geom_vline(xintercept = 0.5, lty="dashed")+
  scale_color_manual(values = c("male older" = "firebrick3", "female older" = "navy")) +
  scale_alpha_manual(values=c(0.30,1), name="FDR < 0.05")+
  theme_classic()+
  labs(x = "Proportion male-older a+sDEGs") + 
  theme(
    axis.text.y= element_text(size=12, color="black"),
    axis.text.x = element_text(size=12, color="black"),
    axis.ticks = element_blank(),
    axis.title.y= element_blank(),
    axis.title.x = element_text(size=14, color="black"),
    legend.position = "none") # 8 w x 6 h

# 8 x 8 (because of legend)
################################################################################
# Fig. 6F: Heatmap of lipid a+sDEG ORA results in quadrant 1 
################################################################################
## See Fig_6_ORA_analysis.R for how ora_result was generated 

ora_result = read_csv("/scratch/klwatki5/R_analysis/bbi/mash_new/purple_genes/ora_results/louvain_region/db_update/ora_result_merged.csv") %>% filter(p.adjust < .05) %>% filter(Count > 1) %>% mutate(logfold_enrich = log10(FoldEnrichment))

ora_result = ora_result %>% mutate(id = paste0(cell_type,"_louvain",louvain,"_",region))

# Filter based on binom test 
sig_dir = read_csv("/scratch/klwatki5/R_analysis/bbi/mash_new/purple_genes/directional_binom_results_enrich_filt.csv") # 216
sig_dir$id = paste0(sig_dir$cell_type,"_louvain", sig_dir$louvain,"_", sig_dir$region)
ora_result = ora_result %>% filter(id %in% sig_dir$id) # 8907

# Define functional categories 
source(file="/scratch/klwatki5/R_analysis/bbi/gene_set_category_master_list_v2.R") # This file contains gene set names from our results, sorted into categories with help from Claude AI

# Assign the categories to the dataframe 
find_pathway_category <- function(pathway_id, category_list) {
  # Check each category
  for (category_name in names(category_list)) {
    if (pathway_id %in% category_list[[category_name]]) {
      return(category_name)
    }
  }
  # If no match found, return "Other"
  return("Other")
}

# Apply the function to add the pathway_category column to ora_result
ora_result <- ora_result %>%
  mutate(pathway_category = sapply(ID, find_pathway_category, category_list = pathway_categories))

# Convert underscores to spaces and capitalize for better readability in plots
ora_result <- ora_result %>%
  mutate(pathway_category_label = gsub("_", " ", pathway_category) %>%
           stringr::str_to_title())

# Add cell type annotations
subtype_labels = read.csv("/scratch/klwatki5/R_analysis/bbi/all_subtypes.csv")
subtype_labels$Cell.class <- gsub("OPC/oligodendrocytes", "oligodendrocytes", subtype_labels$Cell.class)
annotation_map <- subtype_labels[, c("Cell.class", "Louvain", "Putative.annotation")]

# Rename columns to match sig_result for easier merging
colnames(annotation_map) <- c("cell_type", "louvain", "annotation")

# Now merge with sig_result
ora_annotated <- merge(
  ora_result,
  annotation_map,
  by = c("cell_type", "louvain"),
  all.x = TRUE)

ora_annotated = ora_annotated %>% mutate(anno_id = paste(id, annotation, sep = " "))
ora_result = ora_annotated %>% filter(quadrant == 1)

ora_result = ora_result %>% filter(pathway_category =="Lipid_Metabolism_and_Cholesterol")

cell_type_order <- c("astrocytes", "ependymal cells", "microglia", "oligodendrocytes", "OPC",
                     "vascular cells", "basket cells", "cerebellar neurons", 
                     "GABAergic neurons", "glutamatergic neurons", 
                     "medium spiny neurons", "midbrain neurons")

cell_type_palette <- c(
  "medium spiny neurons" = "#5e3c99",
  "cerebellar neurons" = "#44839f",
  "basket cells" = "#4e379a",
  "astrocytes" = "#d53a50",
  "oligodendrocytes" = "#fb8d64",
  "OPC" = "#fbae64",
  "vascular cells" = "#d1538f",
  "microglia" = "#bb471c",
  "ependymal cells" = "#fc6238",
  "glutamatergic neurons" = "#3d9995",
  "GABAergic neurons" = "#4d6aa8",
  "midbrain neurons" = "#65bfb0"
)

# Create region palette based on extracted regions
region_palette = c("NAc"= "#148e34","ACC" = "#5555b0","IPP" = "#6c9cc0","CN" = "#2caa4d",
                   "HIP" = "#84db55","lCb" = "#e6cd69", "mdTN" = "#117867", "dlPFC" = "#345296",
                   "M1" = "#467aa1", "EC" = "#532b70", "MB" = "#fc7462")

#Spread data to create a matrix
heatmap_data <- ora_result %>%
  dplyr::select(ID, anno_id, logfold_enrich) %>%
  pivot_wider(names_from = anno_id, values_from = logfold_enrich, values_fill = 0) %>%
  as.data.frame()

# heatmap_data = matrix_data # if there were duplicates, start here 
heatmap_data <- as.data.frame(heatmap_data)
rownames(heatmap_data) <- heatmap_data$ID
# Remove ID column
heatmap_data$ID <- NULL

# Get unique pathway categories, cell types, and regions
unique_categories <- unique(ora_result$pathway_category_label)
unique_cell_types <- unique(ora_result$cell_type)
unique_regions <- unique(ora_result$region)

# Create row annotation data
row_annotation_data <- ora_result %>%
  dplyr::select(ID, pathway_category_label) %>%
  distinct() %>%
  arrange(pathway_category_label)

# Set row names
rownames(row_annotation_data) <- row_annotation_data$ID

# Create row annotation dataframe for pheatmap
row_anno_df <- data.frame(
  Pathway_Category = row_annotation_data$pathway_category_label,
  row.names = row_annotation_data$ID
)

# Create the col_data frame
col_data <- data.frame(
  anno_id = colnames(heatmap_data),
  stringsAsFactors = FALSE
)

# Process each anno_id individually with a for loop for clarity
cell_types <- vector("character", length(col_data$anno_id))
louvains <- vector("character", length(col_data$anno_id))
regions <- vector("character", length(col_data$anno_id))
annotations <- vector("character", length(col_data$anno_id))

for (i in 1:length(col_data$anno_id)) {
  id <- col_data$anno_id[i]
  
  # Find position of first space after the last underscore
  last_under <- max(gregexpr("_", id)[[1]])
  
  # Look for a space after the last underscore
  temp_str <- substr(id, last_under+1, nchar(id))
  space_pos <- regexpr(" ", temp_str)
  
  # If space found, split the string accordingly
  if (space_pos > 0) {
    # Position in the original string
    actual_space_pos <- last_under + space_pos
    
    # Split into ID part and annotation part
    id_part <- substr(id, 1, actual_space_pos-1)
    annotations[i] <- substr(id, actual_space_pos+1, nchar(id))
    
    # Split the ID part by underscores
    id_pieces <- strsplit(id_part, "_")[[1]]
    
    # Extract the components
    cell_types[i] <- id_pieces[1]
    louvains[i] <- id_pieces[2]
    regions[i] <- id_pieces[3]
  } else {
    # If no space found, just split by underscores
    id_pieces <- strsplit(id, "_")[[1]]
    
    # Extract what we can
    cell_types[i] <- id_pieces[1]
    if (length(id_pieces) > 1) louvains[i] <- id_pieces[2]
    if (length(id_pieces) > 2) regions[i] <- id_pieces[3]
    annotations[i] <- NA
  }
}

# Add the extracted components to col_data
col_data$cell_type <- cell_types
col_data$louvain <- louvains
col_data$region <- regions
col_data$annotation <- annotations

# Create column annotation dataframe for pheatmap
col_anno_df <- data.frame(
  Region = col_data$region,
  Cell_Type = col_data$cell_type,
  row.names = col_data$anno_id
)

# Use bolder, more contrasting colors for pathway categories
category_palette <- c()

# Ensure all categories have colors assigned
missing_categories <- setdiff(unique_categories, names(category_palette))
if(length(missing_categories) > 0) {
  additional_colors <- rainbow(length(missing_categories))
  names(additional_colors) <- missing_categories
  category_palette <- c(category_palette, additional_colors)
}

# Compile all color palettes
annotation_colors <- list(
  Pathway_Category = category_palette,
  Cell_Type = cell_type_palette,
  Region = region_palette
)

# Order rows by pathway category (keep this as before)
ordered_rows <- row_annotation_data$ID
ordered_categories <- row_anno_df[ordered_rows, "Pathway_Category"]
category_changes <- which(ordered_categories[-1] != ordered_categories[-length(ordered_categories)])

# ORDER COLUMNS BY CELL TYPE instead of region
col_data$cell_type_factor <- factor(col_data$cell_type, levels = cell_type_order)
col_data_sorted <- col_data %>% arrange(cell_type_factor)
ordered_cols <- col_data_sorted$anno_id

# Reorder the heatmap data with new column ordering
heatmap_data_ordered <- heatmap_data[match(ordered_rows, rownames(heatmap_data)), 
                                     match(ordered_cols, colnames(heatmap_data))]

# Find where cell types change for creating gaps
ordered_cell_types <- col_data_sorted$cell_type
cell_type_changes <- which(ordered_cell_types[-1] != ordered_cell_types[-length(ordered_cell_types)])

# Create the heatmap with cell type ordering and gaps
pheatmap(
  heatmap_data_ordered,
  color = colorRampPalette(c("navy", "gray90", "firebrick3"))(100),
  breaks = seq(-2, 2, length.out = 101),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = row_anno_df,
  annotation_col = col_anno_df[ordered_cols, ],  # Reorder annotation_col to match new ordering
  annotation_colors = annotation_colors,
  show_rownames = F,
  show_colnames = F,
  annotation_names_row = F,
  annotation_names_col = F,
  fontsize_row = 8,
  fontsize_col = 12,
  gaps_row = category_changes,
  gaps_col = cell_type_changes,  # Now creates gaps between different cell types
  border_color = "NA"
) # 12 x 10 landscape 

################################################################################
# Fig 6F: a+sDEG and AD-DEG overlap
################################################################################
## See Fig_6_AD_DEG.R for how ora_result was generated 

ora_result = read_csv('/scratch/klwatki5/R_analysis/bbi/mash_new/ad_analysis/ora_results/purple_genes/ct_match_mathys/all_merged.csv')

# Prep data for plotting
ora_result = ora_result %>% 
  mutate(log2OR = log2(odds_ratio),
         log2OR_upperci = log2(odds_ratio_upperci),
         log2OR_lowerci = log2(odds_ratio_lowerci)) %>%
  filter(p_adj < 0.05) %>% 
  #mutate(log2OR = pmax(pmin(log2OR, 3), -3)) %>%
  mutate(Term_clean = gsub("_(up|down)$", "", Term)) %>%
  mutate(Term_clean = gsub("_AD", "", Term_clean))

# Create the plot
ggplot(ora_result, aes(x = log2OR, 
                       y = paste(Term, region),
                       color = factor(quadrant))) +
  geom_point(size = 3, alpha = 0.8) +
  geom_linerange(aes(xmin = log2OR_lowerci, xmax = log2OR_upperci),
                 alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("1" = "red4", "2" = "lightblue3", "3" = "pink3", "4" = "blue4"),name = "Quadrant") +
  theme_classic() +
  labs(x = expression(log[2]~odds~ratio)) +
  theme(#axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_blank(),
        legend.position = "none")

