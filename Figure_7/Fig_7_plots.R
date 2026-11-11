library(tidyverse)
library(ggplot2)
library(broom)
library(pheatmap)
library(ggbeeswarm)

################################################################################
# SEX*AGE INTERACTION DEG PLOTS 
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

# Leaving for reference, re-wrote results file with this new column 
# df = df %>% mutate(category = case_when(
#   (interaction_lfsr < 0.05 & quadrant == 1) ~ "male stronger pos",
#   (interaction_lfsr < 0.05 & quadrant == 2) ~ "male stronger neg",
#   (interaction_lfsr < 0.05 & quadrant == 3) ~ "female stronger neg",
#   (interaction_lfsr < 0.05 & quadrant == 4) ~ "female stronger pos",
#   interaction_lfsr > 0.05 ~ "no interaction"
# ))

# Post-MASH DEG results of the additive model 
df = read_csv("/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/4_merged_dataframes/master_merged_new_quadrants.csv")

summary = read_csv("/scratch/klwatki5/R_analysis/bbi/mash_new/interaction/results/4_merged_dataframes/interaction_summary_new_quadrants.csv")

################################################################################
# 7A & 7B: Enrichment of a*sDEGs binom tests
################################################################################
test_interaction_deg_enrichment <- function(df) {
  # Calculate the background proportion of interaction DEGs
  background_prop <- sum(summary$interaction_sig) / sum(summary$total_genes)

  # Function to perform the binomial test for a single row
  test_row <- function(total_genes, interaction_sig) {
    # Skip if data is missing
    if(is.na(total_genes) || is.na(interaction_sig) || total_genes == 0) {
      return(tibble(
        p_value = NA,
        ci_lower = NA,
        ci_upper = NA,
        enrichment = NA,
        enrichment_ci_lower = NA,
        enrichment_ci_upper = NA
      ))
    }

    # Calculate observed proportion
    observed_prop <- interaction_sig / total_genes

    # Calculate enrichment
    enrichment <- observed_prop / background_prop

    # Perform binomial test
    binom_result <- binom.test(x = interaction_sig,
                               n = total_genes,
                               p = background_prop,
                               alternative = "greater")

    # Return results
    return(tibble(
      p_value = binom_result$p.value,
      ci_lower = binom_result$conf.int[1],
      ci_upper = binom_result$conf.int[2],
      enrichment = enrichment,
      enrichment_ci_lower = binom_result$conf.int[1] / background_prop,
      enrichment_ci_upper = binom_result$conf.int[2] / background_prop
    ))
  }

  # Apply the test to each row and bind results
  result_list <- map2(df$total_genes, df$interaction_sig, test_row)
  results <- bind_rows(result_list)

  # Bind to original data
  return(bind_cols(df, results))
}

# Run the binomial tests
binom_result <- test_interaction_deg_enrichment(summary)
# Adjust p-values
binom_result$p_adj <- p.adjust(binom_result$p_value, method = "BH")

binom_sig = binom_result %>% filter(p_adj < 0.05) %>% pull(id)

# Create a summary table
summary_table <- binom_result %>%
  dplyr::select(id, louvain, region, cell_type,
                total_genes, interaction_sig, enrichment,
                p_value, p_adj) %>%
  mutate(log_enrichment = log10(enrichment))

summary_table = summary_table %>% filter(log_enrichment > -10) # remove negative infinity values
background_prop <- sum(summary$interaction_sig) / sum(summary$total_genes)

summary = summary %>% mutate(sig_enriched = case_when(id %in% binom_sig ~ "TRUE",
                                                      TRUE ~ "FALSE")) %>% filter(prop_interaction > 1e-45)

summary = summary %>% group_by(cell_type) %>% mutate(max_ct_prop = max(prop_interaction)) %>% ungroup()
summary = summary %>% group_by(region) %>% mutate(max_reg_prop = max(prop_interaction)) %>% ungroup()

ct_prop_summary = summary %>% group_by(cell_type) %>% summarise(mean(prop_interaction)) %>% mutate(perc = `mean(prop_interaction)`*100) 
reg_prop_summary = summary %>% group_by(region) %>% summarise(mean(prop_interaction)) %>% mutate(perc = `mean(prop_interaction)`*100) 

summary = summary %>% mutate(ct_abbr = case_when(
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

ggplot(summary, aes(y=prop_interaction, x=reorder(ct_abbr, max_ct_prop, decreasing= F), fill=cell_type)) +
  geom_boxplot(outlier.shape = "") +
  geom_point()+
  scale_fill_manual(values = cell_type_colors) +
  labs(
    x= "",
    y="proportion a*sDEGs") +
  theme_classic() +
  theme(axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=12, color="black"),
        axis.title.x = element_text(size=14, color="black"),
        legend.position = "none")+
  coord_flip() # 10 x 5


ggplot(summary, aes(y=prop_interaction, x=reorder(region, max_ct_prop, decreasing= F), fill=region)) +
  geom_boxplot(outlier.shape = "") +
  geom_point()+
  #geom_hline(yintercept = background_prop, lty="dashed", col="red") +
  scale_fill_manual(values = region_colors) +
  labs(
    x= "",
    y="proportion a*sDEGs") +
  theme_classic() +
  theme(axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=12, color="black"),
        axis.title.x = element_text(size=14, color="black"),
        legend.position = "none")+
  coord_flip() # 10 x 5


################################################################################
# 7C: Summary of interaction effects by effect direction
################################################################################
sig_int = df %>% filter(interaction_lfsr < 0.05)

sig_counts <- sig_int %>% 
  group_by(stronger_direction, stronger_in) %>% 
  summarise(sig_counts = n(), .groups = 'drop')

ggplot(sig_counts, aes(x = stronger_direction, y = sig_counts, fill = stronger_in)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("female" = "deepskyblue", "male" = "deepskyblue4")) +
  labs(x="",
       y="")+
  theme_classic()+
  theme(axis.text.y = element_text(size=12, color="black"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position= "none")

################################################################################
# 7D: IQCK gene expression in astrocytes 
################################################################################

### TO RE-MAKE IQCK PLOT, READ IN RAW FILTERED FILE BELOW COMMENTED CODE! 

## Commented-out code shows steps to create the gene expression plot for male vs female across ages. Can be done for any gene in our data. 

# var = readr::read_csv('/scratch/klwatki5/R_analysis/bbi/all_var_meta.csv') # Get ENSEMBL gene name 

# library(ggplot2)
# library(dplyr)
# library(tidyr)
# library(reticulate)
# use_condaenv(condaenv = '~/.conda/envs/TACCO_env/', required = TRUE)
# 
# sc <- import("scanpy")
# pd <- import("pandas")
# np <- import("numpy")

# Load astrocyte anndata file
# adata <- sc$read("/data/CEM/smacklab/data/bican/u01_cell_class/astrocytes.h5ad")
# ct = "astrocytes"

# # Define gene(s) of interest using Ensembl ID
# genes_of_interest <- c("ENSMMUG00000002708") # IQCK

### PERFORM NORMALIZATION AND UMI REGRESSION IN PYTHON
# 
# # Store raw data first
# adata$raw <- adata
# # Normalize to counts per million
# sc$pp$normalize_total(adata, target_sum = 1e6L)
# # Log transform
# sc$pp$log1p(adata)
# # Regress out UMI (assuming 'n_umi' is the column name, adjust if different)
# sc$pp$regress_out(adata, keys = 'n_umi')
# 
# # Updated function to extract normalized expression data
# extract_expression_data_normalized <- function(adata, genes, group_cols = c("age", "sex", "region", "louvain"), 
#                                                use_layer = NULL, use_raw = FALSE) {
#   
#   # Assign to Python environment
#   py$adata <- adata
#   
#   # Check what data layers are available
#   py_run_string("
# import numpy as np
# 
# # Check available layers
# if adata.layers:
#     available_layers = list(adata.layers.keys())
#     print('Available layers:', available_layers)
# else:
#     available_layers = []
#     print('No layers found')
# 
# # Check if raw data is available
# has_raw = adata.raw is not None
# print('Has raw data:', has_raw)
# 
# # Check the main X matrix info
# print('Main X matrix shape:', adata.X.shape)
# print('Main X matrix type:', type(adata.X))
# print('Main X matrix contains normalized/log-transformed data after processing')
# ")
#   
#   available_layers <- py$available_layers
#   has_raw <- py$has_raw
#   
#   cat("Available layers:", paste(available_layers, collapse = ", "), "\n")
#   cat("Has raw data:", has_raw, "\n")
#   
#   # Get available genes
#   py_run_string("available_genes = adata.var_names.tolist()")
#   available_genes <- py$available_genes
#   
#   genes_to_use <- intersect(genes, available_genes)
#   
#   cat("Extracting processed data for genes:", paste(genes_to_use, collapse = ", "), "\n")
#   
#   # Determine which expression matrix to use
#   if (!is.null(use_layer)) {
#     if (use_layer %in% available_layers) {
#       cat("Using layer:", use_layer, "\n")
#       expr_source <- paste0("adata.layers['", use_layer, "']")
#       data_description <- paste("layer:", use_layer)
#     } else {
#       stop(paste("Layer", use_layer, "not found. Available layers:", paste(available_layers, collapse = ", ")))
#     }
#   } else if (use_raw && has_raw) {
#     cat("Using raw data (original counts)\n")
#     expr_source <- "adata.raw.X"
#     data_description <- "raw counts"
#   } else {
#     cat("Using main X matrix (normalized, log-transformed, UMI-regressed)\n")
#     expr_source <- "adata.X"
#     data_description <- "normalized, log-transformed, UMI-regressed"
#   }
#   
#   # Extract expression data gene by gene
#   all_expr_data <- list()
#   
#   for (i in seq_along(genes_to_use)) {
#     gene <- genes_to_use[i]
#     cat("Processing gene:", gene, "\n")
#     
#     # Assign single gene to Python
#     py$current_gene <- gene
#     py$expr_source_str <- expr_source
#     
#     # Get this gene's data from the specified source
#     py_run_string(paste0("
# import numpy as np
# 
# # Get gene index
# gene_idx = adata.var_names.get_loc(current_gene)
# 
# # Extract expression for this gene from specified source
# expr_matrix = ", expr_source, "
# if hasattr(expr_matrix, 'toarray'):
#     expr_vector = expr_matrix[:, gene_idx].toarray().flatten()
# else:
#     expr_vector = expr_matrix[:, gene_idx].flatten()
# 
# # Convert to regular Python array if it's a matrix
# if hasattr(expr_vector, 'A1'):
#     expr_vector = expr_vector.A1
# "))
#     
#     # Get the expression vector
#     expr_vector <- py$expr_vector
#     all_expr_data[[gene]] <- expr_vector
#     
#     # Print some stats about this gene's expression
#     cat("  Expression range:", range(expr_vector), "\n")
#     cat("  Non-zero cells:", sum(expr_vector > 0), "out of", length(expr_vector), "\n")
#   }
#   
#   # Combine into matrix
#   expr_matrix <- do.call(cbind, all_expr_data)
#   colnames(expr_matrix) <- genes_to_use
#   
#   # Get cell metadata
#   obs_data <- py_to_r(adata$obs)
#   
#   # Check if required columns exist
#   missing_cols <- setdiff(group_cols, colnames(obs_data))
#   if (length(missing_cols) > 0) {
#     stop(paste("These columns are missing from cell metadata:", paste(missing_cols, collapse = ", ")))
#   }
#   
#   # Combine expression and metadata
#   expr_df <- as.data.frame(expr_matrix)
#   
#   # Add metadata columns
#   for (col in group_cols) {
#     expr_df[[col]] <- obs_data[[col]]
#   }
#   
#   # Convert to long format
#   long_df <- expr_df %>%
#     pivot_longer(
#       cols = all_of(genes_to_use),
#       names_to = "gene",
#       values_to = "expression"
#     )
#   
#   # Calculate summary statistics
#   summary_df <- long_df %>%
#     group_by(across(all_of(c(group_cols, "gene")))) %>%
#     summarise(
#       mean_expr = mean(expression, na.rm = TRUE),
#       median_expr = median(expression, na.rm = TRUE),
#       sd_expr = sd(expression, na.rm = TRUE),
#       sem_expr = sd(expression, na.rm = TRUE) / sqrt(n()),
#       n_cells = n(),
#       .groups = "drop"
#     )
#   
#   return(list(
#     summary_data = summary_df,
#     raw_data = long_df,
#     available_genes = genes_to_use,
#     data_source = data_description
#   ))
# }
# 
# # Extract processed data (normalized, log-transformed, UMI-regressed)
# result_processed <- extract_expression_data_normalized(adata, genes_of_interest, group_cols = c("age", "sex", "region", "louvain"))
# summary_data_processed <- result_processed$summary_data
# 
# # summary_data_processed = summary_data_processed %>% filter(!region == "HIP") # remove HIP cells from medium spiny neuron plots
# # sig_int = int %>% filter(interaction_lfsr < 0.05)
# sig_clusters = sig_int %>% filter(ensembl_gene_id == genes_of_interest & cell_type == ct) %>% pull(id)
# summary_data_processed$id = paste0(ct,"_louvain",summary_data_processed$louvain,"_",summary_data_processed$region)
# 
# summary_filtered = summary_data_processed %>% filter(id %in% sig_clusters)
# 
# raw_filtered = result_processed$raw_data 
# raw_filtered$id = paste0(ct,"_louvain",raw_filtered$louvain,"_",raw_filtered$region)
# raw_filtered = raw_filtered %>% filter(id %in% sig_clusters) %>% filter(age>1)
# 
# raw_filtered = raw_filtered %>% group_by(louvain,region) %>% mutate(zscore_expr=scale(expression)) %>% ungroup() 
# raw_filtered$zscore_expr <- as.vector(raw_filtered$zscore_expr)


raw_filtered = read_csv('/scratch/klwatki5/R_analysis/bbi/mash_new/iqck_expression_zscored.csv')

iqck = raw_filtered
iqck %>% 
  group_by(age,sex,id) %>% summarize(mean_z=mean(zscore_expr), count=n()) %>% ungroup() %>%
  filter(count>10) %>%
  ggplot(aes(x=age,y=mean_z,col=sex, fill=sex)) + 
  geom_point(alpha=0.6)+
  geom_smooth(method="lm", alpha=0.2)+ 
  scale_color_manual(values = c("F" = "green4", "M" = "purple3")) +
  scale_fill_manual(values = c("F" = "green4", "M" = "purple3")) + 
  theme_classic()+
  labs(y="scaled expression", 
       x="age (years)") + 
  theme(legend.position = "none",
        axis.text = element_text(size=12, color="black"),
        axis.title = element_text(size=14)) # 8 x 6
