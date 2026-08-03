library(dplyr)
library(tidyr)
library(tidyverse)
library(org.Mmu.eg.db)
library(AnnotationDbi)
library(ggplot2)
source("utils/annotation.R")

all_df = read.csv('/scratch/wyang114/outputs/transcriptional/GO_terms/gsea_C5_summary.csv')
load('/scratch/wyang114/outputs/transcriptional/GO_terms/gsea_go_cluster_new.RData')
gene_group = col_annotation%>%filter(Cluster %in% c(8,9))
df_all = read.csv('/scratch/nsnyderm/u01/intermediate_files/nebula/age/all_mashr_summary_filter_missing.csv')[-1]

label_df <- stack(glut_labels)
colnames(label_df) <- c("louvain", "supergroup")

# Preprocess: match descriptions to GO IDs
all_terms <- all_df %>%
  dplyr::select(ID, Description) %>%
  distinct() %>%
  column_to_rownames("Description")

# Loop through all pathways
ct='glutamatergic neurons'
results_list <- list()

for (pathway_desc in row.names(gene_group)) {
  
  # Step 1: Get GO ID
  go_id <- all_terms[pathway_desc, , drop = TRUE]
  
  # Step 2: Get core enriched genes (SYMBOL/ENTREZID) for that pathway
  core_sets <- all_df %>%
    filter(cell_type == "glutamatergic neurons", Description == pathway_desc) %>%
    separate_rows(core_enrichment, sep = "/") %>%
    mutate(core_enrichment = trimws(core_enrichment)) %>%
    filter(core_enrichment != "") %>%
    pull(core_enrichment)
  
  # Step 3: Get groups where the pathway is enriched
  enriched_groups <- all_df %>%
    filter(Description == pathway_desc) %>%
    mutate(id = paste(region, cell_type, louvain, sep = "_")) %>%
    pull(id)
  
  # Step 4: Map GO → ENTREZID → ENSEMBL/SYMBOL
  entrez_ids <- mapIds(org.Mmu.eg.db, keys = go_id, column = "ENTREZID", keytype = "GOALL", multiVals = "list")[[1]]
  ensembl_map <- AnnotationDbi::select(org.Mmu.eg.db, keys = entrez_ids, columns = c("ENSEMBL", "SYMBOL"), keytype = "ENTREZID") %>%
    filter(ENTREZID %in% core_sets)
  
  # Step 5: Extract and filter long-format expression data
  df_long <- df_all %>%
    filter(celltype==ct)%>%
    filter(gene %in% ensembl_map$ENSEMBL) %>%
    mutate(id = paste(region, cell_type, louvain, sep = "_")) %>%
    filter(id %in% enriched_groups)%>%mutate(scaled_effect = age_mashr_pm/age_mashr_psd)
  
  # Step 6: Build matrix and filter for genes expressed in ≥80% of groups
  mat <- df_long %>%
    mutate(gene_id = paste(gene, external_gene_name, sep = "_")) %>%
    dplyr::select(id, gene_id, scaled_effect) %>%
    pivot_wider(id_cols = id, names_from = gene_id, values_from = scaled_effect) %>%
    column_to_rownames("id")
  
  select_col <- colnames(mat)[colSums(!is.na(mat)) > 0.8 * nrow(mat)]
  
  # Step 7: Compute normalized ranks per gene and run KW test on all combinations
  # df_ranked <- df_long %>%
  #   filter(paste(gene, external_gene_name, sep = "_") %in% select_col) %>%
  #   group_by(gene) %>%
  #   mutate(
  #     rank = rank(age_mashr_pm),
  #     norm_rank = (rank - 1) / (n() - 1)
  #   ) %>%
  #   ungroup()
  # kw_test <- kruskal.test(norm_rank ~ id, data = df_ranked)
  
  # # Step 8: Compute normalized ranks per gene and run KW test on region
  # df_ranked <- df_long %>%
  #   filter(paste(gene, external_gene_name, sep = "_") %in% select_col) %>%
  #   group_by(gene) %>%
  #   mutate(
  #     rank = rank(age_mashr_pm),
  #     norm_rank = (rank - 1) / (n() - 1)
  #   ) %>%
  #   ungroup()
  # kw_test_r <- df_ranked %>% group_by(region,louvain)%>%
  #   summarise(mean_rank = mean(norm_rank))%>%
  #   kruskal.test(mean_rank ~ region, data = .)
  
  # Step 9: Compute normalized ranks per gene and run KW test on subclass
  # df_ranked <- df_long %>%
  #   filter(paste(gene, external_gene_name, sep = "_") %in% select_col) %>%
  #   group_by(gene) %>%
  #   mutate(
  #     rank = rank(age_mashr_pm),
  #     norm_rank = (rank - 1) / (n() - 1)
  #   ) %>%
  #   ungroup()
  # kw_test_sc <- df_ranked %>% 
  #   group_by(region,louvain)%>%
  #   summarise(mean_rank = mean(norm_rank))%>%
  #   left_join(label_df, by = "louvain")%>%
  #   kruskal.test(mean_rank ~ subtype, data = .)
  
  #ALTERNATIVE: anova on multiple covariates
  df_long = df_long%>%left_join(label_df,by='louvain')
  aov_res_id <- summary(aov(scaled_effect ~ id, data = df_long))
  
  temp = df_long%>%
    group_by(id,region,supergroup)%>%summarise(mean_scaled_effect=mean(scaled_effect))
  aov_res <- summary(aov(mean_scaled_effect ~ region + supergroup, data = temp))
  aov_res_2 <- summary(aov(scaled_effect ~ region + supergroup, data = df_long))
  
  # # Optional: store results
  # results_list[[pathway_desc]] <- list(
  #   pathway = pathway_desc,
  #   go_id = go_id,
  #   p_value = kw_test$p.value,
  #   kw_statistic = kw_test$statistic,
  #   df = kw_test$parameter,
  #   p_value_r = kw_test_r$p.value,
  #   kw_region_statistic = kw_test_r$statistic,
  #   df_r = kw_test_r$parameter,
  #   p_value_sc = kw_test_sc$p.value,
  #   kw_subclass_statistic = kw_test_sc$statistic,
  #   df_sc = kw_test_sc$parameter,
  #   n_genes = length(select_col)
  # )
  results_list[[pathway_desc]] <- list(
    pathway = pathway_desc,
    go_id = go_id,
    p_value_id = aov_res_id[[1]]["id", "Pr(>F)"],
    df_id = aov_res_id[[1]]["id", "Df"],
    p_value_r = aov_res[[1]]["region", "Pr(>F)"],
    df_r = aov_res[[1]]["region", "Df"],
    p_value_sc = aov_res[[1]]["supergroup", "Pr(>F)"],
    df_sc = aov_res[[1]]["supergroup", "Df"],
    p_value_r_2 = aov_res_2[[1]]["region", "Pr(>F)"],
    df_r_2 = aov_res_2[[1]]["region", "Df"],
    p_value_sc_2 = aov_res_2[[1]]["supergroup", "Pr(>F)"],
    df_sc_2 = aov_res_2[[1]]["supergroup", "Df"],
    n_genes = length(select_col)
  )
  
  
  # #Optional: plot (comment out if looping over many pathways)
  # print(
  #   df_ranked %>%
  #     ggplot(aes(x = id, y = norm_rank)) +
  #     geom_boxplot() +
  #     labs(title = pathway_desc, subtitle = paste0("KW p = ", signif(kw_test$p.value, 3))) +
  #     theme_minimal() +
  #     theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # )
}

# clean_result <- function(x) {
#   data.frame(
#     pathway = x$pathway,
#     go_id = x$go_id,
#     p_value = x$p_value,
#     kw_statistic = as.numeric(x$kw_statistic),
#     df = as.integer(x$df),
#     p_value_r = x$p_value_r,
#     kw_region_statistic = as.numeric(x$kw_region_statistic),
#     df_r = as.integer(x$df_r),
#     p_value_sc = x$p_value_sc,
#     kw_subclass_statistic = as.numeric(x$kw_subclass_statistic),
#     df_sc = ifelse(is.null(x$df_sc), NA, as.integer(x$df_sc)),
#     n_genes = x$n_genes,
#     stringsAsFactors = FALSE
#   )
# }
# 
# results_summary <- do.call(rbind, lapply(results_list, clean_result))

results_summary <- do.call(rbind, lapply(results_list, data.frame))
results_summary$padj = p.adjust(results_summary$p_value_id, method = "BH")
results_summary$padj_r = p.adjust(results_summary$p_value_r_2, method = "BH")
results_summary$padj_sc = p.adjust(results_summary$p_value_sc_2, method = "BH")
hist(results_summary$p_value_id,breaks=50)
hist(results_summary$p_value_r_2,breaks=50)
hist(results_summary$p_value_sc_2,breaks=50)
