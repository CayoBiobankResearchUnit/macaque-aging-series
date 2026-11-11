library(tidyverse)
library(clusterProfiler)
library(org.Mmu.eg.db)
library(msigdbr)

df = read_csv('/data/CEM/smacklab/data/bican/u01_results/mashr_additive_master_merged.csv')

rm(list = setdiff(ls(), "df")) # Remove all but original dataframe to run for up vs down directions 

# Set age-and-AD direction for analysis; run up and down separately 
dir = "down" 

################################################################################
# Prepare AD genes
################################################################################
ref_df = read.csv('/scratch/klwatki5/R_analysis/bbi/mash_new/ad_analysis/sc_all_deg_summary.csv', row.names=1) # External data containing AD DEGs 

ref_df = ref_df %>% mutate(ref_ct_match = case_when(
  ct_label %in% c("Ast", "ASC") ~ "astrocytes",
  ct_label %in% c("Exc","EX") ~ "glutamatergic neurons",
  ct_label %in% c("Inh", "INH") ~ "GABAergic neurons",
  ct_label %in% c("Immune_Mic", "MG") ~ "microglia",
  ct_label %in% c("Oli","olig","ODC") ~ "oligodendrocytes",
  ct_label %in% c("Opc") ~ "OPC",
  ct_label %in% c("Epithelia_End","Epithelia_Per","Epithelia_SMC","endo","PER.END") ~ "vascular cells",
  TRUE ~ ct_label))


ref_df = ref_df %>% mutate(ref_id = paste0(ref_ct_match,"_AD")) %>%
  mutate(direction = case_when(
    LFC > 0 ~ "up",
    LFC < 0 ~ "down"
  )) %>%
  mutate(term = paste(ref_id, direction, sep = "_")) %>%
  filter(dataset == "Mathys_2024")


# Run separately for up and down genes
ref_df = ref_df %>% filter(direction == dir)
ref_df_selected = ref_df

# Convert AD gene symbols to Entrez IDs
ad_genes_symbols <- unique(ref_df_selected$gene)
symbol_to_entrez <- AnnotationDbi::select(org.Mmu.eg.db, 
                                          keys = ad_genes_symbols,
                                          columns = c("ENTREZID", "SYMBOL"),
                                          keytype = "SYMBOL")

# Handle duplicates - keep first mapping only
symbol_to_entrez_unique <- symbol_to_entrez %>%
  group_by(SYMBOL) %>%
  dplyr::slice(1) %>%
  ungroup()

# Merge Entrez IDs back to ref_df - KEEP ALL ORIGINAL COLUMNS
ref_df_entrez <- ref_df_selected %>%
  left_join(symbol_to_entrez_unique, by = c("gene" = "SYMBOL")) %>%
  filter(!is.na(ENTREZID))

print(paste("AD genes with Entrez IDs:", length(unique(ref_df_entrez$ENTREZID))))
print(paste("Number of rows in ref_df_entrez:", nrow(ref_df_entrez)))

create_ref_gene_set_groups <- function(ref_df_entrez) {
  gene_set_groups <- list()
  
  # Get unique terms (which already include direction)
  unique_terms <- unique(ref_df_entrez$term)
  
  print(paste("Creating", length(unique_terms), "gene set groups:"))
  
  for(term_name in unique_terms) {
    # Filter genes for this specific term
    filtered_data <- ref_df_entrez %>%
      filter(term == term_name) %>%
      dplyr::select(term, ENTREZID) %>%
      distinct()
    
    # Rename columns for consistency with downstream analysis
    names(filtered_data) <- c("gs_name", "ncbi_gene")
    
    # Extract metadata from the term for easier access later
    # Parse the term to get dataset, cell_type, region, and direction
    term_parts <- str_split(term_name, "_")[[1]]
    direction <- term_parts[length(term_parts)]  # direction is last part
    ref_id_part <- str_replace(term_name, paste0("_", direction), "")  # everything except direction
    
    gene_set_groups[[term_name]] <- list(
      name = term_name,
      filter = filtered_data,
      ref_id = ref_id_part,
      direction = direction,
      n_genes = nrow(filtered_data)
    )
    
    print(paste("  ", term_name, ":", nrow(filtered_data), "genes"))
  }
  
  return(gene_set_groups)
}

gene_set_groups <- create_ref_gene_set_groups(ref_df_entrez)

# Print summary
print(paste("\n=== Gene Set Summary ==="))
print(paste("Total gene sets created:", length(gene_set_groups)))
for(gene_set_name in names(gene_set_groups)) {
  gene_set <- gene_set_groups[[gene_set_name]]
  n_genes <- nrow(gene_set$filter)
  print(paste("Gene set:", gene_set_name, "- Genes:", n_genes))
}
################################################################################
# Prepare test genes
################################################################################
df_test = df %>% filter(age_lfsr < 0.05 & sexM_lfsr < 0.05) %>% mutate(direction = case_when(
  quadrant %in% c(1,4) ~ "up",
  quadrant %in% c(2,3) ~ "down")) %>%
  mutate(ct_reg_id = paste(cell_type, region, sep="_")) %>%
  mutate(direction_id = paste(ct_reg_id, quadrant, sep="_"))

df_test = df_test %>% filter(direction == dir) %>% filter(cell_type %in% ref_df$ref_ct_match)

# Filter for significant directional enrichment 
direction_summary = read_csv("/scratch/klwatki5/R_analysis/bbi/mash_new/purple_genes/directional_binom_results_enrich_filt.csv")

sig_result = direction_summary %>% filter(p_adj < 0.05)

sig_ct_region_summary = sig_result %>% group_by(ct_reg_id) %>% summarize(mean_fraction = mean(male_older_fraction))
sig_ct_region_summary = sig_ct_region_summary %>%
  mutate(direction = case_when(
    mean_fraction < 0.5 ~ "female older",
    TRUE ~ "male older"))

fo_ids = sig_ct_region_summary %>% filter(direction == "female older") %>% pull(ct_reg_id)
mo_ids = sig_ct_region_summary %>% filter(direction == "male older") %>% pull(ct_reg_id)

mo_test = df_test %>% filter(quadrant %in% c(1,3) & ct_reg_id %in% mo_ids)
fo_test = df_test %>% filter(quadrant %in% c(2,4) & ct_reg_id %in% fo_ids)

sig_ids = rbind(mo_test, fo_test) %>% pull(direction_id)

df_test = df_test %>% filter(direction_id %in% sig_ids)

ct_list = unique(df_test$cell_type) 

df_test = df_test %>% mutate(
  comparison_group = case_when(
    quadrant %in% c(1) ~ "1_m_age_up_asDEGs",
    quadrant %in% c(3) ~ "3_f_age_down_asDEGs",
    quadrant %in% c(2) ~ "2_m_age_down_asDEGs",
    quadrant %in% c(4) ~ "4_f_age_up_asDEGs"
  )
)
################################################################################
# Run test
################################################################################

run_fisher_test <- function(test_genes, background_genes, ad_genes, term_name) {
  # Create contingency table
  in_test_in_ad <- length(intersect(test_genes, ad_genes))
  in_test_not_ad <- length(setdiff(test_genes, ad_genes))
  not_test_in_ad <- length(setdiff(intersect(background_genes, ad_genes), test_genes))
  not_test_not_ad <- length(setdiff(setdiff(background_genes, ad_genes), test_genes))
  
  # Reorganize the matrix to match the alternative hypothesis
  cont_table <- matrix(
    c(in_test_in_ad, not_test_in_ad,      # AD genes: in test, not in test  
      in_test_not_ad, not_test_not_ad),   # Non-AD genes: in test, not in test
    nrow = 2,
    dimnames = list(
      c("In_Test", "Not_Test"),           # Rows: test status
      c("In_AD", "Not_AD")               # Columns: AD status
    )
  )
  
  print("Contingency table:")
  print(cont_table)
  
  # Run Fisher's exact test
  fisher_result <- fisher.test(cont_table, alternative = "two.sided")
  
  # Extract confidence intervals
  lowerci <- fisher_result$conf.int[1]
  upperci <- fisher_result$conf.int[2]
  
  # Calculate fold enrichment
  expected <- (length(test_genes) * length(intersect(background_genes, ad_genes))) / length(background_genes)
  fold_enrichment <- in_test_in_ad / expected
  
  result <- data.frame(
    Term = term_name,
    GeneRatio = paste0(in_test_in_ad, "/", length(test_genes)),
    BgRatio = paste0(length(intersect(background_genes, ad_genes)), "/", length(background_genes)),
    pvalue = fisher_result$p.value,
    odds_ratio = fisher_result$estimate,
    odds_ratio_lowerci = lowerci,
    odds_ratio_upperci = upperci,
    fold_enrichment = fold_enrichment,
    overlapping_genes = in_test_in_ad,
    in_test_not_ad = in_test_not_ad,
    not_test_in_ad = not_test_in_ad,
    not_test_not_ad = not_test_not_ad
  )
  
  return(result)
}

gene_count_summary <- data.frame()

for(ct in ct_list[1:5]){
  
  matching_gene_sets <- gene_set_groups[sapply(gene_set_groups, function(x) {
    ref_ct <- strsplit(x$ref_id, "_")[[1]][1]  # Extract reference cell type from ref_id
    return(ref_ct == ct)  # Return TRUE if it matches current cell type
  })]
  # Create a container for all results from all gene set groups
  all_results = NULL
  
  for(gene_set in matching_gene_sets) {
    
    # Get the gene set name and AD genes
    dataset_name <- gene_set$name
    ad_genes_for_dataset <- unique(gene_set$filter$ncbi_gene)
    
    # Enhanced: Extract gene set metadata for better tracking
    gene_set_direction <- gene_set$direction
    gene_set_ref_id <- gene_set$ref_id
    gene_set_n_genes <- gene_set$n_genes
    
    data_all = NULL
    ct_subset = df_test %>% filter(cell_type == ct)
    
    for(rg in unique(ct_subset$region)){
      # Define background universe as all genes tested in a Louvain-region
      all_genes_bg = df %>%
        filter(cell_type == ct, region == rg) %>%
        pull(ensembl_gene_id) %>%
        unique()
      all_genes_bg <- gsub("\\..*", "", all_genes_bg)
      
      safe_select <- function(db, keys, cols, kt){
        tryCatch({
          AnnotationDbi::select(db, keys=keys, columns=cols, keytype=kt)
        }, error=function(e){
          warning(paste("Failed select for keys:", paste(head(keys), collapse=", ")))
          return(data.frame())
        })
      }
      
      background.df <- safe_select(org.Mmu.eg.db, all_genes_bg, c("ENTREZID","ENSEMBL"), "ENSEMBL")
      background_genes_entrez <- unique(background.df$ENTREZID)
      
      for(comp_group in c("1_m_age_up_asDEGs","4_f_age_up_asDEGs","3_f_age_down_asDEGs","2_m_age_down_asDEGs")) {
        id = paste(ct, rg, comp_group, sep='_')
        
        
        # Filter data for the current region and comparison group
        temp = ct_subset %>% 
          filter(region == rg, 
                 comparison_group == comp_group)
        
        # Log gene counts for reporting
        gene_count = length(unique(temp$ensembl_gene_id))
        gene_count_summary = rbind(gene_count_summary, 
                                   data.frame(cell_type = ct,
                                              region = rg, 
                                              comparison_group = comp_group, 
                                              gene_count = gene_count,
                                              gene_set_name = dataset_name,
                                              gene_set_direction = gene_set_direction))
        
        if(nrow(temp) > 0) {
          gene_group = unique(temp$ensembl_gene_id)
          gene_count = length(gene_group)
          print(paste("Number of significant a+sDEGs:", gene_count))
          
          if(gene_count >= 3) {
            tryCatch({
              # Get EntrezIDs for the genes
              gene.df <- safe_select(org.Mmu.eg.db, gene_group, c("ENTREZID","ENSEMBL"), "ENSEMBL")
              gene.df <- gene.df %>% 
                group_by(ENSEMBL) %>%
                dplyr::slice(1) %>%
                ungroup() %>%
                distinct(ENTREZID, .keep_all = TRUE)
              test_genes_entrez <- unique(gene.df$ENTREZID)
              
              # Run Fisher's test
              result <- run_fisher_test(
                test_genes = test_genes_entrez,
                background_genes = background_genes_entrez,
                ad_genes = ad_genes_for_dataset,
                term_name = dataset_name
              )
              
              # Enhanced: Add more detailed metadata
              result$cell_type = ct 
              result$region = rg 
              result$comparison_group = comp_group
              result$test_direction = dir
              result$gene_set_n_genes = gene_set_n_genes
              result = result %>%
                mutate(quadrant = case_when(
                  comparison_group == "1_m_age_up_asDEGs" ~ 1,
                  comparison_group == "3_f_age_down_asDEGs" ~ 3,
                  comparison_group == "4_f_age_up_asDEGs" ~ 4,
                  comparison_group == "2_m_age_down_asDEGs" ~ 2
                ))
              
              # Create a unique identifier for debugging
              result$test_id = paste(ct, rg, comp_group, sep="_")
              result$n_test_genes = length(test_genes_entrez)
              
              print(paste("Result for", result$test_id, "vs", dataset_name))
              print(paste("  - overlapping genes:", result$overlapping_genes))
              print(paste("  - p-value:", format(result$pvalue, scientific = TRUE, digits = 3)))
              print(paste("  - fold enrichment:", round(result$fold_enrichment, 2)))
              
              # Add to results
              if(is.null(data_all)) {
                data_all = result
              } else {
                data_all = rbind(data_all, result)
              }
              
            }, error = function(e) {
              print(paste("Error in comparison group", comp_group, ":", e$message))
            })
          } else {
            print(paste("Skipping comparison group", comp_group, "- only", gene_count, "genes (need ≥3)"))
          }
        } else {
          print(paste("No genes found for comparison group", comp_group))
        }
      } # End comparison group loop
      
    } # End region loop
    
    
    # Add data_all to all_results after processing all louvain-region-comparison_group combinations for this gene set
    all_results = rbind(all_results, data_all)
    
  } # End gene_set loop
  
  # Write all results to a single CSV file
  if(!is.null(all_results) && nrow(all_results) > 0) {
    all_results <- all_results %>% mutate(p_adj = p.adjust(pvalue, method = "BH"))
    
    # Use the direction from the selected reference genes for filename
    direction_label <- unique(ref_df_selected$direction)[1]
    combined_output_filename = paste0('/scratch/klwatki5/R_analysis/bbi/mash_new/ad_analysis/ora_results/purple_genes/ct_match_mathys/', direction_label, '/', ct, '_ad_ct_match_', direction_label,'.csv')
    write_csv(all_results, combined_output_filename)
    print(paste("All results saved to:", combined_output_filename))
    
    # Enhanced: Summary statistics
    n_significant = sum(all_results$p_adj < 0.05, na.rm = TRUE)
    n_total_tests = nrow(all_results)
    print(paste("Summary for", ct, ":", n_significant, "significant results out of", n_total_tests, "tests"))
    
  } else {
    print(paste("No enrichment results found for cell type:", ct))
  }
  
}
