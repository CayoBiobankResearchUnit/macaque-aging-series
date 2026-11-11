library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(msigdbrdata)

# Read the data
df = read_csv('/scratch/klwatki5/R_analysis/bbi/mash_new/results/master_merged_update.csv')
df$region = factor(df$region, levels=c('ACC','dlPFC','IPP','M1','EC','HIP','CN','NAc','mdTN','MB','lCb'))

# Convert all genes to human orthologs, get NCBI IDs for human genes 
biomart = read.csv('/scratch/klwatki5/ref_data/var_rna_orthologs_human.txt',sep='\t') %>% 
  filter(hsapiens_homolog_ensembl_gene!='' & hsapiens_homolog_orthology_type == "ortholog_one2one")

df = left_join(df,biomart,by=c('ensembl_gene_id'))
df = df %>% filter(ensembl_gene_id != '') %>% filter(hsapiens_homolog_ensembl_gene != "")

gene.df <- AnnotationDbi::select(org.Hs.eg.db,keys = unique(df$hsapiens_homolog_ensembl_gene),columns = c("ENTREZID", "ENSEMBL"),keytype = "ENSEMBL")
df<-left_join(df,gene.df,by=c('hsapiens_homolog_ensembl_gene'='ENSEMBL'))%>%filter(ENTREZID!='')

ct_list = unique(df$cell_type)

# Fix column names to match code 
df$celltype = df$cell_type
df$gene = df$ensembl_gene_id

# Filter for only autosomes
df = df %>% filter(!chr %in% c("X","Y"))

# Filter for significant a+sDEGs first
df_all = df
df = df %>% 
  filter(age_lfsr < 0.05 &   # Significant age effect
           sexM_lfsr < 0.05)  # Significant sex effect

# Create a new variable for the comparison groups
df = df %>% mutate(
  comparison_group = case_when(
    quadrant %in% c(1) ~ "1_m_age_up_asDEGs",
    quadrant %in% c(3) ~ "3_f_age_down_asDEGs",
    quadrant %in% c(2) ~ "2_m_age_down_asDEGs",
    quadrant %in% c(4) ~ "4_f_age_up_asDEGs"
  )
)

# Keep track of gene counts for reporting
gene_count_summary = data.frame()

################################################################################
#Prepare synGO
################################################################################
syngo_annot = read.csv('/data/CEM/smacklab/data/bican/GO_dataset/syngo_annotations.csv')
syngo_genes = read.csv('/data/CEM/smacklab/data/bican/GO_dataset/syngo_genes.csv')
syngo_ontologies <- read_csv("/data/CEM/smacklab/data/bican/GO_dataset/syngo_ontologies.csv")
syngo_ontologies$setSize <- str_count(syngo_ontologies$hgnc_symbol, ",") + 1

##filter geneset with small set size
select_term = syngo_ontologies%>%filter(setSize>10) %>%pull(id) #use it to filter out terms with small setsize
syngo_annot = syngo_annot%>%filter(go_id %in% select_term)

#convert to clusterprofiler compatible format (TERM2GENE = syngo_t2g)
syngo_df = left_join(syngo_annot,syngo_genes,by=c('hgnc_id','hgnc_symbol'))

biomart = read.csv('/scratch/klwatki5/ref_data/var_rna_orthologs_human.txt',sep='\t')
biomart = biomart%>%filter(hsapiens_homolog_ensembl_gene!='' & hsapiens_homolog_orthology_type == "ortholog_one2one")
syngo_df = left_join(syngo_df,biomart,by=c('ensembl_id'='hsapiens_homolog_ensembl_gene'))
colnames(syngo_df)[15:20] = c('human_ensembl_gene','human_entrez_gene','ensembl_gene','human_gene_symbol','human_orthology_type','gene_symbol')
syngo_df = syngo_df%>%filter(gene_symbol!='')
gene.df <- AnnotationDbi::select(org.Hs.eg.db,keys = unique(syngo_df$human_ensembl_gene),columns = c("ENTREZID", "ENSEMBL"),keytype = "ENSEMBL")
syngo_df<-left_join(syngo_df,gene.df,by=c('human_ensembl_gene'='ENSEMBL'))%>%filter(ENTREZID!='')
syngo_t2g = syngo_df %>% dplyr::distinct(go_name, ENTREZID) %>% as.data.frame()
colnames(syngo_t2g)=c('gs_name','entrez_gene')

# Transform syngo_df to match gene_set_groups structure
format_syngo_df <- function(syngo_df) {
  # Create a new dataframe with the necessary columns
  formatted_df <- data.frame(
    gene_symbol = syngo_df$gene_symbol,
    ncbi_gene = as.integer(syngo_df$ENTREZID),  # Use ENTREZID instead of ncbi_gene
    ensembl_gene = syngo_df$human_ensembl_gene,
    db_gene_symbol = syngo_df$human_gene_symbol,
    db_ncbi_gene = as.character(syngo_df$human_entrez_gene),
    db_ensembl_gene = syngo_df$human_ensembl_gene,
    source_gene = syngo_df$gene_symbol,
    stringsAsFactors = FALSE
  )
  
  # Add required constant columns for the gene set info
  formatted_df$gs_id <- "SynGO"
  formatted_df$gs_name <- paste0("SYNGO_", syngo_df$go_name)
  formatted_df$gs_collection <- "SynGO"
  formatted_df$gs_subcollection <- ""
  formatted_df$gs_collection_name <- "SynGO"
  formatted_df$gs_description <- syngo_df$go_name
  formatted_df$gs_source_species <- "HS"
  formatted_df$gs_pmid <- as.character(syngo_df$pubmed_id)
  formatted_df$gs_geoid <- ""
  formatted_df$gs_exact_source <- ""
  formatted_df$gs_url <- ""
  formatted_df$db_version <- "2024.1.Hs"
  formatted_df$db_target_species <- "HS"
  formatted_df$ortholog_taxon_id <- 9544  # Same as in gene_set_groups --> not accurate, should be human 
  formatted_df$ortholog_sources <- "SynGO"
  formatted_df$num_ortholog_sources <- 1
  
  return(formatted_df)
}

# Use this function to create your formatted dataframe
syngo_formatted <- format_syngo_df(syngo_df)

# Now you can create a list element similar to gene_set_groups
syngo_list <- list(
  name = "syngo",
  filter = syngo_formatted,
  category = "SynGO"
)

################################################################################
# Prepare gene sets
################################################################################
gene_list = msigdbrdata(target_species="HS")
gene_list = left_join(gene_list,biomart,by=c('db_ensembl_gene'='hsapiens_homolog_ensembl_gene'))
gene_list = gene_list %>% filter(ensembl_gene_id != '')

# Create a list of gene set groups with their categories
gene_set_groups = list(
  list(name = "hallmark", filter = gene_list %>% filter(gs_collection == "H"), category = "Hallmark"), 
  list(name = "kegg_medicus", filter = gene_list %>% filter(gs_subcollection == "CP:KEGG_MEDICUS"), category = "KEGG_Medicus"),
  list(name = "biocarta", filter = gene_list %>% filter(gs_subcollection == "CP:BIOCARTA"), category = "BioCarta"),
  list(name = "go_bp", filter = gene_list %>% filter(gs_subcollection == "GO:BP"), category = "GO:BP"),
  list(name = "go_mf", filter = gene_list %>% filter(gs_subcollection == "GO:MF"), category = "GO:MF"),
  list(name = "go_cc", filter = gene_list %>% filter(gs_subcollection == "GO:CC"), category = "GO:CC"),
  list(name = "pid", filter = gene_list %>% filter(gs_subcollection == "CP:PID"), category = "PID"),
  list(name = "reactome", filter = gene_list %>% filter(gs_subcollection == "CP:REACTOME"), category = "Reactome"),
  list(name = "wiki", filter = gene_list %>% filter(gs_subcollection == "CP:WIKIPATHWAYS"), category = "WikiPathways")
)

gene_set_groups[[length(gene_set_groups) + 1]] <- syngo_list

################################################################################
# Run the analysis 
################################################################################
for(ct in ct_list[1:12]){
  # Create a container for all results from all gene set groups
  all_results = NULL
  
  # Loop through each gene set group
  for(gene_set in gene_set_groups) {
    print(paste("Processing gene set group:", gene_set$name))
    
    # Create term to gene mapping for current gene set
    msigdbr_t2g = gene_set$filter %>% dplyr::distinct(gs_name, db_ncbi_gene) %>% as.data.frame()
    
    data_all = NULL
    subset = df %>% filter(celltype == ct)
    
    # Skip if no genes for this cell type
    if(nrow(subset) == 0) {
      print(paste("No significant a+sDEGs found for cell type:", ct))
      next
    }
    
    for(lv in unique(subset$louvain)){
      lv_subset = subset %>% filter(louvain == lv)
      
      # Skip if no genes for this louvain cluster
      if(nrow(lv_subset) == 0) {
        print(paste("No significant a+sDEGs found for louvain cluster:", lv, "in cell type:", ct))
        next
      }
      
      for(rg in unique(lv_subset$region)){
        # Get all background genes for this cell type, louvain and region
        # This should be from the original dataset before significance filtering
        all_genes_bg = df_all %>%
          filter(cell_type == ct, louvain == lv, region == rg) %>%
          pull(hsapiens_homolog_ensembl_gene) %>%
          unique()
        
        # Process each comparison group separately
        for(comp_group in c("1_m_age_up_asDEGs","4_f_age_up_asDEGs","3_f_age_down_asDEGs","2_m_age_down_asDEGs")) {
          id = paste(ct, lv, rg, comp_group, sep='_')
          print(id)
          
          # Filter data for the current region and comparison group
          temp = lv_subset %>% 
            filter(region == rg, 
                   comparison_group == comp_group)
          
          # Log gene counts for reporting
          gene_count = length(unique(temp$hsapiens_homolog_ensembl_gene))
          gene_count_summary = rbind(gene_count_summary, 
                                     data.frame(cell_type = ct, 
                                                louvain = lv, 
                                                region = rg, 
                                                comparison_group = comp_group, 
                                                gene_count = gene_count))
          if(nrow(temp) > 0) {
            gene_group = unique(temp$hsapiens_homolog_ensembl_gene)
            gene_count = length(gene_group)
            print(paste("Number of significant a+sDEGs:", gene_count))
            
            if(gene_count >= 3) {  # Only proceed if there are at least 3 genes
              # Get EntrezIDs for the genes
              tryCatch({
                
                # Join gene IDs with the data
                df_gene = temp
                df_gene = na.omit(df_gene)
                
                
                # Get background gene EntrezIDs
                background.df <- AnnotationDbi::select(org.Hs.eg.db, 
                                                       keys = all_genes_bg,
                                                       columns = c("ENTREZID", "ENSEMBL"),
                                                       keytype = "ENSEMBL") %>%
                  distinct(ENTREZID, .keep_all = TRUE)
                
                # Run over-representation analysis
                ora_result = enricher(
                  gene = unique(df_gene$ENTREZID),
                  universe = unique(background.df$ENTREZID),
                  TERM2GENE = msigdbr_t2g,
                  pvalueCutoff = 1,
                  pAdjustMethod = "BH",
                  minGSSize = 10,
                  maxGSSize = 500
                )
                
                if(!is.null(ora_result)) {
                  df_ora = data.frame(ora_result)
                  print(paste("ORA results dimension:", nrow(df_ora), "x", ncol(df_ora)))
                  
                  if(nrow(df_ora) > 0) {
                    # Add metadata
                    df_ora$cell_type = ct
                    df_ora$louvain = lv
                    df_ora$region = rg
                    df_ora$comparison_group = comp_group
                    df_ora$Category = gene_set$category
                    df_ora$gene_set_group = gene_set$name
                    df_ora$num_input_genes = gene_count
                    
                    # Convert GeneRatio to numeric for better comparison
                    df_ora <- df_ora %>%
                      mutate(GeneRatio_numeric = sapply(strsplit(GeneRatio, "/"), function(x) {
                        as.numeric(x[1]) / as.numeric(x[2])
                      }))
                    
                    # Add effect direction for visualization
                    df_ora <- df_ora %>%
                      mutate(quadrant = case_when(
                        comparison_group == "1_m_age_up_asDEGs" ~ 1,
                        comparison_group == "3_f_age_down_asDEGs" ~ 3,
                        comparison_group == "4_f_age_up_asDEGs" ~ 4,
                        comparison_group == "2_m_age_down_asDEGs" ~ 2
                      ))
                    
                    # Append to results
                    data_all = rbind(data_all, df_ora)
                  }
                }
              }, error = function(e) {
                message(sprintf("Error processing %s: %s", id, e$message))
              })
            } else {
              print(paste("Skipping", id, "- insufficient genes (", gene_count, ")"))
            }
          }
        }
      }
    }
    
    # Append this gene set's results to the combined results
    all_results = rbind(all_results, data_all)
  }
  
  #Write all results to a single CSV file
  if(!is.null(all_results) && nrow(all_results) > 0) {
    combined_output_filename = paste0('/scratch/klwatki5/R_analysis/bbi/mash_new/purple_genes/ora_results/louvain_region/db_update/', ct, '.csv')
    write_csv(all_results, combined_output_filename)
    print(paste("All results saved to:", combined_output_filename))
  } else {
    print(paste("No enrichment results found for cell type:", ct))
  }
}

################################################################################
# Combine cell class results into one master file 
################################################################################

# Set directory path
dir_path <- "/scratch/klwatki5/R_analysis/bbi/mash_new/purple_genes/ora_results/louvain_region/db_update/"

# Get all CSV files in the directory
csv_files <- list.files(dir_path, pattern = "\\.csv$", full.names = TRUE)

ora_results <- csv_files %>%
  set_names(tools::file_path_sans_ext(basename(.)) %>% str_remove("_ct_reg$")) %>%
  map_dfr(read_csv, .id = "cell_type")

write_csv(ora_results, "/scratch/klwatki5/R_analysis/bbi/mash_new/purple_genes/ora_results/louvain_region/db_update/ora_result_merged.csv" )

