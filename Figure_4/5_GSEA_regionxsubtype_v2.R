#' @file 4_mashr_summary_clean.R
#' @title GSEA on age-associated DEGs
#' @author Wei Yang
#'
#' @description
#' Gene set enrichment of GSEA of GO(C5), KEGG, hallmark, synGO, SASP(C2) of age-associated DEGs 
#'
#' @section Inputs:
#' \itemize{
#'   \item \code{all_mashr_summary_lfsr<0.05_filter_missing_update.csv}: Combined dataset containing 
#'   both NEBULA and \code{mashr} statistics.
#' }
#'
#' @section Outputs:
#' \itemize{
#'   \item \code{<cell class>_all_gsea_subtypexregion_gsea_<louvain>_thre<pvalueCutoff>.csv}: 
#'   GSEA statistics.
#' }
#'
#' @usage
#' #SBATCH --time=0-01:00:00      #upper bound time limit for job to finish\
#' #SBATCH --mem=200G
#' module load gsl-2.7.1-gcc-11.2.0
#' module load r-4.2.2-gcc-11.2.0
#' module load gcc-11.2.0-gcc-11.2.0
#' Rscript 11_GSEA_regionxsubtype.R $min_louvain $max_louvain $pvalueCutoff
NULL


args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Please provide a cell type as a command line argument.")
}
numeric_args <- as.numeric(args)
ct_num <- numeric_args[1]
min_lv = numeric_args[2]
max_lv <- numeric_args[3]
thre <- numeric_args[4]

source("utils/plotting_theme.R")
source("utils/functions.R")
library(dplyr)
library(tidyverse)
library(ggplot2)
library(stringr)
library(pheatmap)
library(gridExtra)
library(ggpubr)
library(EnhancedVolcano)
library(ggpubr)
library(clusterProfiler)
organism = "org.Mmu.eg.db"
library(organism, character.only = TRUE)
library(RColorBrewer)
library(msigdbr)

set.seed(2025)

# ct_num <- 6
# min_lv = 1
# max_lv <- 30
# thre <- 0.1

df = read.csv('/scratch/wyang114/outputs/nebula/all_mashr_summary_lfsr<0.05_filter_missing_update.csv')[-1]
df$region = factor(df$region,levels=c('ACC','dlPFC','IPP','M1','EC','HIP','CN','NAc','mdTN','MB','lCb'))
ct_list = sort(unique(df$cell_type))
# [1] "astrocytes"            "basket cells"          "cerebellar neurons"    "ependymal cells"       "GABAergic neurons"    
# [6] "glutamatergic neurons" "medium spiny neurons"  "microglia"             "midbrain neurons"      "oligodendrocytes"     
# [11] "OPC"                   "vascular cells"      
print(sprintf('processing:%s %s-%s',ct_list[ct_num],min_lv,max_lv))

#Prepare hallmark set
hallmark_set = msigdbr(species = "Macaca mulatta", category = "H")
msigdbr_t2g_hall = hallmark_set %>% dplyr::distinct(gs_name, entrez_gene) %>% as.data.frame()

##Prepare SASP set
c2_st = msigdbr(species = "Macaca mulatta", category = "C2")
msigdbr_t2g_sasp = c2_st%>%filter(grepl('FAS',gs_name))%>% dplyr::distinct(gs_name, entrez_gene) %>% as.data.frame()

##Prepare synGO
syngo_annot = read.csv('/home/wyang114/data/GO_dataset/syngo_annotations.csv')
syngo_genes = read.csv('/home/wyang114/data/GO_dataset/syngo_genes.csv')
syngo_ontologies <- read_csv("/home/wyang114/data/GO_dataset/syngo_ontologies.csv")
syngo_ontologies$setSize <- str_count(syngo_ontologies$hgnc_symbol, ",") + 1
select_term = syngo_ontologies%>%filter(setSize>10)%>%pull(id)
syngo_annot = syngo_annot%>%filter(go_id %in% select_term)
syngo_df = left_join(syngo_annot,syngo_genes,by=c('hgnc_id','hgnc_symbol'))
biomart = read.csv('~/data/mart_export_human_primate.txt',sep='\t')
biomart = biomart%>%filter(Human.gene.stable.ID!='')
syngo_df = left_join(syngo_df,biomart,by=c('ensembl_id'='Human.gene.stable.ID'))
colnames(syngo_df)[15:19] = c('human_ensembl_gene','human_entrez_gene','ensembl_gene','human_gene_symbol','gene_symbol')
syngo_df = syngo_df%>%filter(gene_symbol!='')
gene.df <- AnnotationDbi::select(org.Mmu.eg.db,keys = unique(syngo_df$ensembl_gene),columns = c("ENTREZID", "ENSEMBL"),keytype = "ENSEMBL")
syngo_df<-left_join(syngo_df,gene.df,by=c('ensembl_gene'='ENSEMBL'))%>%filter(ENTREZID!='')
syngo_t2g = syngo_df %>% dplyr::distinct(go_name, ENTREZID) %>% as.data.frame()
colnames(syngo_t2g)=c('gs_name','entrez_gene')



gsea_folder = '/scratch/wyang114/outputs/transcriptional/GO_terms/'
dir.create(gsea_folder)

for(ct in ct_list[ct_num]){
  print(ct)
  data_all = NULL
  df_subset <- df %>% filter(cell_type == ct)
  louvains = unique(df_subset$louvain)
  if(min_lv==0){
    min_lv = 1  
  }
  if(max_lv==0){
    max_lv = length(louvains)  
  }
  
  #create folder for each subtype
  out_folder= sprintf(paste0(gsea_folder,ct))
  if (!dir.exists(out_folder)) {
  dir.create(out_folder)
  } 
    
  for(lv in louvains[min_lv:max_lv]){
    data_all = NULL
    regions = unique(df_subset%>%filter(louvain==lv)%>%pull(region))
    
    for(rg in regions){
      id = paste(ct,lv,rg,sep='_')
      print(id)
      temp = df_subset%>%filter(louvain==lv & region == rg)

      
      gene.df <- AnnotationDbi::select(org.Mmu.eg.db,keys = unique(temp$ensembl_gene_id),columns = c("ENTREZID", "ENSEMBL"),keytype = "ENSEMBL")
      df_gene = left_join(temp,gene.df,by=c('ensembl_gene_id'='ENSEMBL'))%>%drop_na()
      df_gene = df_gene%>%mutate(sig_score = case_when(
        age_pm >0 ~ -log10(age_lfsr+1e-43),
        age_pm <0 ~ log10(age_lfsr+1e-43)
      ))
      gene_list = set_names(df_gene$sig_score,df_gene$ENTREZID)%>%sort(.,decreasing=TRUE)
      
      #gseGO:32s for 200 genes, 49s for 1042 genes, 2min for 4744 genes
      start_time <- Sys.time()
      compGO = gseGO(gene_list,ont = 'ALL',pvalueCutoff = thre,OrgDb = 'org.Mmu.eg.db',by = 'fgsea', pvalueCutoff =thre) #pAdjustMethod = "BH"
      end_time <- Sys.time()
      print(end_time-start_time)
      df_go = data.frame(compGO)
      print(dim(df_go))

      #gseKEGG
      start_time <- Sys.time()
      compKEGG = gseKEGG(gene_list,organism='mcc',pvalueCutoff = thre)
      end_time <- Sys.time()
      print(end_time-start_time)
      df_kegg = data.frame(compKEGG)
      df_kegg$ONTOLOGY="KEGG"
      print(dim(df_kegg))
      
      #msigdb-SASP/hallmark:15-20s for >4000 genes
      start_time <- Sys.time()
      compGSEA = GSEA(gene = gene_list, TERM2GENE = msigdbr_t2g_hall, pvalueCutoff = thre)
      end_time <- Sys.time()
      print(end_time-start_time)
      df_sasp = data.frame(compGSEA)
      df_sasp$ONTOLOGY="SASP"
      print(dim(df_sasp))

      
      #gsesynGO
      start_time <- Sys.time()
      compGSEA = GSEA(gene = gene_list, TERM2GENE = syngo_t2g, pvalueCutoff = thre)
      end_time <- Sys.time()
      print(end_time-start_time)
      df_syngo = data.frame(compGSEA)
      df_syngo$ONTOLOGY="synGO"
      print(dim(df_syngo))
      
      df_gsea = bind_rows(df_go,df_kegg,df_sasp,df_syngo)
      
      if(dim(df_gsea)[1]!=0){
        df_gsea$cell_type = ct
        df_gsea$louvain = lv
        df_gsea$region=rg
        data_all = rbind(data_all, df_gsea)
          }
        }
    out_file <- sprintf('%s/%s_all_gsea_subtypexregion_gsea_%s_thre%s.csv',out_folder,ct,lv,thre)
    print('saving files...')
    write.table(data_all, out_file, sep = ",", row.names = FALSE, 
                col.names = !file.exists(out_file), append = TRUE)
  }
  
}
