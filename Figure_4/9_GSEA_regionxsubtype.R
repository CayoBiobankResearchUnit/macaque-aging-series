source('/home/wyang114/utils/plotting_theme.R')
source('/home/wyang114/utils/functions.R')
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

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Please provide a cell type as a command line argument.")
}
numeric_args <- as.numeric(args)
ct_num <- numeric_args[1]
min_lv = numeric_args[2]
max_lv <- numeric_args[3]


##version after removing HIP
#df = read.csv('/scratch/nsnyderm/u01/intermediate_files/nebula/age/all_mashr_summary_lfsr<0.05_filter_missing.csv')[-1]
df = read.csv('/scratch/nsnyderm/u01/intermediate_files/nebula/age_olig_split//all_mashr_summary_filter_missing_update.csv')[-1]
df$region = factor(df$region,levels=c('ACC','dlPFC','IPP','M1','EC','HIP','CN','NAc','mdTN','MB','lCb'))
ct_list = unique(df$celltype)#[c(1,4,5,6,7,9,10,11)]
# [1] "astrocytes"            "basket cells"          "cerebellar neurons"    "ependymal cells"      
# [5] "GABAergic neurons"     "glutamatergic neurons" "microglia"             "midbrain neurons"     
# [9] "medium spiny neurons"  "oligodendrocytes"      "vascular cells"   
#print(sprintf('processing:%s %s-%s',ct_list[ct_num],min_lv,max_lv))

##Prepare hallmark set
#hallmark_set = msigdbr(species = "Macaca mulatta", category = "H")
#msigdbr_t2g = hallmark_set %>% dplyr::distinct(gs_name, entrez_gene) %>% as.data.frame()

#c2_st = msigdbr(species = "Macaca mulatta", category = "C2")
#msigdbr_t2g = c2_st%>%filter(grepl('SASP',gs_name))%>% dplyr::distinct(gs_name, entrez_gene) %>% as.data.frame()

##Prepare synGO
# syngo_annot = read.csv('/home/wyang114/data/GO_dataset/syngo_annotations.csv')
# syngo_genes = read.csv('/home/wyang114/data/GO_dataset/syngo_genes.csv')
# syngo_ontologies <- read_csv("/home/wyang114/data/GO_dataset/syngo_ontologies.csv")
# syngo_ontologies$setSize <- str_count(syngo_ontologies$hgnc_symbol, ",") + 1

##filter geneset with small set size
#select_term = syngo_ontologies%>%filter(setSize>10)%>%pull(id)
#syngo_annot = syngo_annot%>%filter(go_id %in% select_term)
#visualize if we can group some of these
# syngo_ontologies$sets <- strsplit(syngo_ontologies$hgnc_symbol, ",")
# jc_mat = matrix(NA, nrow = nrow(syngo_ontologies), ncol = nrow(syngo_ontologies))
# for (i in 1:nrow(syngo_ontologies)) {
#   for(j in 1: nrow(syngo_ontologies)){
#     jc_mat[i,j] <- jaccard(syngo_ontologies$sets[[i]], syngo_ontologies$sets[[j]]) 
#   }
# }
# colnames(jc_mat)=syngo_ontologies$name
# row.names(jc_mat)=syngo_ontologies$name
# pheatmap(jc_mat,fontsize = 1)


# syngo_df = left_join(syngo_annot,syngo_genes,by=c('hgnc_id','hgnc_symbol'))
# biomart = read.csv('~/data/mart_export_human_primate.txt',sep='\t')
# biomart = biomart%>%filter(Human.gene.stable.ID!='')
# syngo_df = left_join(syngo_df,biomart,by=c('ensembl_id'='Human.gene.stable.ID'))
# colnames(syngo_df)[15:19] = c('human_ensembl_gene','human_entrez_gene','ensembl_gene','human_gene_symbol','gene_symbol')
# syngo_df = syngo_df%>%filter(gene_symbol!='')
# gene.df <- AnnotationDbi::select(org.Mmu.eg.db,keys = unique(syngo_df$ensembl_gene),columns = c("ENTREZID", "ENSEMBL"),keytype = "ENSEMBL")
# syngo_df<-left_join(syngo_df,gene.df,by=c('ensembl_gene'='ENSEMBL'))%>%filter(ENTREZID!='')
# syngo_t2g = syngo_df %>% dplyr::distinct(go_name, ENTREZID) %>% as.data.frame()
# colnames(syngo_t2g)=c('gs_name','entrez_gene')

for(ct in ct_list[ct_num]){
  print(ct)
  data_all = NULL
  louvains = unique(df%>%filter(celltype==ct)%>%pull(louvain))
  if(min_lv==0){
    min_lv = 1  
  }
  if(max_lv==0){
    max_lv = length(louvains)  
  }
    
  for(lv in louvains[min_lv:max_lv]){
    regions = unique(df%>%filter(celltype==ct & louvain==lv)%>%pull(region))
    
    for(rg in regions){
      id = paste(ct,lv,rg,sep='_')
      print(id)
      temp = df%>%filter(celltype==ct & louvain==lv & region == rg)
      #temp <- temp %>%
      #  group_by(direction) %>%
      #  arrange(age_mashr_lfsr, .by_group = TRUE) %>%
      #  slice_head(n = 500)

      
      gene.df <- AnnotationDbi::select(org.Mmu.eg.db,keys = unique(temp$gene),columns = c("ENTREZID", "ENSEMBL"),keytype = "ENSEMBL")
      df_gene = left_join(temp,gene.df,by=c('gene'='ENSEMBL'))%>%drop_na()
      df_gene = df_gene%>%mutate(sig_score = case_when(
        age_mashr_pm >0 ~ -log10(age_mashr_lfsr+1e-43),
        age_mashr_pm <0 ~ log10(age_mashr_lfsr+1e-43)
      ))
      gene_list = set_names(df_gene$sig_score,df_gene$ENTREZID)%>%sort(.,decreasing=TRUE)
      
      #gseGO:32s for 200 genes, 49s for 1042 genes, 2min for 4744 genes
      start_time <- Sys.time()
      compGO = gseGO(gene_list,ont = 'ALL',pvalueCutoff = 0.05,OrgDb = 'org.Mmu.eg.db',by = 'fgsea') #pAdjustMethod = "BH"
      end_time <- Sys.time()
      print(end_time-start_time)
      df_gsea = data.frame(compGO)
      print(dim(df_gsea))
      
      ##gseKEGG
      #start_time <- Sys.time()
      #compKEGG = gseKEGG(gene_list,organism='mcc',pvalueCutoff = 0.05)
      #end_time <- Sys.time()
      #print(end_time-start_time)
      #df_gsea = data.frame(compKEGG)
      #print(dim(df_gsea))
      
      # # #msigdb-hallmark:15-20s for >4000 genes
      # start_time <- Sys.time()
      # compGSEA = GSEA(gene = gene_list, TERM2GENE = msigdbr_t2g, pvalueCutoff = 0.05)
      # end_time <- Sys.time()
      # print(end_time-start_time)
      # df_gsea = data.frame(compGSEA)
      # print(dim(df_gsea))

      
      # #gsesynGO
      # start_time <- Sys.time()
      # compGSEA = GSEA(gene = gene_list, TERM2GENE = syngo_t2g, pvalueCutoff = 0.05)
      # end_time <- Sys.time()
      # print(end_time-start_time)
      # df_gsea = data.frame(compGSEA)
      # print(dim(df_gsea))

      if(dim(df_gsea)[1]!=0){
        df_gsea$cell_type = ct
        df_gsea$louvain = lv
        df_gsea$region=rg
          }
      data_all = rbind(data_all, df_gsea)
        }
  }
  print('saving files...')
  write.csv(data_all, sprintf('/scratch/wyang114/outputs/transcriptional/GO_terms/subtype_region_gsea_C5/%s_all_go_subtypexregion_gsea_%s_%s.csv',ct,min_lv,max_lv))
}


