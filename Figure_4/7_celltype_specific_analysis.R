source('/home/wyang114/utils/plotting_theme.R')
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

#version after removing HIP
df = read.csv('/scratch/nsnyderm/u01/intermediate_files/nebula/age/all_mashr_summary_lfsr<0.05_filter_missing.csv')[-1]
df$region = factor(df$region,levels=c('ACC','dlPFC','IPP','M1','EC','HIP','CN','NAc','mdTN','MB','lCb'))
df_nebula = read.csv('/scratch/nsnyderm/u01/intermediate_files/nebula/age/nebula_NBGMM_filtered_rminfants.csv')
df_all = read.csv('/scratch/nsnyderm/u01/intermediate_files/nebula/age/all_mashr_summary_filter_missing.csv')[-1]

reg_total = df%>%group_by(cell_type,region)%>%summarise(n=n())%>%group_by(cell_type)%>%summarise(n=n())
#region-speicific GO
count_list = list()
data_all = NULL
perc = 0.6
#wba (genes change in >60% of the regions)
ct_list = unique(df$cell_type)[c(1,4,5,6,7,9,10,11)]
for(ct in ct_list){
  n_total = reg_total%>%filter(cell_type==ct)%>%pull(n)
  for(dir in c('up','down')){
      id = paste0(ct,'_',dir)
      print(id)
      subset = df%>%filter(cell_type==ct)
      temp = subset%>%filter(direction==dir)%>%group_by(region,gene)%>%summarise(mean_pm = mean(age_mashr_pm))
      perc_region = temp%>%group_by(gene)%>%summarise(perc_region=n()/n_total)
      
      gene_group = perc_region%>%filter(perc_region>=perc)%>%pull(gene)%>%unique
      count_list[id]=length(gene_group)
      
      if(length(gene_group)!=0){
        universe = df_nebula%>%filter(cell_type==ct)%>%pull(gene)%>%unique()
        gene.df <- AnnotationDbi::select(org.Mmu.eg.db,keys = gene_group,columns = c("ENTREZID", "ENSEMBL"),keytype = "ENSEMBL")
        df_gene = left_join(temp,gene.df,by=c('gene'='ENSEMBL'))%>%drop_na()
        all.gene.df = AnnotationDbi::select(org.Mmu.eg.db,keys = universe,columns = c("ENTREZID", "ENSEMBL"),keytype = "ENSEMBL")
        
        #change it to enrichKEGG and enrichGO
        compKEGG = enrichKEGG(df_gene$ENTREZID, organism='mcc',universe = all.gene.df$ENTREZID,pvalueCutoff=0.05)
        df_KEGG = data.frame(compKEGG)
        if(dim(df_KEGG)[1]!=0){
          df_KEGG$ONTOLOGY = 'KEGG'   
          df_KEGG$Description = gsub(" - Macaca mulatta \\(rhesus monkey\\)", "", df_KEGG$Description)
        }
        compGO = enrichGO(df_gene$ENTREZID, OrgDb='org.Mmu.eg.db', ont='ALL',universe = all.gene.df$ENTREZID,pvalueCutoff=0.05)
        df_GO = data.frame(compGO)
        df_combine = rbind(df_KEGG,df_GO)
        if(dim(df_combine)[1]!=0){
          df_combine <- df_combine %>%
            mutate(GeneRatio = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))%>%
            mutate(Description=paste0(ONTOLOGY,'-',Description))
          df_combine$cell_type = ct
          df_combine$direction = dir
        }

        data_all = rbind(data_all, df_combine)
      }
  }
}

data_all$ct_dir = paste0(data_all$cell_type,'_',data_all$direction)
data_all%>%group_by(ct_dir)%>%summarise(n_term=n())%>%ggplot(.,aes(ct_dir,n_term))+geom_bar(stat='identity')
write.csv(data_all, paste0('~/output/transcriptional/celltype_region_patterns/all_wba_go_',perc,'.csv'))

perc=0.6
data_all = read.csv(paste0('~/output/transcriptional/celltype_region_patterns/all_wba_go_',perc,'.csv'))[-1]

dir = 'up'
go_filter = data_all%>%filter(direction==dir)%>%filter(p.adjust<0.1)
go_filter%>%group_by(Description,cell_type)%>%summarise(n=n())%>%
  group_by(Description)%>%summarise(n_celltype=n())%>%pull(n_celltype)%>%hist(main='number of cell type per term')
temp =go_filter%>%group_by(Description,cell_type)%>%summarise(n=n())%>%
  group_by(Description)%>%summarise(n_celltype=n())
share_term = temp%>%filter(n_celltype==2)%>%pull(Description) #12
go_share = data_all%>%filter(Description %in% share_term)


top_terms <- go_filter %>%
  filter(direction==dir)%>%
  filter(p.adjust<0.1)%>%
  arrange(p.adjust, desc(GeneRatio)) %>%
  group_by(cell_type) %>%
  slice_min(order_by = row_number(), n = 5)%>%
  rbind(go_share)

# Create dotplot
all_ct = unique(df$celltype)[c(1,10,7,11,4,6,5,9)]
ggplot(top_terms, aes(x = factor(cell_type,levels =all_ct), y = Description)) +
    geom_point(aes(size = GeneRatio, color = p.adjust)) +
    scale_color_gradient(low = "red", high = "blue", name = "Adjusted P-value") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8)
    ) +
    labs(
      x = "Cluster",
      y = "GO Term",
      size = "Gene Ratio",
      title = paste0(dir,'_region_GO')
    )

#region-specific patterns
all_ct = unique(df$cell_type)[c(1,4,5,6,7,9,10,11)]
regions = c('ACC','dlPFC','IPP','M1','EC','HIP','CN','NAc','mdTN','MB','lCb')
df$region = factor(df$region, levels=regions)
all_go = NULL
k_bests = c(3,3,3,4,5,5,7,5)
for(i in 1:length(all_ct)){
  
  mode = 'continuous_pm'
  ct = all_ct[i]
  
  if(mode =='continuous_pm'){
    #continuous_pm
    mat = df%>%filter(celltype==ct)%>%group_by(region,gene)%>%summarise(mean_pm = mean(age_mashr_pm))%>%
      pivot_wider(names_from=region,values_from=mean_pm)%>%column_to_rownames('gene')
    mat[is.na(mat)]=0
  }else{
    #binary_pm
    mat = df%>%filter(celltype==ct)%>%group_by(region,gene)%>%summarise(mean_pm = mean(age_mashr_pm))%>%
      mutate(binary_mean_pm=case_when(mean_pm > 0 ~ 1,
                                      mean_pm < 0 ~ -1))%>%
      ungroup%>%
      dplyr::select(region,gene,binary_mean_pm)%>%
      pivot_wider(names_from=region,values_from=binary_mean_pm)%>%column_to_rownames('gene')
    mat[is.na(mat)]=0
  }

  
  #calculate the optimal number of clusters with 30 different metrics
  # library(NbClust)
  # start = Sys.time()
  # nb = NbClust(data = mat, diss = NULL, distance = "euclidean",
  #         min.nc = 3, max.nc = 15, method = "ward.D")
  # end = Sys.time()
  # print(end-start)
  # k_frequencies <- table(nb$Best.nc[1, ])
  # pdf(paste0('~/output/transcriptional/celltype_region_patterns/',mode,'/',ct,'_k_frequencies.pdf'), width = 7, height = 5)
  # barplot(k_frequencies,
  #         main = "Frequency of Optimal k across Indices",
  #         xlab = "Number of Clusters (k)",
  #         ylab = "Frequency",
  #         col = "skyblue",
  #         border = "black")
  # dev.off()
  # k_best = names(k_frequencies)[which.max(k_frequencies)]%>%as.numeric()
  ##Using pre-existing k-bests
  k_best = k_bests[i]
  
  #heatmap
  color_palette <- colorRampPalette(c("blue", "white", "red"))(50)
  breaks <- seq(-0.03, 0.03, length.out = 51)
  ph <- pheatmap(mat,clustering_method = "ward.D",show_rownames = FALSE,cluster_cols = FALSE)
  hc <- ph$tree_row
  row_clusters <- cutree(hc, k = k_best)
  row_group_colors <- data.frame(Group = factor(row_clusters))
  row_group_colors$Group%>%table
  rownames(row_group_colors) <- rownames(mat)

  # pheatmap(mat,cluster_rows = TRUE,cluster_cols = FALSE,show_rownames = FALSE,
  #          clustering_method = 'ward.D',
  #          #breaks=breaks,
  #          color = color_palette,
  #          annotation_row = row_group_colors,
  #          treeheight_row = 20,
  #          treeheight_col = 20,main = ct,
  #          filename=paste0('~/output/transcriptional/celltype_region_patterns/',mode,'/',ct,'_heatmap_wardD_cl',k_best,'.pdf'))
  
  
  #GO term for each cluster (p.adjust<0.2,top 10 terms per cluster among GO and KEGG)
  row_group_colors = row_group_colors%>%rownames_to_column('gene')
  gene_group = row_group_colors$gene%>%unique
  universe = df_nebula%>%filter(cell_type==ct)%>%pull(gene)%>%unique()
  gene.df <- AnnotationDbi::select(org.Mmu.eg.db,keys = gene_group,columns = c("ENTREZID", "ENSEMBL"),keytype = "ENSEMBL")
  df_gene = left_join(row_group_colors,gene.df,by=c('gene'='ENSEMBL'))%>%drop_na()
  all.gene.df = AnnotationDbi::select(org.Mmu.eg.db,keys = universe,columns = c("ENTREZID", "ENSEMBL"),keytype = "ENSEMBL")
  compKEGG = compareCluster(geneClusters=ENTREZID~Group, data = df_gene, fun = "enrichKEGG",
                            organism='mcc',universe = all.gene.df$ENTREZID,
                            pvalueCutoff=0.2)
  df_KEGG = data.frame(compKEGG)
  df_KEGG$ONTOLOGY = 'KEGG'
  df_KEGG$Description = gsub(" - Macaca mulatta \\(rhesus monkey\\)", "", df_KEGG$Description)
  compGO =compareCluster(geneClusters=ENTREZID~Group, data = df_gene, fun = "enrichGO",OrgDb='org.Mmu.eg.db',
                         ont='ALL',universe = all.gene.df$ENTREZID,pvalueCutoff=0.2)
  df_GO = data.frame(compGO)
  df_combine = rbind(df_KEGG,df_GO)
  df_combine <- df_combine %>%
    mutate(GeneRatio = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))%>%
    mutate(Description=paste0(ONTOLOGY,'-',Description))
  ##filter out GO term with >90 characters
  #df_combine <- df_combine[nchar(df_combine$Description) <= 90, ]
  
  ##(optional)add pvalue filter
  #df_combine <- df_combine[df_combine$p.adjust<0.1, ]
  
  # if(!is.null(df_combine)){
  #   # Select top terms for each cluster (e.g., top 10 terms per cluster)
  #   top_terms <- df_combine %>%
  #     group_by(Group) %>%
  #     slice_min(order_by = p.adjust, n = 10)  # Adjust 'n' to control top terms per cluster
  #   
  #   # Create dotplot
  #   ggplot(top_terms, aes(x = Group, y = Description)) +
  #     geom_point(aes(size = GeneRatio, color = p.adjust)) +
  #     scale_color_gradient(low = "red", high = "blue", name = "Adjusted P-value") +
  #     theme_minimal() +
  #     theme(
  #       axis.text.x = element_text(angle = 45, hjust = 1),
  #       axis.text.y = element_text(size = 8)
  #     ) +
  #     labs(
  #       x = "Cluster",
  #       y = "GO Term",
  #       size = "Gene Ratio",
  #       title = paste0(ct,'_region_GO')
  #     )
  #  ggsave(paste0('/home/wyang114/output/transcriptional/celltype_region_patterns/',mode,'/',ct,'_region_GO.pdf'),width = 7,height = 7)
  # }
  
  df_combine$cell_type = ct
  all_go = rbind(all_go, df_combine)

  ##number of regions
  # n_region=ncol(mat)
  # mat <- mat %>%rownames_to_column('gene')%>%
  #   rowwise() %>%
  #   mutate(count_gt_zero = sum(!is.na(c_across(2:(n_region+1))))) %>%
  #   filter(count_gt_zero>0)%>%
  #   ungroup()
  # mat%>%ggplot(.,aes(x=count_gt_zero))+geom_histogram(binwidth = 0.5)+
  #   stat_bin(binwidth = 0.5, geom = "text", aes(label = ifelse(..count.. > 0, paste0(round(..count.. / sum(..count..) * 100, 1), "%"), "")), 
  #            vjust = -0.5,size=2)+
  #   scale_x_continuous(breaks = seq(1, 11, 1))+xlab('num_regions')+ylab('n_DEG')+
  #   scale_y_continuous(expand = expansion(mult = c(0, 0.11)))+ggtitle(ct)
  # ggsave(paste0('~/output/transcriptional/celltype_region_patterns/',ct,'_nregion_histogram.pdf'),width = 3,height = 2)
}

write.csv(all_go,'~/output/transcriptional/celltype_region_patterns/all_go_summary.csv')

all_go = read.csv('~/output/transcriptional/celltype_region_patterns/all_go_summary.csv')[-1]
all_go$celltype_Group = paste0(all_go$cell_type,'_',all_go$Group)
go_filter = all_go%>%filter(p.adjust<0.1)
go_filter%>%group_by(Description,cell_type)%>%summarise(n=n())%>%
  group_by(Description)%>%summarise(n_celltype=n())%>%pull(n_celltype)%>%hist(main='number of cell type per term')
temp =go_filter%>%group_by(Description,cell_type)%>%summarise(n=n())%>%
  group_by(Description)%>%summarise(n_celltype=n())
share_term = temp%>%filter(n_celltype>6)%>%pull(Description) #12
go_share = all_go%>%filter(Description %in% share_term)


###plot
top_terms <- go_filter %>%
  group_by(celltype_Group) %>%
  slice_min(order_by = p.adjust, n = 1)%>%  
  rbind(go_share)

# Create dotplot
ggplot(top_terms, aes(x = celltype_Group, y = Description)) +
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  scale_color_gradient(low = "red", high = "blue", name = "Adjusted P-value") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8)
  ) +
  labs(
    x = "Cluster",
    y = "GO Term",
    size = "Gene Ratio",
    title = paste0('all_region_GO')
  )
