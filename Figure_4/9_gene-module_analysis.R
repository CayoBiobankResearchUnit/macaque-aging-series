source("utils/plotting_theme.R")
source("utils/annotation.R")
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
library(NbClust)
regions = c('ACC','dlPFC','IPP','M1','EC','HIP','CN','NAc','mdTN','MB','lCb')

df = read.csv('/scratch/nsnyderm/u01/intermediate_files/nebula/age/all_mashr_summary_lfsr<0.05_filter_missing.csv')[-1]
df_all = read.csv('/scratch/nsnyderm/u01/intermediate_files/nebula/age/all_mashr_summary.csv')[-1]

################
#test on binary (does not scale up well)
################
ct = unique(df$cell_type)[1]
subset_mat = df%>%filter(celltype==ct)%>%mutate(id=paste0(region,'_',subtype))%>%
  pivot_wider(id_cols = 'id',names_from='gene',values_from = 'age_mashr_pm')%>%
  column_to_rownames('id')
temp=subset_mat/abs(subset_mat)
temp[is.na(temp)]=0

cor_matrix <- cor(temp,method='spearman') 
pheatmap(cor_matrix,show_rownames = F, show_colnames = F,clustering_method = 'ward.D')

######################
#by subtype x region
######################
condition_cutoff=10

for(ct in unique(df$cell_type)[1]){
  for(dir in c('up','down')){
  subset_mat = df%>%filter(direction==dir & celltype==ct)%>%
    pivot_wider(id_cols = c('region','subtype'),names_from='gene',values_from = 'age_mashr_pm')
  #use all coefficients
  genes = colnames(subset_mat)
  subset_mat = df_all%>%filter(direction==dir & celltype==ct)%>%
    pivot_wider(id_cols = c('region','subtype'),names_from='gene',values_from = 'age_mashr_pm')
  subset_mat = subset_mat[,genes]
  subset_mat[is.na(subset_mat)]=0
  
  png(sprintf('%s_%s',ct,dir), width = 550, height = 400)
  hist(colSums(subset_mat[,3:dim(subset_mat)[2]]!=0),breaks = 100,
       xlab = 'num of subtype x region', ylab='freq of genes',main=sprintf('%s_%s',ct,dir))
  dev.off()
  
  for(cc in condition_cutoff){
    shared_gene = colnames(subset_mat)[colSums(subset_mat!=0)>cc]
    raw_mat = subset_mat[shared_gene[3:length(shared_gene)]]
    cor_matrix <- cor(raw_mat,method='spearman')
    breaks <- seq(-1, 1, length.out = 51)
    color_palette <- colorRampPalette(c("blue", "white", "red"))(50)
    
    #optimal hierarchical clustering
    dist_mat = dist(cor_matrix)
    nb = NbClust(data = t(raw_mat) , distance = NULL, diss = dist_mat,
                 min.nc = 3, max.nc = 6, method = "ward.D")
    k_frequencies <- table(nb$Best.nc[1, ])
    k_set = names(k_frequencies)[which.max(k_frequencies)]%>%as.numeric()
    row_clust <- hclust(dist(cor_matrix),method = 'ward.D')
    row_annotation <- data.frame(Cluster = factor(cutree(row_clust, k = k_set)))
    rownames(row_annotation) <- rownames(cor_matrix)
    
    pheatmap(cor_matrix,show_rownames = F, show_colnames = F,clustering_method = 'ward.D',
             border_color = 'grey80',
             cutree_rows = k_set,cutree_cols = k_set,
             annotation_row = row_annotation,
             main=sprintf('microglia_by_regionxsubtype_%s_%s_clu%s',dir,cc,k_set),
             breaks=breaks,color = color_palette,
             filename = sprintf('/home/wyang114/output/transcriptional/transcriptional_summary/gene-gene_correlation/%s_by_regionxsubtype_%s_%s_clu%s.pdf',ct,dir,cc,k_set),
             width=5.5,height=4
    )
             
    write.csv(row_annotation,
              sprintf('/home/wyang114/output/transcriptional/transcriptional_summary/gene-gene_correlation/%s_by_regionxsubtype_%s_%s_clu%s',ct,dir,cc,k_set))
    
    
    }
  }
}

#HC
library(dynamicTreeCut)
library(viridis)
coef_mat <- df%>%filter(celltype=='astrocytes')%>%
  filter(direction=='up')%>%dplyr::select(gene,region,subtype,age_mashr_pm)%>%
  pivot_wider(id_cols = 'gene', names_from = c('region','subtype'),values_from = 'age_mashr_pm')%>%
  column_to_rownames('gene')
coef_mat[is.na(coef_mat)]=0
binary_matrix=matrix(as.numeric(coef_mat>0), nrow = nrow(coef_mat), ncol = ncol(coef_mat))
row.names(binary_matrix)=row.names(coef_mat)
colnames(binary_matrix)=colnames(coef_mat)
dist_matrix <- dist(binary_matrix)
hc <- hclust(dist_matrix,method = 'ward.D')
dynamic_clusters <- cutreeDynamic(dendro = hc, distM = as.matrix(dist_matrix), 
                                  deepSplit = 4, pamRespectsDendro = FALSE)
annotation <- data.frame(Cluster = factor(dynamic_clusters))
rownames(annotation) <- rownames(binary_matrix)
pheatmap(binary_matrix,
         cluster_rows = hc,               # Use hierarchical clustering for rows
         cluster_cols = FALSE,            # Avoid clustering columns (optional)
         show_rownames = F,show_colnames = F,
         color = viridis(100),
         annotation_row = annotation)     # Add row annotation)
hist(dynamic_clusters)
#matrix of cluster x (region x subtype)
one_hot <- model.matrix(~ Cluster - 1, data = annotation)%>%t()%>%as.matrix()
group_mat = one_hot %*% binary_matrix 
enrich_mat = group_mat/ rowSums(one_hot)
expected_perc = colSums(binary_matrix)/nrow(binary_matrix)
enrich_mat = t(t(enrich_mat)/expected_perc)
col_annotation = data.frame('id'=colnames(enrich_mat))%>%mutate('region'=str_split_fixed(id,'_',3)[,1],
                                                           'celltype'= str_split_fixed(id,'_',3)[,2],
                                                           'subtype'=str_split_fixed(id,'_',2)[,2])%>%
  column_to_rownames('id')
pheatmap(enrich_mat,show_colnames = F,clustering_method = 'ward.D',
         annotation_col = col_annotation)

df_enrich = enrich_mat%>%t()%>%data.frame()%>%rownames_to_column('id')%>%
  pivot_longer(cols=2:17,names_to = 'factor',values_to = 'enrich_score')%>%
  mutate('region'=str_split_fixed(id,'_',3)[,1],
         'celltype'= str_split_fixed(id,'_',3)[,2],
         'subtype'=str_split_fixed(id,'_',2)[,2])

df_enrich$region = factor(df_enrich$region,levels=regions)
df_enrich%>%ggplot(.,aes(subtype,enrich_score,fill=factor))+geom_boxplot()+facet_wrap(~factor,scales = 'free_y')
df_enrich%>%ggplot(.,aes(region,enrich_score,fill=factor))+geom_boxplot()+facet_wrap(~factor,scales = 'free_y')


#NMF
if (!requireNamespace("NMF", quietly = TRUE)) {
  install.packages("NMF")
}
library(NMF)
set.seed(123)  # For reproducibility
coef_mat <- df%>%filter(celltype=='astrocytes')%>%
  filter(direction=='up')%>%dplyr::select(gene,region,subtype,age_mashr_pm)%>%
                    pivot_wider(id_cols = 'gene', names_from = c('region','subtype'),values_from = 'age_mashr_pm')%>%
                    column_to_rownames('gene')
coef_mat[is.na(coef_mat)]=0
pheatmap(coef_mat,show_rownames = F, show_colnames = F,clustering_method = 'ward.D')
binary_matrix=matrix(as.numeric(coef_mat>0), nrow = nrow(coef_mat), ncol = ncol(coef_mat))
row.names(binary_matrix)=row.names(coef_mat)
colnames(binary_matrix)=colnames(coef_mat)
pheatmap(binary_matrix,show_rownames = F, show_colnames = F,clustering_method = 'ward.D')
rank_k <- 15
nmf_result <- nmf(binary_matrix, rank = rank_k, method = "brunet", nrun = 10)
W <- basis(nmf_result)  # Dimensions: 240 x rank_k
H <- coef(nmf_result)  # Dimensions: rank_k x 150

row.names(H)=seq(1,10,1)
pheatmap(H,show_colnames = F,cluster_rows = F)
colnames(W)=seq(1,10,1)
pheatmap(W,show_rownames = F,cluster_cols = F,clustering_method = 'ward.D',scale = 'column')

df_W = W%>%data.frame()%>%rownames_to_column('gene')%>%
  pivot_longer(cols=2:16,names_to = 'factor',values_to = 'loadings')
select_genes = df_W%>%filter(factor=='X10')%>%top_n(20)%>%pull(gene)
pheatmap(coef_mat[select_genes,],show_rownames = F)

df_H = H%>%t()%>%data.frame()%>%rownames_to_column('id')%>%
  pivot_longer(cols=2:16,names_to = 'factor',values_to = 'nmf_score')%>%
  mutate('region'=str_split_fixed(id,'_',3)[,1],
         'celltype'= str_split_fixed(id,'_',3)[,2],
         'subtype'=str_split_fixed(id,'_',2)[,2])


df_H$region = factor(df_H$region,levels=regions)
df_H%>%ggplot(.,aes(subtype,nmf_score,fill=factor))+geom_boxplot()+facet_wrap(~factor,scales = 'free_y')
df_H%>%ggplot(.,aes(region,nmf_score,fill=factor))+geom_boxplot()+facet_wrap(~factor,scales = 'free_y')

write.csv(as.data.frame(H),'/home/wyang114/output/transcriptional/nmf/upaDEG_all_10f/H_mat.csv',row.names = TRUE)
write.csv(as.data.frame(W),'/home/wyang114/output/transcriptional/nmf/upaDEG_all_10f/W_mat.csv',row.names = TRUE)
