##Plot pseudotime for microglia
##Wei Yang
##Input: Monocle object with pseudotime calculation from pseudotime_microglia.R

library(reticulate)
use_condaenv(condaenv = '~/.conda/envs/TACCO_env/', required = TRUE)
library("anndata")
library(Seurat)
library(monocle3)
library(tidyverse)
library(dplyr)
library(data.table)


###############################################
#Create Seurat object and find markers for partitions
data_dir = "/scratch/nsnyderm/u01/intermediate_files/cell-class_h5ad_update/"
files = list.files(data_dir)
file = files[grep('Res1_microglia_new.h5ad',files)]
data <- read_h5ad(paste0(data_dir,file))
var_info <- data$var
var_info  %>% write.table('/scratch/wyang114/data/var_info.txt',
                       quote=F, row.names=F, sep='\t')
obj <- CreateSeuratObject(counts = t(data$X), meta.data = data$obs)
obj = obj[,!obj$louvain %in% c(2,7,11)]
obj<- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)

cds = readRDS('/scratch/wyang114/outputs/pseudotime/pseudotime_microglia_2partitions.rds')
all(row.names(pData(cds))==row.names(obj@meta.data))
obj = obj[,row.names(pData(cds))]
obj@active.ident = partitions(cds)
###############################################

###############################################
#Take 5000 variable genes in partition 1 
#Bin time/pseudotime to 100 bins and calculate mean expression for each bin
cds = readRDS('/scratch/wyang114/outputs/pseudotime/pseudotime_microglia_2partitions.rds')
pData(cds)%>%data.frame%>%group_by(partition,louvain)%>%summarise(n=n())%>%ungroup()%>%
  group_by(partition)%>%
  mutate(freq=n/sum(n))%>%ggplot(.,aes(x=as.character(partition),y=freq,fill =louvain))+
  geom_bar(stat='identity')+scale_fill_brewer(palette = 'Paired')+theme_classic()+
  xlab('group')

obj$partition = partitions(cds)
#obj$pseudotime_p1 = pData(cds)$pseudotime_p1
#obj$pseudotime_p2 = pData(cds)$pseudotime_p2
obj$pseudotime = pData(cds)$pseudotime
obj@meta.data[obj@meta.data$partition==1,'pseudotime.bin'] = cut(obj@meta.data[obj@meta.data$partition==1,'pseudotime'], breaks=50)
obj@meta.data%>%data.frame()%>%filter(partition==1)%>%
  group_by(pseudotime.bin,partition)%>%summarise(mean_age=mean(age))%>%
  ggplot(.,aes(pseudotime.bin,mean_age,group=partition))+
  geom_point()+
  #geom_line()+
  theme_classic()+
  ylim(c(5,12))+
  theme(axis.text.x = element_text(angle = 90))

subset = obj[,obj$partition==1]
subset = FindVariableFeatures(subset,nfeatures = 5000)
subset = ScaleData(subset)
saveRDS(subset,'/scratch/wyang114/outputs/pseudotime/microglia_p1_seurat_v2.RDS')
var_genes <- VariableFeatures(subset)
mat = subset@assays[['RNA']]$data[var_genes,]

rna <- subset@meta.data %>% rownames_to_column('cell')%>%select(cell, pseudotime) %>% 
  merge(data.frame(t(mat), check.names=F),
        by.x='cell', by.y='row.names') %>%
  pivot_longer(all_of(var_genes)) %>%
  mutate(new_bin=cut_number(pseudotime, 100)) %>%
  data.frame
rna %>% saveRDS('/scratch/wyang114/outputs/pseudotime/microglia_5kvargenes_p1_v2.Rds')

 rna2 <- rna %>%
  group_by(name, new_bin) %>%
  summarize(mean_value=mean(value),
            mean_ptime=mean(pseudotime), .groups='drop') %>%
  data.frame
rna2 %>% saveRDS(paste0('/scratch/wyang114/outputs/pseudotime/microglia_5kvargenes_mean_p1_v2.Rds'))

# Check the standard deviation of mean expression across bins
rna_wide <- rna2 %>%
  pivot_wider(id_cols=c('new_bin', 'mean_ptime'),
              names_from=name, values_from=mean_value)
df_mat <- rna_wide %>% select(-new_bin, -mean_ptime) %>% as.matrix
rownames(df_mat) <- rna_wide$new_bin
new_vars <- apply(df_mat, 2, sd)
###############################################

###############################################
#Perform smoothing on already averaged data by a moving time average with 15 windows
rna <- readRDS(paste0('/scratch/wyang114/outputs/pseudotime/microglia_5kvargenes_mean_p1_v2.Rds'))
n_adjacent_windows <- 15 # 15 looks good
total_windows <- sort(unique(rna$new_bin))
n_windows <- length(total_windows)
start_edge_idx <- 1:n_adjacent_windows
end_edge_idx <- (length(total_windows)-(n_adjacent_windows-1)):length(total_windows)
start_edge_bins <- total_windows[start_edge_idx]
end_edge_bins <- total_windows[end_edge_idx]
used_srnas <- rbindlist(lapply(1:n_windows, function(i) {
  print(paste0('working on bin: ', total_windows[i]))
  
  # decide which windows to average over
  if (i %in% start_edge_idx) {# beginning edge
    place <- 'first'
    start_idx <- 1; end_idx <- i+n_adjacent_windows
  }
  if (i %in% end_edge_idx) { # end edge
    place <- 'end'
    start_idx <- i-n_adjacent_windows; end_idx <- n_windows
  }
  if (!(i %in% start_edge_idx) & !(i %in% end_edge_idx)) { # in between
    place <- 'middle'
    start_idx <- i-n_adjacent_windows; end_idx <- i+n_adjacent_windows
  }
  
  print(paste0('averaging over: ', start_idx, ' to ', end_idx))
  central_bin_ptime <- mean((rna %>% filter(new_bin == total_windows[i]))$mean_ptime)
  # otherwise return the average values across the bins
  return(rna %>% filter(new_bin %in% total_windows[start_idx:end_idx]) %>% 
           group_by(name) %>%
           summarize(mean_value=mean(mean_value, trim=0.1),
                     mean_ptime=central_bin_ptime,
                     .groups='drop') %>% 
           mutate(new_bin=total_windows[i]) %>%
           select(name, new_bin, mean_value, mean_ptime))
}))
###############################################

###############################################
#Select top variable genes and scaling
rna_wide <- used_srnas %>%
  pivot_wider(id_cols=c('new_bin', 'mean_ptime'),
              names_from=name, values_from=mean_value)
df_mat <- rna_wide %>% select(-new_bin, -mean_ptime) %>% as.matrix
rownames(df_mat) <- rna_wide$new_bin

#Select top2000 genes with large stdev along the time axis
new_vars <- apply(df_mat, 2, sd)
#temp = data.frame('name'=colnames(df_mat),'sd'=new_vars)
#temp$external_gene_name = var_info[temp$name,'external_gene_name']
df_mat <- scale(df_mat[,rev(order(new_vars))[1:2000]])


# z-score scaling of expression for each gene across time
scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / stats::sd(x, na.rm=TRUE)
}
rna_wide <- used_srnas %>%
  filter(name %in% names(rev(sort(new_vars)))[1:2000]) %>%
  group_by(name) %>% mutate(mean_value=scale_this(mean_value)) %>% ungroup %>%
  pivot_wider(id_cols=c('new_bin', 'mean_ptime'),
              names_from=name, values_from=mean_value)
df_mat <- rna_wide %>% select(-new_bin, -mean_ptime) %>% as.matrix
rownames(df_mat) <- rna_wide$new_bin

###############################################
##### Eucleadian distance
k_clusters = 7
km.result = kmeans(t(df_mat),centers=k_clusters,nstart=25)
clusts <- data.frame(name=colnames(df_mat), gene_cluster=km.result$cluster[colnames(df_mat)])
clusts %>% write.table('/scratch/wyang114/outputs/pseudotime/microglia_p1_gene_clusters_kmeans7_2k_v2.txt',
                       quote=F, row.names=F, sep='\t')
clusts <- fread('/scratch/wyang114/outputs/pseudotime/microglia_p1_gene_clusters_kmeans7_2k_v2.txt');k_clusters = 7
###############################################
##### tyring slightly fancier clustering strategy, dynamic time warping
require(dtw); require(dtwclust) # for time clustering
k_clusters <- 7 #determined from k_metric_comparison
clusts0 <- tsclust(t(df_mat), type='partitional',
                   k=7, # to check different values
                   distance='dtw_basic', centroid='pam', trace=T)
saveRDS(clusts0, '/scratch/wyang114/outputs/pseudotime/dtw_clust.Rds') # saved already
clusts0 <- readRDS('/scratch/wyang114/outputs/pseudotime//dtw_clust.Rds')
clusts <- data.frame(name=colnames(df_mat), gene_cluster=clusts0[[1]]@cluster) %>%
  # for now using same previous order, but might need to reorder
  #merge(data.frame(gene_cluster=1:k_clusters, new_clusters=c(3,1,2,4))) %>%
  identity()
clusts %>% write.table('/scratch/wyang114/outputs/pseudotime/micrglia_p1_gene_clusters.txt',
                       quote=F, row.names=F, sep='\t')
clusts <- fread('/scratch/wyang114/outputs/pseudotime/micrglia_p1_gene_clusters.txt');

# Make plots of clustering metrics with the dtw_clust object
names(clusts0) <- paste0('k_', 4:10)
res_mat <- sapply(clusts0, cvi, type='internal')
saveRDS(res_mat, '/scratch/wyang114/outputs/pseudotime/evaluation_metrics.Rds')
res_mat <- readRDS('/scratch/wyang114/outputs/pseudotime/evaluation_metrics.Rds')
res_df <- res_mat %>% data.frame %>% tibble::rownames_to_column('metric') %>%
  pivot_longer(!metric, names_to='k') %>%
  mutate(k=as.integer(gsub('k_', '', k)))
p <- res_df %>% filter(metric!='SF') %>% # didn't work well
  ggplot(aes(x=k, y=value)) + geom_point(size=0.4) + geom_line() +
  facet_wrap(.~metric, scales='free') +
  ylab('Metric value') + xlab('K')
ggsave('~/figures/k_metric_comparison.png',
       p + geom_vline(xintercept=4, linetype='dashed'), height=5, width=5.5)
# saving as supp figure
ggsave('~/figures/k_metric_comparison.pdf',
       p + geom_vline(xintercept=4, linetype='dashed', size=0.3),
       height = 100, width = 100,
       units='mm', useDingbats=FALSE)
###############################################

###############################################
##Plot expression of clusters (z-score scaling. Shall we try min-max instead of z-score?)
rnasc <- used_srnas %>% # sliding window values
  filter(name %in% colnames(df_mat)) %>% group_by(name) %>%
  mutate(mean_value=scale_this(mean_value)) %>% ungroup %>%
  merge(clusts) %>% identity()

# # Let's make a version with an averged pattern
some_colors <- RColorBrewer::brewer.pal(n=k_clusters, name='Set1')
names(some_colors) <- as.character(1:k_clusters)
# Try smooth fit - mgcv::gam for gam - generalized additive models
p1 <- rnasc %>% 
  ggplot(aes(x=mean_ptime, y=mean_value, color=factor(gene_cluster),
             group=factor(gene_cluster))) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
              size=0.5, se=F) +
  scale_color_manual("Cluster", values=some_colors) +
  guides(color=guide_legend(override.aes=list(size=1.5), keywidth=unit(5, 'mm'))) +
  theme(legend.position='top', panel.grid.major.x=element_line(size=0.1,
                                                               color=scales::alpha('black', 0.2))) +
  #scale_x_continuous(breaks=5:18, limits=c(4.5, 18.5), expand=c(0,0)) +
  xlab('pseudotime') + ylab("Scaled expression") + NULL

var_info = read.table('/scratch/wyang114/data/var_info.txt', sep='\t',header = TRUE)
genes = var_info%>%filter(external_gene_name %in% c('CPM','SPP1'))%>%pull(ensembl_gene_id)
p2 <- rnasc %>% filter(name %in% genes)%>%
  ggplot(aes(x=mean_ptime, y=mean_value,group = name, color=name)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
              size=0.5, se=F) +
  #scale_color_manual("Cluster", values=some_colors) +
  guides(color=guide_legend(override.aes=list(size=1.5), keywidth=unit(5, 'mm'))) +
  theme(legend.position='top', panel.grid.major.x=element_line(size=0.1,
                                                               color=scales::alpha('black', 0.2))) +
  xlab('pseudotime') + ylab("Scaled expression") + NULL

# make a heatmap
require(patchwork); require(ComplexHeatmap)
some_colors2 <- some_colors
names(some_colors2) <- as.character(1:k_clusters)
ha <- rowAnnotation(Clusters=as.character(clusts$gene_cluster),
                    col = list(Clusters=some_colors2), show_annotation_name = F, show_legend=F)
some_genes <- c('CPM','SPP1','APOE','TPRG1','GPNMB')
var_info <- fread('/scratch/wyang114/data/var_info.txt')
var_info <- var_info%>%column_to_rownames('ensembl_gene_id')
clusts$external_gene_name = var_info[clusts$name,'external_gene_name']
some_genes <- some_genes[some_genes %in% clusts$external_gene_name]
some_genes_idx <- sapply(some_genes, function(g) which(g == clusts$external_gene_name))
ha2 <- rowAnnotation(
  foo = anno_mark(at = some_genes_idx, labels = some_genes,
                  link_width = unit(4, "mm"),
                  labels_gp = grid::gpar(fontfamily='Helvetica', fontsize=6)))
p <- Heatmap(t(df_mat)[clusts$name,], cluster_columns=F,
             show_column_names=F, right_annotation = ha2,
             column_names_gp = grid::gpar(fontsize = 1, fontfamily="Helvetica"),
             name='Expression', row_split=as.character(clusts$gene_cluster),
             row_gap = unit(1, "mm"), left_annotation=ha,
             use_raster = F, show_row_names=F, cluster_rows=F)
mm_to_inches_conversion <- 0.0393701
cairo_pdf('~/figures/mgl_p1_heatmap_pseudotime_kmeans7_updated_pseudotime.pdf',
          height=100*mm_to_inches_conversion, width=110*mm_to_inches_conversion)
print(p); dev.off()


#plot genes
subset = readRDS('/scratch/wyang114/outputs/pseudotime/microglia_p1_seurat.RDS')
dim(subset)
temp = subset@meta.data
temp2 = pData(cds)%>%data.frame()
temp[!row.names(temp) %in% row.names(temp2),]%>%summary()
subset = subset[,colnames(subset) %in% row.names(temp2)] #remove cells that are re-assigned
dim(subset)
subset$update_louvain = temp2[row.names(subset@meta.data),'louvain']
var_info = read.table('/scratch/wyang114/data/var_info.txt', sep='\t',header = TRUE)
genes = var_info%>%filter(external_gene_name %in% c('CPM','SPP1'))%>%pull(ensembl_gene_id)

cor(subset@meta.data$pseudotime_p1,subset@meta.data$age,method = 'spearman')

subset@meta.data%>%data.frame()%>%filter(partition==1)%>%
  mutate(pseudotime.bin=as.numeric(cut(pseudotime_p1,breaks=50,labels=seq(1,50,1))))%>%
  group_by(pseudotime.bin)%>%summarise(mean_age=mean(age))%>%
  ggplot(.,aes(pseudotime.bin,mean_age))+
  geom_point()+
  theme_classic()+
  ylim(c(0,21))+
  theme(axis.text.x = element_text(angle = 90))+
  ylab('mean_age(years)')

subset@meta.data%>%data.frame()%>%filter(partition==1)%>%
  mutate(pseudotime.bin=as.numeric(cut(pseudotime_p1,breaks=10,labels=seq(1,10,1))))%>%
  group_by(pseudotime.bin,update_louvain)%>%summarise(n=n())%>%ungroup%>%
  group_by(pseudotime.bin)%>%
  mutate(freq=n/sum(n))%>%
  ggplot(.,aes(pseudotime.bin,y=freq,fill=update_louvain))+geom_bar(stat='identity')+
  scale_fill_brewer(palette = 'Paired')+theme_classic()+
  xlab('group')



rna_df <- subset@meta.data %>%rownames_to_column('cell')%>%
  select(cell, pseudotime_p1) %>%
  merge(data.frame(t(subset@assays[['RNA']]$scale.data[genes,]),
                   check.names=F), by.x='cell', by.y='row.names')%>%
  mutate(CPM=ENSMMUG00000000946,SPP1=ENSMMUG00000008793)%>%identity()%>%
  mutate(region=str_split_fixed(cell,'-',4)[,4])%>%select(cell,pseudotime_p1,CPM,SPP1,region)
min_max <- function(x) {(x-min(x))/(max(x)-min(x))}
df_mean <- rna_df %>% select(pseudotime_p1, CPM,SPP1) %>% 
  pivot_longer(!c('pseudotime_p1')) %>% group_by(name) %>%
  mutate(bin=cut_interval(pseudotime_p1, 100)) %>% ungroup %>% 
  group_by(name, bin) %>%
  summarize(mean_age=mean(pseudotime_p1),
            mean_value=mean(value), .groups='drop') %>%
  group_by(name) %>% 
  mutate(mean_value=min_max(mean_value)) %>% ungroup
df_mean %>%
  ggplot(aes(x=mean_age, y=mean_value,
             color=name, group=name)) +
  geom_smooth(method='loess', se=F, size=0.5) +
  scale_color_manual("Gene", values=c('#D95F02', '#7570B3')) +
  guides(color=guide_legend(override.aes = list(size = 0.8)),
         linetype=guide_legend(override.aes = list(size = 0.5))) +
  theme_classic()+
  xlab('pseudotime_p1') + ylab('Scaled value')

#plot regions
rna_df%>%ggplot(.,aes(x=pseudotime_p1,fill=region))+geom_density(alpha=0.5)+facet_wrap(~ region)
rna_df%>%select(pseudotime_p1,region)%>%
  mutate(bin=cut_interval(pseudotime_p1, 50))%>%
  group_by(bin,region)%>%summarise(count=n())%>%
  group_by(bin)%>%
  mutate(freq=count/sum(count))%>%
  ggplot(.,aes(x=bin,y=freq,fill=region))+geom_bar(stat = 'identity')+theme_classic()
