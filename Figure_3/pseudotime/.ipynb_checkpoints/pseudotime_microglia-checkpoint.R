library(reticulate)
use_condaenv(condaenv = '~/.conda/envs/TACCO_env/', required = TRUE)
library("anndata")
library(Matrix)
library(tidyverse)
library(dplyr)
library(monocle3)
library(ggridges)
theme_set(theme_minimal())
source('/home/wyang114/utils/plotting_theme.R')

data_dir = "/scratch/nsnyderm/u01/intermediate_files/cell-class_h5ad_update//"
files = list.files(data_dir)
file = files[grep('Res1_microglia_new.h5ad',files)]
data <- read_h5ad(paste0(data_dir,file))
expression_matrix = as(data$X, "CsparseMatrix")  #convert dgRmatrix from h5ad to dgCmatrix for monocle input
cell_metadata = data$obs
gene_annotation = data$var
cds <- new_cell_data_set(t(expression_matrix),
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
#Select microglia cells only
cds = cds[,!pData(cds)$louvain %in% c(11,2,7)]

cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds,umap.min_dist = 0.1,umap.n_neighbors = 15L,cores = 4)
cds <- cluster_cells(cds, reduction_method = "UMAP", resolution = 1e-5) #partitions
cds <- learn_graph(cds)
pData(cds)$partition = partitions(cds)
#pData(cds)$age_bins = cut(pData(cds)$age,breaks = c(0,3,9,15,21))  #setting to vertex with the most cells from young animals
get_earliest_principal_node <- function(cds, time_bin=c(5,9)){
  cell_ids_1 <- which(pData(cds)[, "louvain"] %in% time_bin  & 
                      pData(cds)[, "partition"] == 1)
  cell_ids_2 <- which(pData(cds)[, "louvain"] %in% time_bin  & 
                        pData(cds)[, "partition"] == 2)
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    c(igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                    (which.max(table(closest_vertex[cell_ids_1,]))))],
      igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                    (which.max(table(closest_vertex[cell_ids_2,]))))]
    )
  root_pr_nodes
}
cds <- order_cells(cds,root_pr_nodes=get_earliest_principal_node(cds))

plot_cells(cds,
           show_trajectory_graph = FALSE,
           color_cells_by = "pseudotime",
           #group_cells_by = "louvain",
           label_cell_groups=TRUE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=10,
           rasterize = TRUE)

pseudotime_tbl <-
  monocle3::pseudotime(cds) %>%
  tibble::enframe("sample_id", "pseudotime")
all(pseudotime_tbl$sample_id==row.names(pData(cds)))
pData(cds)$pseudotime = pseudotime_tbl$pseudotime

#plot mean time on binned pseudotime
pData(cds)$partition = partitions(cds)
pData(cds)$pseudotime.bin = NA
pData(cds)[pData(cds)$partition==1,'pseudotime.bin'] = cut(pData(cds)[pData(cds)$partition==1,
                                                                      'pseudotime'], breaks=50)
pData(cds)[pData(cds)$partition==2,'pseudotime.bin'] = cut(pData(cds)[pData(cds)$partition==2,
                                                                      'pseudotime'], breaks=50)
pData(cds)%>%data.frame()%>%
  group_by(pseudotime.bin,partition)%>%summarise(mean_age=mean(age))%>%
  ggplot(.,aes(pseudotime.bin,mean_age,group=partition,color=partition))+
  geom_line()+theme_classic()+
  ylim(c(5,12))+
  theme(axis.text.x = element_text(angle = 90))

pData(cds)%>%data.frame()%>%filter(partition==1)%>%
  mutate(pseudotime.bin=as.numeric(cut(pseudotime,breaks=50,labels=seq(1,50,1))))%>%
  group_by(pseudotime.bin)%>%summarise(mean_age=mean(age))%>%
  ggplot(.,aes(pseudotime.bin,mean_age))+
  geom_point()+
  theme_classic()+
  ylim(c(0,21))+
  theme(axis.text.x = element_text(angle = 90))+
  ylab('mean_age(years)')
subset = cds[,cds$partition==1]
cor(subset$pseudotime,subset$age,method = 'spearman')

pData(cds)%>%data.frame()%>%filter(partition==1)%>%
  group_by(pseudotime.bin,louvain)%>%summarise(louvain_count=n())%>%
  ggplot(.,aes(pseudotime.bin,louvain_count,group=louvain,color=louvain))+
  geom_line()+theme_classic()+
  #ylim(c(5,12))+
  theme(axis.text.x = element_text(angle = 90))+facet_wrap(~louvain)



saveRDS(cds,'/scratch/wyang114/outputs/pseudotime/pseudotime_microglia_2partitions.rds')

cds = readRDS('/scratch/wyang114/outputs/pseudotime/pseudotime_microglia_2partitions.rds')
cds <- preprocess_cds(cds, num_dim = 50,norm_method='log')
mean_age_df = read.csv('/scratch/wyang114/outputs/pseudotime/microglia_mean_age.csv',row.names = 1)
row.names(mean_age_df)=mean_age_df$self_index
pData(cds)$nn_mean_age = mean_age_df[colnames(cds),'mean_knn_age']

#hist(pData(cds)$nn_mean_age)
#hist(pData(cds)$age,col=rgb(1,0,0,0.5),add=T)
#pData(cds)%>%data.frame%>%ggplot(.,aes(x=nn_mean_age,fill=region))+geom_histogram()+
#  scale_fill_brewer(palette = 'Set3')+scale_y_log10()

#Test with module score
# Subset expression data
rowData(cds) = rowData(cds)[1:3]
pick_genes <- c("APOE","GPNMB","CPM","SPP1")
homeo_genes <- c('HEXB', 'CST3', 'CX3CR1', 'CTSD', 'CSF1R', 'CTSS', 'SPARC', 
                 'TMSB4X', 'P2RY12', 'C1QB')
tremi <- c('TYROBP', 'CTSB', 'APOE', 'B2M', 'FTH1')

#normalize matrix
mat = exprs(cds)
norm_mat = t(t(mat) / colSums(mat)) * 10000
norm_mat = log(norm_mat+1) #sparse to dense

select_genes = homeo_genes
rowData(cds)$gene_short_name = rowData(cds)$external_gene_name
gene_id = rownames(cds)[rowData(cds)$external_gene_name %in% select_genes]
cell_id = colnames(cds)[colData(cds)$partition==1]
select_genes_mat <- norm_mat[gene_id,cell_id]
row.names(select_genes_mat) = rowData(cds)[gene_id,'external_gene_name']
select_gene_meta = pData(cds)[cell_id,]

plot_data <- data.frame(
  #pseudotime.bin = select_gene_meta$pseudotime.bin
  pseudotime.bin =cut(select_gene_meta$nn_mean_age,breaks=10)
  #pseudotime.bin = cut(select_gene_meta$age,breaks=10)
)
for(gene in select_genes){
  expression <- select_genes_mat[gene, ]
  plot_data[[gene]] = expression
}
plot_data%>%ggplot(.,aes(x=pseudotime.bin))+geom_histogram(stat='count')+
  scale_y_log10()
plot_data%>%
  pivot_longer(cols =1:length(select_genes)+1, names_to = 'gene',values_to = 'expression')%>%
  group_by(pseudotime.bin,gene)%>%
  summarise(mean_expr=mean(expression))%>%
  ggplot(.,aes(x=pseudotime.bin,y=mean_expr))+
  geom_point()+
  #geom_smooth(method = "loess", color = "blue", size = 1) +  # Smoothed line
  theme_classic() +theme(axis.text.x = element_blank())+
  labs(x = "pseudotime", y = "Expression", title = 'homeostatic')+facet_wrap(~gene,scale='free_y')


bins=50
pData(cds)$pseudotime.bin = NA
pData(cds)[pData(cds)$partition==1,'pseudotime.bin'] = cut(pData(cds)[pData(cds)$partition==1,
                                                                      'pseudotime'], breaks=bins)
pData(cds)[pData(cds)$partition==2,'pseudotime.bin'] = cut(pData(cds)[pData(cds)$partition==2,
                                                                      'pseudotime'], breaks=bins)
pData(cds)%>%data.frame()%>%
  mutate(pseudotime.bin = cut(nn_mean_age,breaks=bins))%>%
  group_by(pseudotime.bin,partition)%>%summarise(mean_age=mean(age))%>%
  ggplot(.,aes(pseudotime.bin,mean_age,group=partition,color=partition))+
  geom_line()+theme_classic()+
  #ylim(c(5,12))+
  ggtitle(sprintf('bins:%s',bins))+
  theme(axis.text.x = element_text(angle = 90))


# #Use umap coordinates from scanpy
# reducedDim(cds, type = "PCA") = data$obsm$X_pca
# cds@reduce_dim_aux[['PCA']][['model']][['prop_var_expl']] <- data$uns$pca$variance_ratio
# cds@reduce_dim_aux[['PCA']][['model']][['svd_v']] <- data$uns$pca$variance
# plot_pc_variance_explained(cds)
# cds@int_colData@listData$reducedDims$UMAP <- data$obsm$X_umap
# plot_cells(cds)
# cds$clusters <- cds$louvain
# cds <- cluster_cells(cds, reduction_method = "UMAP", resolution = 1e-5)
# rownames(cds@principal_graph_aux$UMAP$dp_mst) <- NULL
# colnames(cds@int_colData@listData$reducedDims$UMAP) <- NULL
# plot_cells(cds, color_cells_by = "partition")
# plot_cells(cds, color_cells_by = "louvain")


