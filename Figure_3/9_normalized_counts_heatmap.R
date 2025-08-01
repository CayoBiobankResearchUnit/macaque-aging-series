source('/home/wyang114/utils/plotting_theme.R')
library(monocle3)
library(hooke)
library(purrr)
library(tidyverse)
library(Matrix)
library(RColorBrewer)
require(ComplexHeatmap)
set.seed(2024)

###pooling all results from different cell type
folders =list.dirs("output/hooke_result/")
counts_all = NULL
for(folder in folders[c(2:9,14:length(folders))]){
  ct = basename(tail(folder, 1))
  print(ct)
  counts = read_csv(paste0("output/hooke_result/",ct,"/",ct,"_hooke_normalized_counts.csv"))[-1]
  counts_all <- rbind(counts_all,counts)
}
rm(counts)
write.csv(counts_all,'/home/wyang114/output/hooke_result_final/all_hooke_normalized_counts.csv')

counts_all=read.csv('/home/wyang114/output/hooke_result_final/cell-subtype_full_linear_ctrl-sex_counts.csv')[-1]
count_table = read.table(file = '/home/wyang114/output/hooke_result_final/count_table.csv',sep=',',row.names = 1,header = TRUE)

counts_old=read.csv('/home/wyang114/output/hooke_result_final/cell-subtype_full_linear_ctrl-sex_age>3_count.csv')[-1]

##plot loess
temp = color_map[color_map$identity %in% counts_all$region,]
col_vec = temp$color
names(col_vec) = temp$identity
select_region=c('ACC','dlPFC','IPP','M1','EC','HIP','NAc')
pdf('/home/wyang114/astr1_full.pdf',width=6,height=3)
counts_all%>%filter(celltype=='astrocytes',cell_group=='1')%>%
  filter(region %in% select_region)%>%
  ggplot(.,aes(x=age,y=num_cells))+geom_point(size=0.5)+
  geom_smooth(method='loess')+
  ylab('normalized count')+
  #scale_color_brewer(palette = 'Paired')+
  #scale_color_manual(values = col_vec)+
  ggtitle('astrocytes 1 full lifespan')+facet_wrap(~region,scale='free_y',ncol=4)
dev.off()
pdf('/home/wyang114/astr1_full_old.pdf',width=6,height=3)
counts_old%>%filter(celltype=='astrocytes',cell_group=='1')%>%
  filter(region %in% select_region)%>%
  ggplot(.,aes(x=age,y=num_cells))+geom_point(size=0.5)+
  geom_smooth(method='loess')+
  ylab('normalized count')+
  #scale_color_brewer(palette = 'Paired')+
  #scale_color_manual(values = col_vec)+
  ggtitle('astrocytes 1 (age>3)')+facet_wrap(~region,scale='free_y',ncol=4)
dev.off()

#try using 10 intervals
bins = 10
counts_all$age_bins = cut_interval(counts_all$age,n=bins)
mean_age = counts_all%>%group_by(age_bins)%>%summarise(mean_age=mean(age))

select_subtype = count_table%>%filter(n>=50)%>%pull('id')

temp = counts_all%>%
  group_by(cell_group,region,celltype,age_bins)%>%
  summarise(mean_age = mean(num_cells))%>%
  pivot_wider(.,names_from=age_bins,values_from = mean_age)%>%
  mutate(id = paste(region,celltype,cell_group,sep='_'))%>%
  filter(id %in% select_subtype)%>%
  filter(!(region=='MB' & celltype=='basket-cells'))
#First try z-score
data_norm = apply(temp[4:(4+bins-1)],1,scale)%>%t()%>%data.frame()
colnames(data_norm)= colnames(temp[4:(4+bins-1)])
row.names(data_norm) = temp$id
#data_norm = cbind(temp[1:3],data_norm)

#######Ward.D2 hierarchical clustering
library(pheatmap)
data_subset = data_norm
data_subset = data_norm[grep('astrocytes_0',row.names(data_norm)),]
global_min = min(data_norm)
global_max = max(data_norm)
breaks <- seq(global_min, global_max, length.out = 51)
pal = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(50)
pheatmap(t(data_subset),show_rownames = T, show_colnames = T,cluster_rows = F,breaks = breaks,
         clustering_method = "ward.D2",color = pal)


data_norm_subtract =apply(data_subset,2,function(x){x-data_subset[,1]})
global_min = -4
global_max = 4
breaks <- seq(global_min, global_max, length.out = 51)
pheatmap(t(data_norm_subtract),show_rownames = T, show_colnames = T,cluster_rows = F,breaks = breaks,
         clustering_method = "ward.D2",color = pal)



#######k-means clustering 
km.result = kmeans(data_norm_subtract,centers=4,nstart=25)
#km.result = kmeans(clean_df_all,centers=9,nstart=25)
cluster_samples <- lapply(1:max(km.result$cluster), function(cluster) {
  which(km.result$cluster == cluster)
})
g_list=NULL
for(i in 1:max(km.result$cluster)){
  #samples = cluster_samples[[i]][grep('opc-olig',names(cluster_samples[[i]]))]  #for subtype specific clustering
  summary=data_norm_subtract[cluster_samples[[i]],]%>%t()%>%as.data.frame()%>%rownames_to_column(var='age')%>%
    pivot_longer(cols=-age,names_to = 'id',values_to = 'z-score')
  summary$age = factor(summary$age,levels=colnames(data_norm_subtract))
  g <- ggplot(summary,aes(x=age,y=`z-score`,group=id))+geom_line()+theme_classic()+
    theme(legend.position = "none")+ggtitle(paste0('cluster_',i))
  if(i %in% c(1,3,6,7)){
    mean_color='pink'
  }else{mean_color='lightblue'}
  g=g + stat_summary(fun=mean, geom="line", lwd = 0.5, color=mean_color, aes(group=1))
  g_list[[i]]=g
}
grid.arrange(grobs = g_list, ncol = 3,top=textGrob("kmeans"))

#######loess (gives similar results as spline)
counts_all = counts_all%>%mutate(id=paste(region,celltype,cell_group,sep='_'))
ids = counts_all$id%>%unique()
subset = counts_all%>%filter(id==ids[12])
loess_fit <- loess(num_cells ~ age, data =subset, span = 0.75)
# Predict values on a grid of time points (optional, for smoother plotting)
time_grid <- seq(min(subset$age), max(subset$age), length.out = 10)
predicted_values <- predict(loess_fit, newdata = data.frame(age = time_grid))
# Plot the original data points
plot(subset$age, subset$num_cells, main = "LOESS Fit with Data Points",
     xlab = "Time Points", ylab = "Values", pch = 16, col = "black")
# Add the LOESS fit line
lines(time_grid, predicted_values, col = "blue", lwd = 2)

residuals <- residuals(loess_fit)
rss <- sum(residuals^2)
tss <- sum((subset$num_cells - mean(subset$num_cells))^2)
r_squared <- 1 - (rss / tss)
predicted

#######time-warping
require(dtw); require(dtwclust) # for time clustering
#data_norm = data_norm[2:10] #if excluding new-borns
k_clusters <- 9
clusts0 <- tsclust(data_norm, type='partitional',
                   k=4:9, # to check different values
                   distance='dtw_basic', centroid='pam', trace=T)
saveRDS(clusts0, '~/temp/dtw_clust_mean_age.Rds')
clusts0 = readRDS('~/temp/dtw_clust_mean_age.Rds')
clusts <- data.frame(name=row.names(data_norm), gene_cluster=clusts0[[6]]@cluster) %>%
  # for now using same previous order, but might need to reorder
  merge(data.frame(gene_cluster=1:k_clusters)) %>%
  identity()

clusts$region = str_split_fixed(clusts$name,'_',2)[,1]
clusts$subtype = str_split_fixed(clusts$name,'_',2)[,2]
df_linear <- read_csv("/home/wyang114/output/hooke_result/cell_abundance_subtype_full_linear_ctrl-sex_mash_age>3.csv")[-1]%>%
  data.frame()
df_all = left_join(df_linear,clusts,by=c('region','cell_subtype'='subtype'))
df_all = df_all%>%mutate(significance = case_when(age_mashr_lfsr<0.05 ~ 'significant',
                                                                age_mashr_lfsr>0.05 ~ 'not significant'))
df_all%>%
  mutate(gene_cluster=as.character(gene_cluster))%>%
  ggplot(.,aes(gene_cluster,fill=significance))+geom_bar(width=0.5)+
  ylab('num of region-subtypes')+xlab('k-means cluster')+
  coord_flip()


#metrics
names(clusts0) <- paste0('k', 4:9)
res_mat <- sapply(clusts0, cvi, type='internal')
saveRDS(res_mat, '~/temp/evaluation_metrics.Rds')
res_mat <- readRDS('~/temp/evaluation_metrics.Rds')
res_df <- res_mat %>% data.frame %>% tibble::rownames_to_column('metric') %>%
  pivot_longer(!metric, names_to='k') %>%
  mutate(k=as.integer(gsub('k', '', k)))
p <- res_df %>% filter(metric!='SF') %>% # 7 seems to work the best
  ggplot(aes(x=k, y=value)) + geom_point(size=0.4) + geom_line() +
  facet_wrap(.~metric, scales='free') +
  ylab('Metric value') + xlab('K')
ggsave('~/temp/k_metric_comparison.png',
       p + geom_vline(xintercept=4, linetype='dashed'), height=5, width=5.5)

#average pattern
some_colors <- RColorBrewer::brewer.pal(n=k_clusters, name='Paired')
names(some_colors) <- as.character(1:k_clusters)

# Try smooth fit - mgcv::gam for gam - generalized additive models
p1 <- data_norm %>% rownames_to_column('name')%>%
  left_join(.,clusts[1:2])%>%column_to_rownames('name')%>%
  group_by(gene_cluster)%>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))%>%
  pivot_longer(cols = 2:10,names_to = 'age_bins',values_to = 'mean_value')%>%
  left_join(.,mean_age[2:10,])%>%
  ggplot(aes(x=as.numeric(mean_age), y=mean_value, color=factor(gene_cluster),
             group=factor(gene_cluster))) + geom_line()+
  #geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs",k = 3),
  #            size=0.5, se=F) +
  scale_color_manual("Cluster", values=some_colors) +
  guides(color=guide_legend(override.aes=list(size=1.5), keywidth=unit(5, 'mm'))) +
  theme(legend.position='top', panel.grid.major.x=element_line(size=0.1,
                                                               color=scales::alpha('black', 0.2))) +
  facet_wrap(~gene_cluster)+
  #scale_x_continuous(breaks=5:18, limits=c(4.5, 18.5), expand=c(0,0)) +
  xlab('mean_age(yrs)') + ylab("scaled abundance") + NULL
p1
ggsave('figures/mesoderm_cluster_fit.pdf', p1, height = 37,
       width = 83, units='mm', useDingbats=FALSE)

#make heatmap
some_colors2 <- some_colors
names(some_colors2) <- as.character(1:k_clusters)
ha <- rowAnnotation(Clusters=as.character(clusts$gene_cluster),
                    col = list(Clusters=some_colors2), show_annotation_name = F, show_legend=F)
p <- Heatmap(as.matrix(data_norm)[clusts$name,1], cluster_columns=F,
             show_column_names=F, 
             column_names_gp = grid::gpar(fontsize = 1, fontfamily="Helvetica"),
              row_split=as.character(clusts$gene_cluster),
             row_gap = unit(1, "mm"), left_annotation=ha,show_heatmap_legend = F,
             use_raster = F, show_row_names=F, cluster_rows=F)+
    Heatmap(as.matrix(data_norm)[clusts$name,2:10], cluster_columns=F,
            show_column_names=F, 
            column_names_gp = grid::gpar(fontsize = 1, fontfamily="Helvetica"),
            name='z-score', row_split=as.character(clusts$gene_cluster),
            row_gap = unit(1, "mm"), 
            use_raster = F, show_row_names=F, cluster_rows=F)
             
mm_to_inches_conversion <- 0.0593701
cairo_pdf('~/temp/data_norm_clu9_heatmap_age<3.pdf',
          height=70*mm_to_inches_conversion, width=70*mm_to_inches_conversion)
print(p); dev.off()
