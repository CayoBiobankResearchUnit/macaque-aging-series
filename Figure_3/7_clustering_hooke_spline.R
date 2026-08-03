#' @file 7_clustering_hooke_spline.R
#' @author Wei Yang
#'
#' @description
#' Generates panels for Figure 3I
#'
#' @usage
#' # Open in RStudio and execute lines sequentially (Cmd+Enter / Ctrl+Enter)
#' # Or render interactively:
#' source("7_clustering_hooke_spline.R", echo = TRUE)
NULL

source("utils/plotting_theme.R")
library(tidyverse)
library(dplyr)
library(gridExtra)
library(grid)
library(readr)
library(pheatmap)

###pooling all results from different cell type with >=50 cells without infants (1131)
folders =list.dirs("output/hooke_result/")
clean_df_all = NULL
for(folder in folders[c(2:9,14:length(folders))]){
  ct = basename(tail(folder, 1))
  print(ct)
  pred <- read_csv(paste0("output/hooke_result/",ct,"/",ct,"_all_prediction_bootstrap_update.csv"))
  pred$id <- paste0(pred$region,'_',ct,'_',pred$cell_group)
  df <- pred%>%select(log_abund,id,age)%>%
    pivot_wider(id_cols=age,names_from =id,values_from = log_abund)%>%
    tibble::column_to_rownames('age')
  zscore_df <- df %>%
    mutate(across(everything(), function(x){(x - mean(x)) / sd(x)}))
  clean_df <- select(zscore_df, where(~all(!is.na(.))))  #check cases with NA
  clean_df_all = rbind(clean_df_all,t(clean_df))
}
###clustering non-linear trajectories (1374)
#Remove groups with cell counts <100
df_linear <- read_csv("output/hooke_result_final/cell-subtype_full_linear_ctrl-sex_mash_age>3.csv")[-1]%>%
  data.frame()
df_linear$id = paste(df_linear$region,df_linear$cell_subtype,sep='_')
select_id = df_linear%>%filter(count>=100)%>%pull(id)
clean_df_all = clean_df_all[select_id,] #1045
#Remove msn in HIP 
clean_df_all = clean_df_all[-grep('HIP_medium-spiny-neurons',row.names(clean_df_all)),] #1039
clean_df_all = as.data.frame(clean_df_all)
write_csv(clean_df_all,file = '/home/wyang114/output/hooke_result_final/clean_zscore_df_all_spline_boostrap_update.csv')
clean_df_all = read.csv('/home/wyang114/output/hooke_result_final/clean_zscore_df_all_spline_boostrap_update.csv',row.names = 1)
colnames(clean_df_all) = seq(0,21,1)

#k-means clustering and plotting, centers could be 2:50
set.seed(2024)

km.result = kmeans(clean_df_all,centers=9,nstart=25)
cluster_samples <- lapply(1:max(km.result$cluster), function(cluster) {
  which(km.result$cluster == cluster)
})
g_list=NULL
for(i in 1:max(km.result$cluster)){
  summary=clean_df_all[cluster_samples[[i]],]%>%t()%>%as.data.frame()%>%rownames_to_column(var='age')%>%
    pivot_longer(cols=-age,names_to = 'id',values_to = 'z-score')
  summary$age = as.numeric(summary$age)
  g <- ggplot(summary,aes(x=age,y=`z-score`,group=id))+geom_line()+theme_classic()+
    theme(legend.position = "none")+ggtitle(paste0('cluster_',i))+ylim(c(-3,3))
  if(i %in% c(1,3,5)){
    mean_color='purple'
  }else{mean_color='green'}
  g=g + stat_summary(fun=mean, geom="line", lwd = 0.5, color=mean_color, aes(group=1))
  g_list[[i]]=g
}
grid_output = grid.arrange(grobs = g_list, ncol = 9,top=textGrob("kmeans"))
pdf('/home/wyang114/output/hooke_result_final/non-linear_trajectory.pdf',width = 12, height=1.5)  # Specify file name and size
grid.draw(grid_output)  # Draw the grid into the PDF
dev.off()

clean_df_meta = data.frame('id' = row.names(clean_df_all), 
                           'region' = str_split_fixed(row.names(clean_df_all),'_',3)[,1],
                           'cell_class' = str_split_fixed(row.names(clean_df_all),'_',3)[,2],
                           'cell_subtype' = str_split_fixed(row.names(clean_df_all),'_',2)[,2]
                           )%>%column_to_rownames('id')
clean_df_meta$group = ifelse(clean_df_meta$cell_class %in% c('astrocytes','ependymal-cells',
                                                           'microglia','opc-olig','vascular-cells'),'non-neurons','neurons')
for(i in 1:max(km.result$cluster)){
  clean_df_meta[cluster_samples[[i]],'cluster']=i
}

#clean_df_meta$cluster = factor(clean_df_meta$cluster,levels=c(1,3,8,9,5,7,2,6,4))
clean_df_meta%>%group_by(cluster)%>%summarise(n=n())%>%mutate(cluster=cluster)%>%
  ggplot(.,aes(cluster,n))+geom_bar(stat='identity',width=0.5,
                                                            fill=c(rep('lightblue',5),rep('pink',4)))+
  ylab('num of region-subtypes')+xlab('k-means cluster')
write.csv(clean_df_meta,'~/output/hooke_result_final/clean_df_meta_update_age>3.csv',row.names = TRUE)


#Plotting individuals
pattern = 'astrocytes_1'
subset = clean_df_all[grep(pattern,row.names(clean_df_all)),]
subset%>%rownames_to_column('id')%>%pivot_longer(cols = 2:dim(clean_df_all)[2]+1, names_to = 'time',values_to = 'count')%>%
  mutate(region=str_split_fixed(id,'_',2)[,1])%>%
  ggplot(.,aes(as.numeric(time),count,group=region,color=region))+geom_line()+ylim(c(-2,2))+
  theme_classic()


#add information about linear changes
clean_df_meta = read.csv('~/output/hooke_result_final/clean_df_meta_update_age>3.csv',row.names = 1)
df_linear <- read.csv("output/hooke_result_final/cell-subtype_full_linear_ctrl-sex_mash_age>3.csv",row.names=1)%>%
  data.frame()
row.names(df_linear) = paste(df_linear$region,df_linear$cell_subtype,sep = '_')
row.names(clean_df_meta) = paste(clean_df_meta$region,clean_df_meta$cell_subtype,sep = '_')
clean_df_meta$age_mashr_lfsr = df_linear[row.names(clean_df_meta),'age_mashr_lfsr']
clean_df_meta$age_mashr_pm = df_linear[row.names(clean_df_meta),'age_mashr_pm']
clean_df_meta$age_mashr_psd = df_linear[row.names(clean_df_meta),'age_mashr_psd']
clean_df_meta$count = df_linear[row.names(clean_df_meta),'count']
clean_df_meta = clean_df_meta%>%mutate(significance = case_when(age_mashr_lfsr<0.05 ~ 'significant',
                                                                age_mashr_lfsr>0.05 ~ 'not significant'))
clean_df_meta%>%filter(!is.na(significance))%>%
  mutate(cluster=as.character(cluster))%>%
  ggplot(.,aes(cluster,fill=significance))+geom_bar(width=0.5)+
  scale_fill_manual(values = c('not significant'='grey','significant'='pink'))+
  ylab('num of region-subtypes')+xlab('k-means cluster')
ggsave('/home/wyang114/output/hooke_result_final/non-linear_significance_barplot.pdf',width=7,height=2)

temp = clean_df_meta%>%select(c(1:3,5))%>%pivot_wider(names_from = region,values_from = cluster)
write_csv(temp,'~/output/hooke_result_final/sig_cluster_spline.csv')

########################################################################
#num of subtypes (neurons/non-neurons) in each cluster
clean_df_meta = read.csv('~/output/hooke_result_final/clean_df_meta_update.csv')
all_df = clean_df_meta%>%group_by(group)%>%summarise(n=n())%>%
  mutate(proportion = n / sum(n))%>%
  mutate(cluster='all')
clean_df_meta%>%group_by(cluster,group)%>%summarise(n=n())%>%
  mutate(proportion = n / sum(n),cluster = as.character(cluster))%>%rbind(all_df)%>%
  ggplot(.,aes(cluster,proportion,group=group,fill=group))+geom_bar(stat='identity')+
  scale_fill_manual(values = c('#C4E1F6','#FFBD73'))+
  ylab('num of region-subtypes')+xlab('k-means cluster')

#num of subtypes (region) in each cluster
palette = c('NAc'='#148e34','ACC'='#5555b0','IPP'='#6c9cc0','CN'='#2caa4d','HIP'='#84db55',
            'lCb'='#e6cd69','mdTN'='#117867','dlPFC'='#345296','EC'='#532b70','MB'='#fc7462')
clean_df_meta%>%group_by(cluster,region)%>%summarise(n = n()) %>%
  mutate(proportion = n / sum(n))%>%
  ggplot(.,aes(as.character(cluster),proportion,group=region,fill=region))+geom_bar(stat='identity',color = "black")+theme_classic()+
  ylab('num of region-subtypes')+xlab('k-means cluster')+scale_fill_manual(values=palette)

clean_df_meta%>%group_by(region,cluster)%>%summarise(n = n()) %>%
  mutate(proportion = n / sum(n),cluster=as.character(cluster))%>%
  ggplot(.,aes(region,proportion,group=cluster,fill=cluster))+
  geom_bar(stat='identity',color = "black")+
  scale_fill_brewer(palette='Pastel1')
  ylab('num of region-subtypes')+xlab('k-means cluster')

#################################
#region enrichment
meta = clean_df_meta%>%filter(cell_class=='microglia')
region_mat = meta%>%group_by(region,cluster)%>%summarise(n = n()) %>%
  mutate(proportion = n / sum(n))%>%pivot_wider(id_cols='region',names_from='cluster',values_from='proportion')%>%
  column_to_rownames('region')
all_prop = meta%>%group_by(cluster)%>%summarise(n=n())%>%
  mutate(proportion = n / sum(n))
region_mat = region_mat[,all_prop$cluster]
myBreaks <- c(seq(0,3,length.out=50))
region_mat[is.na(region_mat)]=0
sweep(region_mat[,as.character(seq(1,9,1))], 2, all_prop$proportion, `/`)%>%t()%>%
  pheatmap(cluster_rows = FALSE,cluster_cols = TRUE,
           color=colorRampPalette(c("#219ebc","white","#900603"))(50),
           breaks=myBreaks)
#pheatmap(t(region_mat),cluster_rows = FALSE,cluster_cols = FALSE)
temp = meta%>%select(c(region,cell_subtype,cluster))%>%pivot_wider(names_from = cell_subtype,values_from = cluster)
temp = temp%>%column_to_rownames('region')
pheatmap(temp,cluster_rows = FALSE,cluster_cols = FALSE,display_numbers = TRUE,number_format = "%d")

#cell class enrichment
class_mat = clean_df_meta%>%group_by(cell_class,cluster)%>%summarise(n = n()) %>%
  mutate(proportion = n / sum(n))%>%pivot_wider(id_cols='cell_class',names_from='cluster',values_from='proportion')%>%
  column_to_rownames('cell_class')
all_prop = clean_df_meta%>%group_by(cluster)%>%summarise(n=n())%>%
  mutate(proportion = n / sum(n))
class_mat = class_mat[,all_prop$cluster]
class_mat = class_mat[c('astrocytes','ependymal-cells','microglia','opc-olig','vascular-cells',
                       "GABAergic-neurons","basket-cells","cerebellar-neurons","medium-spiny-neurons",
                       "midbrain-neurons"),]
pheatmap(t(class_mat),cluster_rows = FALSE,cluster_cols = FALSE)


clean_df_meta%>%group_by(cluster,cell_subtype)%>%summarise(n = n()) %>%
  pivot_wider(id_cols='cluster',names_from='cell_subtype',values_from='n')%>%
  column_to_rownames('cluster')%>%mutate_all(~ if_else(is.na(.), 0, .))%>%t()%>%
  pheatmap(cluster_cols = FALSE)
