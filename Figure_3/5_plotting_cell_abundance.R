#' @file 5_plotting_cell_abundance.R
#' @author Wei Yang
#'
#' @description
#' Generates panels for Figure 3
#'
#' @usage
#' # Open in RStudio and execute lines sequentially (Cmd+Enter / Ctrl+Enter)
#' # Or render interactively:
#' source("5_plotting_cell_abundance.R", echo = TRUE)
NULL

source('/home/wyang114/utils/plotting_theme.R')
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyverse)
library(pheatmap)

class_level = c('astrocytes','oligodendrocytes','oligodendrocyte precursor cells',
                'microglia','ependymal cells','vascular cells',
                'glutamatergic neurons','GABAergic neurons','medium spiny neurons',
                'cerebellar neurons','basket cells','midbrain neurons')
region_level=c('ACC','dlPFC','IPP','M1','EC','HIP','CN','NAc','mdTN','MB','lCb')
data.dir = '/home/wyang114/output/hooke_result_final/'

############cell class############
#dotplot cell class log-fold change across full life/aged animals 
cell_class_abundance <- readr::read_csv(paste0(data.dir,"cell-class_full_linear_ctrl-sex_age>3.csv"))[-1]
cell_class_abundance$significance = ifelse(cell_class_abundance$age_padjust<0.05,1,0)
cell_class_abundance$significance = ifelse(cell_class_abundance$age_coef<0,
                                           -cell_class_abundance$significance,
                                           cell_class_abundance$significance)%>%as.factor()
cell_class_abundance$celltype=factor(cell_class_abundance$celltype,
                                     levels=class_level)
cell_class_abundance = cell_class_abundance %>%
  filter(!(celltype=='medium spiny neurons' & region=='HIP'))


cell_class_abundance%>%ggplot(.)+
  geom_point(aes(x=age_coef,y=celltype,color='NS'),size = 3)+
  geom_point(data = cell_class_abundance%>%filter(significance!=0),
             aes(x=age_coef,y=celltype,color=significance,),size = 3)+
  coord_cartesian(c(-0.15,0.15))+
  xlab('log-fold change/year')+
  scale_color_manual(values = c('#C4E1F6','#FFBD73','grey'))+
  geom_vline(xintercept=-1,linetype=3,colour = 'grey')+
  geom_vline(xintercept=1,linetype=3,colour = 'grey')+
  theme_classic()
ggsave(filename = paste0('/home/wyang114/output/hooke_result_final/cell-class_dotplot.pdf'),width = 6,height=3)

#heatmap (600 x 300)
cell_class_abundance$heatmap_coef = cell_class_abundance$age_coef
#cell_class_abundance$heatmap_coef = ifelse(cell_class_abundance$significance==0, 0, cell_class_abundance$heatmap_coef)

myBreaks <- c(seq(-0.05, 0, length.out=ceiling(50/2) + 1), 
              seq(0.05/50, 0.05, length.out=floor(50/2)))
mat = cell_class_abundance%>%pivot_wider(id_cols=celltype,names_from = region,values_from = heatmap_coef)%>%
  column_to_rownames('celltype')
p_value_matrix = cell_class_abundance%>%pivot_wider(id_cols=celltype,names_from = region,values_from = age_padjust)%>%
  column_to_rownames('celltype')
significance_matrix <-  ifelse(is.na(p_value_matrix), "",
                               ifelse(p_value_matrix < 0.001, "***",
                                      ifelse(p_value_matrix < 0.01, "**",
                                             ifelse(p_value_matrix < 0.05, "*", ""))))
pheatmap(mat[,region_level],color=colorRampPalette(c("#219ebc","white","#900603"))(50),
         cluster_rows = F, cluster_cols = F,breaks=myBreaks,
         display_numbers = significance_matrix[,region_level],fontsize_number = 13,
         filename = paste0('/home/wyang114/output/hooke_result_final/cell-class_heatmap.pdf'),width = 5.5,height=3)

#Plot log cell abundance 
cell_class_count = readr::read_csv(paste0(data.dir,"cell-class_full_linear_ctrl-sex_counts.csv"))[-1]
subset = cell_class_count%>%filter(cell_group %in% c('oligodendrocytes','oligodendrocyte precursor cells'))
subset$region = factor(subset$region,levels=region_level)
subset$age_group = cut(subset$age,breaks = c(0,3,9,15,21))
subset%>%ggplot(.,aes(region,num_cells,groups=age_group,color=age_group))+geom_boxplot(outlier.size = 1)+
  facet_wrap(~cell_group,nrow=2,scales = 'free_y')+ylab('log_abundance') #400x400
ggsave(filename = paste0('/home/wyang114/output/hooke_result_final/log_abundance_olig_opc.pdf'),width = 3,height=3)

#log abundance correlation
df = cell_class_count%>%mutate(id = paste(region,cell_group,sep='_'))%>%
  select(c(animal,id,num_cells))%>%
  pivot_wider(names_from = animal,values_from = num_cells)%>%column_to_rownames('id')
df[is.na(df)]=0
df = t(df)%>%scale()%>%t()
correlation_matrix <- cor(df, method = "spearman")
temp = cell_class_count%>%group_by(animal,age)%>%summarise(n=n())%>%arrange(age)
mat=correlation_matrix[temp$animal,temp$animal]
row.names(mat) = temp$age
colnames(mat) = temp$age
breaks <- seq(-0.4, 0.4, length.out = 100)
pheatmap(mat,cluster_rows = F, cluster_cols = F,fontsize = 5,breaks=breaks)

############cell subtype############
#dotplot cell subtype log-fold change across aging with cell-class x region
class_level = c('astrocytes','opc-olig',
                'microglia','ependymal cells','vascular cells',
                'glutamatergic neurons','GABAergic neurons','medium spiny neurons',
                'cerebellar neurons','basket cells','midbrain neurons')
fdr_threshold=0.05
cell_subtype_abundance <- readr::read_csv("/home/wyang114/output/hooke_result_final/cell-subtype_full_linear_ctrl-sex_mash_age>3.csv")
cell_subtype_abundance = cell_subtype_abundance%>%
                        filter(!(celltype=='medium-spiny-neurons' & region=='HIP'))
cell_subtype_abundance = cell_subtype_abundance%>%
  filter(count>=100) #1039

sum(cell_subtype_abundance$age_mashr_lfsr<0.05)#101


cell_subtype_abundance$celltype = ifelse(cell_subtype_abundance$celltype!='opc-olig',
                                    gsub('-',' ',cell_subtype_abundance$celltype),
                                    cell_subtype_abundance$celltype)
cell_subtype_abundance$cell_subtype = paste0(cell_subtype_abundance$celltype,'_',cell_subtype_abundance$louvain)
non_neurons = c('astrocytes','opc-olig','ependymal-cells','vascular-cells','microglia')
cell_subtype_abundance$group = ifelse(cell_subtype_abundance$celltype %in% non_neurons, 'non_neurons','neurons')
cell_subtype_abundance$significance = as.numeric(cell_subtype_abundance$age_mashr_lfsr<fdr_threshold)
cell_subtype_abundance$significance = ifelse(cell_subtype_abundance$age_mashr_pm<0,
                                             -cell_subtype_abundance$significance,
                                             cell_subtype_abundance$significance)
cell_subtype_abundance$significance=as.factor(cell_subtype_abundance$significance)
cell_subtype_abundance$region=factor(cell_subtype_abundance$region,
                                     levels=region_level)

cell_subtype_abundance$celltype=factor(cell_subtype_abundance$celltype,
                                       levels=class_level)

cell_subtype_abundance%>%
  ggplot(.)+
  geom_point(aes(x=age_mashr_pm,y=interaction(region,celltype),color='NS'),size = 1.5)+
  geom_point(data = cell_subtype_abundance%>%filter(significance!=0),
             aes(x=age_mashr_pm,y=interaction(region,celltype),color=significance),size = 1.5)+
  coord_cartesian(c(-0.4,0.4))+
  xlab('log-fold change/year')+
  ylab('')+
  scale_color_manual(values = c('darkgreen','purple','grey'))+
  geom_vline(xintercept=-1,linetype=3,colour = 'grey')+
  geom_vline(xintercept=1,linetype=3,colour = 'grey')+
  theme_classic()
ggsave(filename = paste0('/home/wyang114/output/hooke_result_final/cell-subtype_dotplot_gt100.pdf'),width = 6.1,height=10)


#Summarise number of subtypes in each region
cell_subtype_abundance%>%filter(significance!=0)%>%
  group_by(region,group)%>%
  summarise(n=n())%>%
  ggplot(.,aes(x=factor(region,levels=region_level),y=n,fill=group))+geom_bar(stat='identity')+
  scale_fill_manual(values=c('#C4E1F6','#FFBD73'))+
  ylim(0,13)+
  xlab('regions')+
  ylab('number of subtypes')
ggsave(filename = paste0('/home/wyang114/output/hooke_result_final/cell-subtype_region_barplot_gt100.pdf'),width = 4,height=2)

#Summarise number of regions showing significant change for each subtype
cell_subtype_abundance$cell_subtype = paste0(cell_subtype_abundance$celltype,'_',cell_subtype_abundance$louvain)
temp = cell_subtype_abundance%>%filter(significance!=0)%>%group_by(group,cell_subtype)%>%summarise(n=n())
temp%>%ggplot(.,aes(x=n,fill=group))+geom_bar(position = position_dodge2(preserve = "single"))+
  ylim(0,12)+
  scale_fill_manual(values=c('#C4E1F6','#FFBD73'))+
  scale_x_continuous(n.breaks=11)+
  ylim(0,13)+
  xlab('number of regions')+
  ylab('number of subtypes')
ggsave(filename = paste0('/home/wyang114/output/hooke_result_final/cell-subtype_num_barplot_gt100.pdf'),width = 4,height=2)

#heatmap
mat = cell_subtype_abundance%>%
  pivot_wider(id_cols=cell_subtype,names_from = region,values_from = age_mashr_pm)%>%
  column_to_rownames('cell_subtype')
res.lfsr = cell_subtype_abundance%>%
  pivot_wider(id_cols=cell_subtype,names_from = region,values_from = age_mashr_lfsr)%>%
  column_to_rownames('cell_subtype')
mat[res.lfsr>0.05]=0
# #Set missing values to NA
# coef_mat = df%>%mutate(cellsubtype = paste0(celltype,'_',louvain))%>%
#   pivot_wider(id_cols = cell_subtype, names_from = region, values_from = age_coef)%>%
#   column_to_rownames('cell_subtype')
#mat[is.na(coef_mat)]=NA
mat = mat[,c('CN','NAc')]
mat = mat[rowSums(mat!=0,na.rm = T)>0,]
#mat = mat[select.row,]
myBreaks <- c(seq(-0.1, 0, length.out=ceiling(50/2) + 1), 
              seq(0.1/50, 0.1, length.out=floor(50/2)))
pheatmap(mat[sort(row.names(mat)),region_level],show_rownames = TRUE,cluster_cols = FALSE,cluster_rows=FALSE,
         color=colorRampPalette(c("#219ebc","white","#900603"))(50),
         border_color=NA,breaks=myBreaks,fontsize=18,
         filename = '~/output/hooke_result_final/all_subtype_mashr_data_all_age>3_gt100.pdf',width = 12,height = 20)
pheatmap(mat,show_rownames = TRUE,cluster_cols = FALSE,cluster_rows=FALSE,
         color=colorRampPalette(c("#219ebc","white","#900603"))(50),
         border_color=NA,breaks=myBreaks,fontsize=5)



#log abundance correlation
cell_subtype_count = readr::read_csv(paste0(data.dir,"cell-subtype_full_linear_ctrl-sex_counts.csv"))[-1]
cell_subtype_count = cell_subtype_count%>%mutate(subtype = paste(celltype,cell_group,sep='_'))
#Plot log cell abundance 
subset = cell_subtype_count%>%filter(subtype %in% c('microglia_4','microglia_5'))
subset$region = factor(subset$region,levels=region_level)
subset$age_group = cut(subset$age,breaks = c(0,3,9,15,21))
subset%>%ggplot(.,aes(region,num_cells,groups=age_group,color=age_group))+geom_boxplot(outlier.size = 1)+
  facet_wrap(~subtype,nrow=1,scales = 'free_y')+ylab('log_abundance') #400x400
ggsave(filename = paste0('/home/wyang114/output/hooke_result_final/log_abundance_microglia_subtype.pdf'),width = 4.5,height=2.5)



# ###stacked barplot of region/celltype/detection method distribution
# #add information of non-linear pattern
# spline_df = read.csv('/home/wyang114/output/hooke_result_final/clean_df_meta_update.csv')
# df_all = left_join(cell_subtype_abundance,spline_df,by=c('region','cell_subtype'))
# cluster_class = data.frame('cluster'=seq(1,9,1),
#                            'spline_cat'=c('dev depletion','expansion', 'depletion','age expansion',
#                                           'age depletion','expansion','dev expansion','depletion',
#                                           'depletion'))
# ##sig regions
# all = df_all%>%group_by(region)%>%
#   summarise(n=n())%>%mutate(freq=n/sum(n))%>%mutate(spline_cat='z-all')
# temp = left_join(df_all,cluster_class, by='cluster')%>%
#   filter(significance!=0)%>%
#   group_by(spline_cat,region)%>%summarise(n=n())%>%mutate(freq=n/sum(n))
# temp %>%
#   ggplot(aes(x=spline_cat,y=n)) +
#   geom_col() + 
#   coord_flip() + 
#   labs(x="Aging trajectory", y="n")+
#   ylab('number of significant subtypes')
# rbind(temp,all) %>%
#   ggplot(aes(x=spline_cat,y=freq, fill=region)) +
#   geom_col() + 
#   coord_flip() + 
#   labs(x="Aging trajectory", y="perc of subclusters")+
#   scale_fill_brewer(palette = 'Paired')
#   
# ##sig celltypes
# all = df_all%>%group_by(cell_class)%>%
#   summarise(n=n())%>%mutate(freq=n/sum(n))%>%mutate(spline_cat='z-all')
# temp = left_join(df_all,cluster_class, by='cluster')%>%
#   filter(significance!=0)%>%
#   group_by(spline_cat,cell_class)%>%summarise(n=n())%>%mutate(freq=n/sum(n))
# temp %>%
#   ggplot(aes(x=spline_cat,y=n)) +
#   geom_col() + 
#   coord_flip() + 
#   labs(x="Aging trajectory", y="n")
# rbind(all,temp) %>%
#   ggplot(aes(x=spline_cat,y=freq, fill=cell_class)) +
#   geom_col() + 
#   coord_flip() + 
#   labs(x="Aging trajectory", y="perc of subclusters")+
#   scale_fill_brewer(palette = 'Spectral')
# 
# sig_class = data.frame('sig_class'=c('adult-only','full_only','both'),'adult_significance'=c(1,0,1),'significance'=c(0,1,1))
# cell_subtype_abundance = left_join(cell_subtype_abundance,sig_class,by=c('adult_significance','significance'))
# temp = cell_subtype_abundance%>%filter(combined_significance!=0)%>%
#   group_by(spline_cat,sig_class)%>%summarise(n=n())%>%mutate(freq=n/sum(n))
# temp %>%
#   ggplot(aes(x=spline_cat,y=freq, fill=sig_class)) +
#   geom_col() + 
#   coord_flip() + 
#   labs(x="Aging trajectory", y="perc of subclusters")+
#   scale_fill_brewer(palette = 'Set2')+theme_classic()