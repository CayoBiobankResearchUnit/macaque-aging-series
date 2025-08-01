source('/home/wyang114/utils/plotting_theme.R')
source('/home/wyang114/utils/annotation.R')
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


regions = c('ACC','dlPFC','IPP','M1','EC','HIP','CN','NAc','mdTN','MB','lCb')

###nebula outputs (filter it the same way as mashr inputs)
df_nebula = read.csv('/scratch/klwatki5/R_analysis/bbi/nebula/infants_removed/combined_results.csv') #all nebula outputs: 7686734, 1323 subtype x region (including 13 MB-basket and 16 HIP-msn)

####Use 30 donors and 300 cells as cutoff. (filter all basket cells in MB)
cell_count = read.csv('/scratch/klwatki5/R_analysis/bbi/nebula/infants_removed/cell_counts_updated.csv') #1486 (including 13 MB-basket and 16 HIP-msn)
cell_count$identity = paste0(cell_count$cell_type,'_louvain',cell_count$louvain,'_',cell_count$region)
all_identity = cell_count$identity%>%unique()
e.filter.others = cell_count%>%filter(!(cell_type %in% c('basket cells','medium spiny neurons')))%>%
  filter(unique_animal_count<30 | cell_count < 300)%>%pull(identity) #removing 530 region x subtypes
e.filter.basket = cell_count%>%filter(cell_type=='basket cells')%>%
  filter(region=='MB'| unique_animal_count<30 | cell_count < 300)%>%pull(identity) #removing 13 region x subtypes
e.filter.msn = cell_count%>%filter(cell_type=='medium spiny neurons')%>%
  filter(region=='HIP'| unique_animal_count<30 | cell_count < 300)%>%pull(identity)#removing 23 region x subtypes
e.filter = c(e.filter.basket,e.filter.msn,e.filter.others)

###############split into celltype-based datasets
#filtering for NA, n_cells, n_animals, MSN in HIP, basket in MB (844 subtype x region) - 5146072
df_nebula_filter = df_nebula%>%
  filter(!id %in% e.filter)%>%         # filter out region x subtype lower than threshold
  filter(!is.na(p_age))                # filter out genes with NA in p_age or p_sex
write.csv(df_nebula_filter,'/scratch/nsnyderm/u01/intermediate_files/nebula/age/nebula_NBGMM_filtered_rminfants_update.csv') 

###mashr summary
data.dir = '/scratch/klwatki5/R_analysis/bbi/mashr/results/1_raw_matrices'
output.dir = '/scratch/nsnyderm/u01/intermediate_files/nebula/age_olig_split/'
files = list.files(data.dir)
files = files[grep('age_DEG_results',files)]
files = files[files != 'oligodendrocytes_mash_age_DEG_results.RData']
df = NULL
for(file in files){
  load(file.path(data.dir,file))
  ct = str_split_fixed(file,'_',2)[,1]
  print(ct)
  n_col = dim(lfsr)[2]
  df_lfsr = lfsr%>%data.frame()%>%rownames_to_column('gene')%>%pivot_longer(cols = 2:(n_col+1), names_to = 'id', values_to = 'age_mashr_lfsr')
  df_pm = pm%>%data.frame()%>%rownames_to_column('gene')%>%pivot_longer(cols = 2:(n_col+1), names_to = 'id', values_to = 'age_mashr_pm')
  df_psd = psd%>%data.frame()%>%rownames_to_column('gene')%>%pivot_longer(cols = 2:(n_col+1), names_to = 'id', values_to = 'age_mashr_psd')
  df_combined = Reduce(function(x, y) merge(x, y, by = c("gene","id")), list(df_lfsr, df_pm, df_psd))
  #write.csv(df_combined, file.path(output.dir,sprintf('%s_mashr_all_summary.csv',ct)))
  #df_combined = df_combined%>%filter(age_mashr_lfsr<0.05)
  #write.csv(df_combined, file.path(output.dir,sprintf('%s_mashr_filtered_summary.csv',ct)))
  df = rbind(df,df_combined)
}
df$region = str_split_fixed(df$id,'_',3)[,3]
df$louvain =gsub('louvain','',str_split_fixed(df$id,'_',3)[,2])
df$celltype = gsub('\\.',' ',str_split_fixed(df$id,'_',3)[,1])
df$subtype = paste0(df$celltype,'_',df$louvain)
df$id = gsub('\\.',' ',df$id)
df = df%>%mutate(direction = case_when(age_mashr_pm > 0 ~ 'up', age_mashr_pm < 0 ~ 'down'))
write.csv(df,file.path(output.dir,'all_mashr_summary_update.csv')) #all mashr results: 9351513 --> 9289774


###filter to valid tests
df = read.csv(file.path(output.dir,'all_mashr_summary_update.csv'))
df_nebula_filter = read.csv('/scratch/nsnyderm/u01/intermediate_files/nebula/age/nebula_NBGMM_filtered_rminfants_update.csv') #5146072
df_filter = left_join(df, df_nebula_filter,by=c('gene','id'))%>%filter(!is.na(logFC_age)) #5146072 (15117 unique)
write.csv(df_filter,file.path(output.dir,'all_mashr_summary_filter_missing_update.csv')) 
df_filter=read.csv(file.path(output.dir,'all_mashr_summary_filter_missing_update.csv'))
df_sig = df_filter%>%filter(age_mashr_lfsr<0.05) #313672 --> 313576 (9372)
write.csv(df_sig,file.path(output.dir,'all_mashr_summary_lfsr<0.05_filter_missing_update.csv'))

#calculate stats
df = read.csv('/scratch/nsnyderm/u01/intermediate_files/nebula/age_olig_split/all_mashr_summary_filter_missing_update.csv')[-1]
# num of sig vs num of test on cell type
temp_all = df%>%group_by(gene)%>%summarise(n_test = n_distinct(celltype))
temp_sig = df%>%filter(age_mashr_lfsr<0.05)%>%group_by(gene)%>%summarise(n_sig = n_distinct(celltype))
temp = left_join(temp_sig,temp_all,by='gene')
write.csv(temp,'/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/aDEG_test-vs-sig.csv')
temp = read.csv('/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/aDEG_test-vs-sig.csv')[-1]

pdf('/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/aDEG_test-vs-sig.pdf',width = 4,height = 5)
temp %>%count(n_sig, n_test) %>%group_by(n_sig) %>%mutate(perc = n / sum(n)) %>%ggplot(aes(x = n_test, y = perc)) +
  geom_bar(stat = "identity") +
  facet_wrap(~n_sig, scales = "free_y",ncol=3) +xlab('Number of cell classes tested') +
  ggtitle('Number of cell classes detected as aDEGs') +scale_y_continuous(labels = scales::percent_format())#+
  #scale_x_continuous(breaks = 0:11)
dev.off()

# of genes that are cell-class specific
temp%>%filter(n_sig==1)%>%nrow
temp%>%filter(n_sig==1)%>%nrow/nrow(temp)
median(temp$n_sig)
#percentage of cell class specific
unique_gene = temp%>%filter(n_sig==1)%>%pull(gene)
temp = df %>%filter(age_mashr_lfsr<0.05)%>%
  mutate(unique = ifelse(gene %in% unique_gene, 1, 0)) %>%
  group_by(celltype, unique) %>%
  summarise(count = n_distinct(gene), .groups = "drop") %>%
  group_by(celltype) %>%
  mutate(freq = count / sum(count))%>%
  filter(unique==1)
write.csv(temp,'/home/wyang114/output/transcriptional/transcriptional_summary/cell-class_specific_perc.csv')
# number of sig/total
n_sig = df%>%filter(age_mashr_lfsr<0.05)%>%pull(gene)%>%unique()%>%length()
n_total = df%>%pull(gene)%>%unique()%>%length()
n_sig/n_total




###plots
#n_DEG: celltype by regions
df = read.csv('/scratch/nsnyderm/u01/intermediate_files/nebula/age/all_mashr_summary_lfsr<0.05_filter_missing_update.csv')[-1]
ct_color = ct_color%>%filter(identity!='oligodendrocyte precursor cells')
pdf('/home/wyang114/output/transcriptional/transcriptional_summary/olig_split/boxplot_deg_by_celltype.pdf',width = 4,height = 6)
df%>%group_by(celltype,louvain,direction)%>%summarise(n_DEG=n())%>%
  left_join(.,ct_color,by=c('celltype'='identity'))%>%
  ggplot(.,aes(x=reorder(celltype, n_DEG, FUN = mean), y=n_DEG, fill=color, 
               alpha = factor(direction,levels=c('down','up'))))+
  geom_boxplot()+
  scale_fill_identity()+
  scale_y_log10()+               
  xlab('cell type')+coord_flip()
dev.off()

order = c('glutamatergic neurons','oligodendrocytes', 'cerebellar neurons',
          'astrocytes','medium spiny neurons','GABAergic neurons','basket cells',
          'microglia','vascular cells','midbrain neurons','ependymal cells')

#logFC by celltype (add celltype or region specific color to make it prettier)
pdf('/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/boxplot_logFC_by_celltype.pdf',width = 4,height = 6)
ct_color = ct_color%>%filter(identity!='oligodendrocyte precursor cells')
df%>%mutate(abs_age_pm=abs(age_mashr_pm))%>%
  left_join(.,ct_color,by=c('celltype'='identity'))%>%
  ggplot(.,aes(x=factor(celltype,levels=rev(order)), y=abs_age_pm, fill=color, 
               alpha = factor(direction,levels=c('down','up')
               )))+
  geom_boxplot()+scale_fill_identity()+scale_y_log10()+
  xlab('cell type')+ylab('adjusted logFC/year')+coord_flip()
dev.off()
#######supplemental##############
# #Percentage barplot on up and down regulated aDEGs
# df%>%group_by(celltype,direction)%>%summarise(n_DEG=n_distinct(gene))%>%
#   group_by(celltype)%>%mutate(percentage = n_DEG/sum(n_DEG)*100)%>%
#   ggplot(.,aes(x=factor(celltype,levels=rev(order)),y=percentage,fill=direction))+geom_bar(stat='identity',position='stack')+
#   geom_hline(yintercept = 50,linetype = 'dotted')+xlab('celltype')+
#   coord_flip()+scale_fill_manual(values = c('#C4E1F6','#FFBD73'))
# 
# df%>%group_by(region,direction)%>%summarise(n_DEG=n_distinct(gene))%>%
#   group_by(region)%>%mutate(percentage = n_DEG/sum(n_DEG)*100)%>%
#   ggplot(.,aes(x=factor(region,levels=rev(regions)),y=percentage,fill=direction))+geom_bar(stat='identity',position='stack')+
#   geom_hline(yintercept = 50,linetype = 'dotted')+xlab('region')+
#   coord_flip()+scale_fill_manual(values = c('#C4E1F6','#FFBD73'))

#distribution of up aDEG in all 844 subtype x region
#compute percentage of upDEGs in each subcluster x region
df%>%group_by(direction)%>%summarise(n_DEG=n_distinct(gene))

temp = df%>%group_by(id,direction)%>%summarise(n_DEG=n_distinct(gene))%>%
  group_by(id)%>%mutate(percentage = n_DEG/sum(n_DEG)*100)%>%
  ungroup%>%
  complete(id,direction,fill = list(percentage = 0))%>%
  mutate(region=str_split_fixed(id,'_',3)[,3])%>%
  filter(direction=='up')
base = temp$percentage%>%mean()

pdf('/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/hist_perc_up-aDEG.pdf',width = 6,height = 4)
df%>%group_by(id,direction)%>%summarise(n_DEG=n_distinct(gene))%>%
  group_by(id)%>%mutate(percentage = n_DEG/sum(n_DEG)*100)%>%
  ungroup%>%
  complete(id,direction,fill = list(percentage = 0))%>%
  mutate(region=str_split_fixed(id,'_',3)[,3])%>%
  filter(direction=='up')%>%
  ggplot(.,aes(x=percentage))+geom_histogram()+geom_vline(xintercept=base,color='red',linetype='dotted')+
  ylab('num of subtype x regions')+facet_wrap(~region)+ggtitle('percentage of aDEGs')
dev.off()

#more than expected number of subclusterxregion with >57.5% of upDEG
base_perc = sum(temp$percentage>base)/nrow(temp)
up_enrichment = temp%>%group_by(region)%>%summarise(n_total=n(),n_greater=sum(percentage>base),freq=n_greater/n_total)%>%
  rowwise() %>%
  mutate(
    pval = binom.test(n_greater, n_total, p = base_perc)$p.value,
    conf_low = binom.test(n_greater, n_total, p = base_perc)$conf.int[1],
    conf_high = binom.test(n_greater, n_total, p = base_perc)$conf.int[2]
  ) %>%
  ungroup()%>%filter(pval<0.05)
#####################



###check effects of cluster size and seq depth
#check effects of sequencing depth and cluster size
cell_obs = read.csv('/scratch/nsnyderm/u01/intermediate_files/cell-class_h5ad_update/all_regions_obs.csv') ##need to update the cell count after infants
cell_obs = cell_obs%>%filter(age>3)
cell_obs$cell_class_annotation = ifelse(cell_obs$cell_class_annotation=='oligodendrocyte precursor cells','oligodendrocytes',cell_obs$cell_class_annotation)
#test num
df_nebula=read.csv('/scratch/nsnyderm/u01/intermediate_files/nebula/age/nebula_NBGMM_filtered_rminfants_update.csv')[-1]
df_test = df_nebula%>%group_by(id)%>%summarise(test_count=n())
df_test$celltype = str_split_fixed(df_test$id,'_',3)[,1] #id: celltype_louvainX_region
df_test$louvain=gsub('louvain','',str_split_fixed(df_test$id,'_',3)[,2])%>%as.numeric()
df_test$region=str_split_fixed(df_test$id,'_',3)[,3]
df_test%>%ggplot(.,aes(celltype,test_count))+geom_boxplot()

#Check linear relationship between n_DEG and mean UMI and cluster size
df=read.csv('/scratch/nsnyderm/u01/intermediate_files/nebula/age_olig_split/all_mashr_summary_lfsr<0.05_filter_missing_update.csv')[-1]
df_deg =  df%>%group_by(celltype,louvain,region)%>%summarise(n_DEG=n())
df_umi = cell_obs%>%mutate(celltype = cell_class_annotation)%>%group_by(celltype,louvain,region)%>%summarise(mean_umi=mean(n_umi))
df_size = cell_obs%>%mutate(celltype = cell_class_annotation)%>%group_by(celltype,louvain,region)%>%summarise(cluster_size = n())

df_combined <- Reduce(function(x, y) merge(x, y, by = c("celltype","louvain","region")), 
                      list(df_deg, df_umi, df_size,df_test))
write.csv(df_combined,'/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/region_cluster_stat.csv')


#try normalization
df_combined$n_DEG_norm = df_combined$n_DEG/df_combined$test_count
df_combined = df_combined%>%mutate(log_mean_umi = log(mean_umi), log_cluster_size = log(cluster_size),
                                   log_nDEG_norm = log(n_DEG_norm*10000))
df_combined%>%ggplot(.,aes(x=log_mean_umi,y=log_nDEG_norm))+geom_point()+
  geom_smooth(method = "glm", formula = y ~ x, se = TRUE)+geom_point(aes(color=celltype))+
  scale_color_brewer(palette = 'Paired')
df_combined%>%ggplot(.,aes(x=log_cluster_size,y=log_nDEG_norm))+geom_point()+
  geom_smooth(method = "glm", formula = y ~ x, se = TRUE)+geom_point(aes(color=celltype))+
  scale_color_brewer(palette = 'Paired')
df_combined%>%ggplot(.,aes(x=test_count,y=log_nDEG_norm))+geom_point()+
  geom_smooth(method = "glm", formula = y ~ x, se = TRUE)+facet_wrap(~celltype)

# df_combined%>%ggplot(.,aes(x=log_cluster_size,y=log_mean_umi))+geom_point()+
#   geom_smooth(method = "glm", formula = y ~ x, se = TRUE)+geom_point(aes(color=direction))
# df_combined%>%pivot_wider(id_cols = 'id',names_from = 'direction',values_from = 'n_DEG_norm')%>%
#   mutate(across(c('up', 'down'), ~replace(., is.na(.), 0)))%>%
#   ggplot(.,aes(x=up,y=down))+geom_point()+xlim(0,0.4)+ylim(0,0.4)+
#   xlab('upregulated DEG(normalized)')+
#   ylab('downregulated DEG(normalized)')

#Use residuals
library(MASS)
lm_model = lm(data = df_combined, formula = n_DEG_norm ~cluster_size+mean_umi)
summary(lm_model)

lm_model = lm(data = df_combined, formula = log_nDEG_norm ~log_cluster_size+log_mean_umi)
summary(lm_model)
glm_model = glm(data = df_combined, formula = log_nDEG_norm ~log_cluster_size+log_mean_umi,
                family = gaussian(link = "identity")) 
summary(glm_model)
AIC(lm_model, glm_model)
BIC(lm_model, glm_model)


lm_pred <- predict(lm_model, newdata = df_combined)
glm_pred <- predict(glm_model, newdata = df_combined, type = "response")
lm_mse <- mean((df_combined$log_nDEG_norm - lm_pred)^2)
print(lm_mse)
glm_mse <- mean((df_combined$log_nDEG_norm - glm_pred)^2)
print(glm_mse)

weights <- 1 / (fitted(glm_model)^2)  # Example weight function
glm_wls <- glm(log_nDEG_norm ~log_cluster_size+log_mean_umi, data = df_combined, weights = weights)
AIC(lm_model, glm_wls)
BIC(lm_model, glm_wls)
wls_pred <- predict(glm_wls, newdata = df_combined, type = "response")
wls_mse <- mean((df_combined$log_nDEG_norm - wls_pred)^2)
print(wls_mse)

loess_model <- loess(log_nDEG_norm ~log_cluster_size+log_mean_umi, data = df_combined)
loess_pred <- predict(loess_model, newdata = df_combined)
loess_mse <- mean((df_combined$log_nDEG_norm - loess_pred)^2)
print(loess_mse)
loess_rss <- sum((df_combined$log_nDEG_norm - loess_pred)^2)
tss <- sum((df_combined$log_nDEG_norm - mean(df_combined$log_nDEG_norm))^2)
loess_rsq <- 1 - loess_rss / tss
print(loess_rsq)
# plot loess
x1='log_mean_umi'
x2='log_cluster_size'
x2_mean <- mean(df_combined[[x2]])
x1_range <- seq(min(df_combined[[x1]]), max(df_combined[[x1]]), length.out = 100)
predict_data <- data.frame(x1 = x1_range, x2 = rep(x2_mean, 100))
colnames(predict_data)=c(x1,x2)
y_pred <- predict(loess_model, newdata = predict_data)
pdf('/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/loess_fit_DEG_log-enrichment_umi.pdf',width = 6,height = 4)
ggplot(df_combined, aes(x = log_mean_umi, y = log_nDEG_norm, color=celltype)) +
  geom_point() +  # Original data points
  scale_color_brewer(palette = 'Paired')+
  geom_line(data = predict_data, aes(x = log_mean_umi, y = y_pred), color = "red", size = 1) +  # LOESS fit
  labs(title = "LOESS Fit with cluster size and normalized nDEG",
       x = "log(mean umi)",
       y = "log(nDEG/num gene tested * 10E4)") +
  theme_classic()
dev.off()
cat("Residual Sum of Squares (RSS):", loess_rss, "\n")
cat("Mean Squared Error (MSE):", loess_mse, "\n")

x1='log_cluster_size'
x2='log_mean_umi'
x2_mean <- mean(df_combined[[x2]])
x1_range <- seq(min(df_combined[[x1]]), max(df_combined[[x1]]), length.out = 100)
predict_data <- data.frame(x1 = x1_range, x2 = rep(x2_mean, 100))
colnames(predict_data)=c(x1,x2)
y_pred <- predict(loess_model, newdata = predict_data)
pdf('/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/loess_fit_DEG_log-enrichment_clustersize.pdf',width = 6,height = 4)
ggplot(df_combined, aes(x = log_cluster_size, y = log_nDEG_norm, color=celltype)) +
  geom_point() +  # Original data points
  scale_color_brewer(palette = 'Paired')+
  geom_line(data = predict_data, aes(x = log_cluster_size, y = y_pred), color = "red", size = 1) +  # LOESS fit
  labs(title = "LOESS Fit with cluster size and normalized nDEG",
       x = "log(cluster size)",
       y = "log(nDEG/num gene tested * 10E4)") +
  theme_classic()
dev.off()

df_combined$residuals<- residuals(loess_model)
write.csv(df_combined,'/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/region_cluster_stat_loess.csv')


df_combined = read.csv('/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/region_cluster_stat_loess.csv')[-1]
pdf('/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/loess_residual_nDEG.pdf',width = 6,height = 4)
df_combined%>%left_join(.,ct_color,by=c('celltype'='identity'))%>%
  ggplot(.,aes(x=reorder(celltype, -residuals, FUN = median), y=residuals,fill=color))+
  geom_boxplot(outlier.size = 0.5)+
  scale_fill_identity()+
  xlab('celltype')
dev.off()
pdf('/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/loess_residual_nDEG_regionxcelltype.pdf',width = 4,height = 7)
df_combined%>%left_join(.,ct_color,by=c('celltype'='identity'))%>%
  ggplot(.,aes(x=interaction(factor(region,regions),celltype), y=residuals, fill=color))+
  geom_boxplot(outlier.size = 0.5)+
  xlab('region')+scale_fill_identity()+coord_flip()
dev.off()


#heatmap
mat = df_combined%>%dplyr::select(celltype,region,residuals)%>%
  group_by(celltype,region)%>%summarise(mean_residuals=mean(residuals))%>%
  ungroup()%>%
  pivot_wider(names_from = region,values_from = mean_residuals)%>%
  column_to_rownames('celltype')
#mat[is.na(mat)]=0
regions = c('dlPFC','M1','CN','NAc','IPP','lCb','EC','HIP','MB','ACC','mdTN')
order = c('glutamatergic neurons','GABAergic neurons',
           'medium spiny neurons','oligodendrocytes',
          'astrocytes','microglia','vascular cells','ependymal cells')
color_palette <- colorRampPalette(c("blue", "white", "red"))(50)
breaks <- seq(-2,2, length.out = 51)
pdf('/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/loess_residual_heatmap_ordered.pdf',width = 6,height = 4)
pheatmap(mat[order,regions],cluster_rows = F,cluster_cols = T,scale='row',
         breaks = breaks, color=color_palette,clustering_method = 'ward.D2',
         cellwidth = 20,cellheight = 20,treeheight_col = 5)
dev.off()


df_combined%>%arrange(desc(residuals))%>%mutate(id = paste(celltype, paste0('louvain',louvain), region, sep = "_"))%>%
  left_join(.,ct_color,by=c('celltype'='identity'))%>%
  head(20)%>%
  ggplot(.,aes(x=reorder(id, -residuals),y=residuals,fill=color))+geom_bar(stat='identity')+scale_fill_identity()

###############DEGs overlap among regions - spearman's correlation
df_all=read.csv('/scratch/nsnyderm/u01/intermediate_files/nebula/age_olig_split/all_mashr_summary_filter_missing_update.csv')[-1]
cell_obs = read.csv('/scratch/nsnyderm/u01/intermediate_files/cell-class_h5ad_update/all_regions_obs.csv')
cell_count = cell_obs %>% mutate(cell_class_annotation=
                      case_when(cell_class_annotation=='oligodendrocyte precursor cells'~ 'oligodendrocytes',
                                cell_class_annotation!='oligodendrocyte precursor cells'~ cell_class_annotation))%>%
  group_by(region,cell_class_annotation,louvain) %>%
  summarise(count = n(), .groups = "drop", n_animal=n_distinct(animal_id)) %>%
  group_by(cell_class_annotation,louvain) %>%
  mutate(percent = count / sum(count) * 100)%>%
  filter(percent>1 & count>300 & n_animal>30)%>%
  mutate(id=paste(cell_class_annotation,paste0('louvain',louvain),region,sep='_'))
##everything (cell type x region)
mat = df_all%>%group_by(gene,region,celltype)%>%
  filter(id %in% cell_count$id)%>%
  summarise(mean_pm=mean(age_mashr_pm))%>%
  mutate(id = paste0(region,'_',celltype))%>%
  pivot_wider(id_cols=id,names_from = gene,values_from = mean_pm)%>%
  column_to_rownames('id')

cor_matrix <- cor(t(mat), use = "pairwise.complete.obs", method = "spearman")
write.csv(cor_matrix,'/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/celltypexregion_correlation.csv')

#look at nebula only
df_all=read.csv('/scratch/nsnyderm/u01/intermediate_files/nebula/age/nebula_NBGMM_filtered_rminfants_update.csv')[-1]
df_all$celltype = str_split_fixed(df_all$id,'_',2)[,1]
df_all$region = str_split_fixed(df_all$id,'_',3)[,3]
mat = df_all%>%group_by(gene,region,celltype)%>%
  filter(id %in% cell_count$id)%>%
  summarise(mean_lfc=mean(logFC_age))%>%
  mutate(id = paste0(region,'_',celltype))%>%
  pivot_wider(id_cols=id,names_from = gene,values_from = mean_lfc)%>%
  column_to_rownames('id')

cor_matrix <- cor(t(mat), use = "pairwise.complete.obs", method = "spearman")
write.csv(cor_matrix,'/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/celltypexregion_correlation_nebula.csv')

row_annotation <- data.frame(region=str_split_fixed(colnames(cor_matrix),'_',2)[,1],
                             celltype = str_split_fixed(colnames(cor_matrix ),'_',2)[,2])
rownames(row_annotation) <- colnames(cor_matrix)

regions = c('ACC','dlPFC','IPP','M1','EC','HIP','CN','NAc','mdTN','MB','lCb')
df_region = region_color%>%filter(identity %in% regions)
region_colors <- setNames(df_region$color, df_region$identity)
df_ct<- color_map%>%filter(identity %in% unique(df_all$celltype))  
ct_colors <- setNames(df_ct$color, df_ct$identity) 

annotation_colors <- list(
  region = region_colors,
  celltype = ct_colors
)
library(viridis)
pheatmap(cor_matrix,annotation_row = row_annotation,annotation_colors = annotation_colors,
         clustering_method = 'ward.D',
         color=viridis(100))

#Violin/box/beeswarm plots of spearman’s rho with x-axis including: 
#1) each cell type against itself among regions 
#2) glia-neurons 
#3) GABA-glut 
#4) each region among cell types 
cor_matrix=read.csv('/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/celltypexregion_correlation.csv')
cor_matrix = cor_matrix%>%data.frame()%>%column_to_rownames('X')
temp = cor_matrix %>% data.frame()%>%rownames_to_column('sender_id')%>%
  pivot_longer(cols = 2:ncol(.), names_to = 'receiver_id', values_to = 'correlation')
temp$receiver_id = gsub('\\.',' ',temp$receiver_id)
temp$sender_ct = str_split_fixed(temp$sender_id,'_',2)[,2]
temp$sender_rg = str_split_fixed(temp$sender_id,'_',2)[,1]
temp$receiver_ct = str_split_fixed(temp$receiver_id,'_',2)[,2]
temp$receiver_rg = str_split_fixed(temp$receiver_id,'_',2)[,1]
temp = temp%>%filter(sender_id!=receiver_id)
#1) each cell type against itself among regions 
temp1 = temp%>%filter(sender_ct==receiver_ct)
temp1$category = paste0(temp1$sender_ct)
temp1$sender_group=NA
temp1$receiver_group=NA
#2) glia-neurons 
glia = c('astrocytes','microglia','oligodendendrocytes','vascular cells','ependymal cells')
temp$sender_group = ifelse(temp$sender_ct %in% glia, 'glia','neuron')
temp$receiver_group = ifelse(temp$receiver_ct %in% glia, 'glia','neuron')
temp2 = temp%>%filter(sender_group=='neuron',receiver_group=='glia')
temp2$category = 'glia-neurons'
#3) gaba-glut
temp3 = temp%>%filter(sender_ct=='glutamatergic neurons',receiver_ct=='GABAergic neurons')
temp3$category = 'glut-GABA'
#4) each region among cell class
temp4 <- temp %>%filter(sender_rg==receiver_rg)%>%
  mutate(A_min = pmin(sender_id, receiver_id),A_max = pmax(sender_id, receiver_id)) %>%
  distinct(A_min, A_max, .keep_all = TRUE) %>%
  dplyr::select(-A_min, -A_max)
temp4$category = paste0(temp4$sender_rg)
pdf('/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/boxplot_cor_celltypexregion.pdf',width = 6,height = 3)
orders = c('glia-neurons','glut-GABA',
           'MB','HIP','CN','NAc','ACC','mdTN','M1','lCb','dlPFC','EC','IPP',
           'ependymal cells','vascular cells','microglia','GABAergic neurons','astrocytes','glutamatergic neurons',
           'oligodendrocytes','medium spiny neurons')

#cor_summary = rbind(temp1,temp2)%>%rbind(.,temp3)%>%rbind(.,temp4)
cor_summary%>%
  ggplot(.,aes(factor(category, levels=orders),correlation))+
  geom_violin() +
  stat_summary(fun = mean, geom = "point", shape = 16, size = 1, color = "black") +
  #geom_jitter(size=0.1)+
  xlab("") +
  ylim(c(-0.1,1))+
  ylab("Spearman's rho") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5),plot.margin = margin(t = 10, r = 20, b = 20, l = 10))+
  ggtitle('Correlation among mean pm of cell class x region')
dev.off()


write.csv(cor_summary,'/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/cellclassxregion_correlation_summary.csv')
cor_summary = read.csv('/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/cellclassxregion_correlation_summary.csv')[-1]

#glut cortical-subcortical neurons
temp5 = temp%>%filter(sender_ct %in% c('glutamatergic neurons') & 
                        receiver_ct %in% c('glutamatergic neurons'))%>%
  mutate(A_min = pmin(sender_id, receiver_id),A_max = pmax(sender_id, receiver_id)) %>%
  distinct(A_min, A_max, .keep_all = TRUE) %>%
  dplyr::select(-A_min, -A_max)
cortical=c('ACC','dlPFC','M1','IPP','EC')
subcortical = c('HIP','CN','NAc','mdTN')
temp5$sender_group = ifelse(temp5$sender_rg %in% cortical, 'cortical','subcortical')
temp5$receiver_group = ifelse(temp5$receiver_rg %in% cortical, 'cortical','subcortical')
temp5$category = ifelse(temp5$sender_group!= temp5$receiver_group, 'cortical-subcortical',
                        ifelse(temp5$sender_group=='cortical','cortical-cortical',
                               'subcortical-subcortical'))
temp5$sub_cat = paste0(temp5$sender_rg,'-',temp5$receiver_rg)
temp5%>%
  ggplot(.,aes(fct_reorder(category, correlation, .fun = mean),correlation))+
  #geom_boxplot()+
  geom_violin() +
  stat_summary(fun = mean, geom = "point", shape = 16, size = 1, color = "black") +
  #geom_jitter(size=0.1)+
  xlab("") +
  ylim(c(-0.1,1))+
  ylab("Spearman's rho") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5))+
  ggtitle('glutamatergic neurons')

#gaba cortical-subcortical
temp6 = temp%>%filter(sender_ct %in% c('GABAergic neurons') & 
                        receiver_ct %in% c('GABAergic neurons'))%>%
  mutate(A_min = pmin(sender_id, receiver_id),A_max = pmax(sender_id, receiver_id)) %>%
  distinct(A_min, A_max, .keep_all = TRUE) %>%
  dplyr::select(-A_min, -A_max)
cortical=c('ACC','dlPFC','M1','IPP','EC')
subcortical = c('HIP','CN','NAc','mdTN')
temp6$sender_group = ifelse(temp6$sender_rg %in% cortical, 'cortical','subcortical')
temp6$receiver_group = ifelse(temp6$receiver_rg %in% cortical, 'cortical','subcortical')
temp6$category = ifelse(temp6$sender_group!= temp6$receiver_group, 'cortical-subcortical',
                        ifelse(temp6$sender_group=='cortical','cortical-cortical',
                               'subcortical-subcortical'))
temp6$sub_cat = paste0(temp6$sender_rg,'-',temp6$receiver_rg)
temp6%>%
  ggplot(.,aes(fct_reorder(category, correlation, .fun = mean),correlation))+
  #geom_boxplot()+
  geom_violin() +
  stat_summary(fun = mean, geom = "point", shape = 16, size = 1, color = "black") +
  #geom_jitter(size=0.1)+
  xlab("") +
  ylim(c(-0.1,1))+
  ylab("Spearman's rho") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5))+
  ggtitle('GABAergic neurons')

#lcb vs astrocytes
sub_mat = df_all%>%
  filter(celltype %in% c('cerebellar neurons','basket cells','astrocytes'))%>%
  group_by(gene,region,cell_type)%>%
  filter(id %in% cell_count$id)%>%
  mutate(id = paste0(region,'_',cell_type,'_',louvain))%>%
  pivot_wider(id_cols=id,names_from = gene,values_from = age_mashr_pm)%>%
  column_to_rownames('id')
cor_sub_matrix <- cor(t(sub_mat), use = "pairwise.complete.obs", method = "spearman")
write.csv(cor_sub_matrix,'/home/wyang114/output/transcriptional/transcriptional_summary/lcb-astrocytes.csv')
cor_sub_matrix = read.csv('/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/lcb-astrocytes.csv')[-1]
row.names(cor_sub_matrix)=colnames(cor_sub_matrix)
temp = cor_sub_matrix  %>% data.frame()%>%rownames_to_column('sender_id')%>%
  pivot_longer(cols = 2:ncol(.), names_to = 'receiver_id', values_to = 'correlation')
temp$receiver_id = gsub('\\.',' ',temp$receiver_id)
temp$sender_ct = str_split_fixed(temp$sender_id,'_',2)[,2]
temp$sender_rg = str_split_fixed(temp$sender_id,'_',2)[,1]
temp$receiver_ct = str_split_fixed(temp$receiver_id,'_',2)[,2]
temp$receiver_rg = str_split_fixed(temp$receiver_id,'_',2)[,1]
temp = temp%>%filter(sender_id!=receiver_id)%>% mutate(A_min = pmin(sender_id, receiver_id),A_max = pmax(sender_id, receiver_id)) %>%
  distinct(A_min, A_max, .keep_all = TRUE) %>%
  dplyr::select(-A_min, -A_max)
temp$sender_group = ifelse(temp$sender_ct=='astrocytes_9','Bergmann',
                           ifelse(grepl('cerebellar.neurons',temp$sender_ct),'cerebellar',
                                  ifelse(grepl('basket.cells',temp$sender_ct),'basket','astrocytes')))
temp$receiver_group = ifelse(temp$receiver_ct=='astrocytes_9','Bergmann',
                           ifelse(grepl('cerebellar.neurons',temp$receiver_ct),'cerebellar',
                                  ifelse(grepl('basket.cells',temp$receiver_ct),'basket','astrocytes')))
temp$sub_cat = paste0(temp$sender_group,'-',temp$receiver_group)
temp$category = ifelse(temp$sub_cat %in% c('Bergmann-basket','basket-Bergmann'),'Bergmann-basket',
                       ifelse(temp$sub_cat %in% c('Bergmann-cerebellar','cerebellar-Bergmann'),'Bergmann-cerebellar',
                              ifelse(temp$sub_cat %in% c('astrocytes-Bergmann'),'Bergmann-astrocytes',
                                     ifelse(temp$sub_cat %in% c('astrocytes-astrocytes'),'astrocytes-astrocytes',
                                            ifelse(temp$sub_cat %in% c('cerebellar-basket','basket-cerebellar'),'cerebellar-basket',NA)))))
temp = temp%>%filter(!is.na(category))
cor_summary = read.csv('/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/cellclassxregion_correlation_summary.csv')[-1]
regions = c('ACC','dlPFC','IPP','M1','EC','HIP','CN','NAc','mdTN','MB','lCb')

orders = c('MB','HIP','CN','NAc','ACC','mdTN','M1','lCb','dlPFC','EC','IPP',
           'Bergmann-basket','cerebellar-basket','Bergmann-cerebellar','Bergmann-astrocytes','astrocytes-astrocytes')

pdf('/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/boxplot_cor_celltypexregion_lcb.pdf',width = 7,height = 3)
cor_summary%>%filter(category %in% regions)%>%bind_rows(.,temp)%>%
  ggplot(.,aes(factor(category, levels=orders),correlation))+
  #geom_boxplot()+
  geom_violin() +
  stat_summary(fun = mean, geom = "point", shape = 16, size = 1, color = "black") +
  #geom_jitter(size=0.1)+
  xlab("") +
  ylim(c(-0.1,1))+
  ylab("Spearman's rho") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5),plot.margin = margin(t = 10, r = 20, b = 20, l = 10))+
  ggtitle('cerebellum programs')
dev.off()


###########update to here#########
#####specific and shared DEG summary
#upset plot
#library(UpSetR)
#library(ComplexUpset)
library(ComplexHeatmap)
#calculate number of cell type each DEG was found
df= read.csv('/scratch/nsnyderm/u01/intermediate_files/nebula/age_olig_split/all_mashr_summary_lfsr<0.05_filter_missing_update.csv')[-1]
pdf('/home/wyang114/output/transcriptional/transcriptional_summary/olig-split//barplot_nDEG_celltype.pdf',width = 3.1,height = 3.5)
order = c('glutamatergic neurons','medium spiny neurons','astrocytes','cerebellar neurons',
          'oligodendrocytes','basket cells','GABAergic neurons','midbrain neurons','vascular cells',
          'microglia','ependymal cells')
df %>%
  group_by(celltype) %>%
  summarise(n_DEG = n_distinct(gene)) %>%
  mutate(celltype = factor(celltype, levels = rev(order))) %>%
  ggplot(aes(x = celltype, y = n_DEG)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n_DEG), hjust = -0.2, size = 3) +  # move label to the right
  xlab("Cell type") +
  ylab("Number of DEGs") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),limits = c(0,8000)) +
  coord_flip() 
dev.off()

  
#cell type specific DEGs
mat = df%>%group_by(celltype,gene)%>%summarise(mean_pm = mean(age_mashr_pm))%>%
  pivot_wider(names_from=celltype,values_from=mean_pm)%>%column_to_rownames('gene')
mat[is.na(mat)]=0
mat <- mat %>%rownames_to_column('gene')%>%
  rowwise() %>%
  mutate(count_gt_zero = sum(c_across(2:12) != 0)) %>%
  ungroup()
pdf('/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/num_ct_histogram.pdf',width = 8,height = 2)
mat%>%ggplot(.,aes(x=count_gt_zero))+geom_histogram(binwidth = 0.5)+
  stat_bin(binwidth = 0.5, geom = "text", aes(label = ifelse(..count.. > 0, ..count.., "")), vjust = -0.5,size=3)+
  scale_x_continuous(breaks = seq(1, 11, 1))+xlab('num_cell type')+ylab('n_DEG')+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) 
dev.off()

mat%>%filter(count_gt_zero==11)
write.csv(mat,'/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/DEG_by_celltype.csv')



#upset plot
mat = read.csv('/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/DEG_by_celltype.csv')[-1]
upset_mat = mat%>%
  rowwise()%>% mutate(across(2:12, ~ ifelse(. != 0, 1, 0)))%>%
  mutate(across(2:12, ~ as.integer(.)))%>%
  #rowwise()%>% mutate(n_celltype = sum(!is.na(c_across(2:12))))%>%
  filter(count_gt_zero >0)%>%
  ungroup()%>%data.frame()
upset_mat%>%ggplot(.,aes(x=count_gt_zero))+geom_histogram(binwidth = 0.5)+
  stat_bin(binwidth = 0.5, geom = "text", aes(label = ifelse(..count.. > 0, ..count.., "")), vjust = -0.5,size=3)+
  scale_x_continuous(breaks = seq(1, 11, 1))+xlab('num_cell type')+ylab('n_DEG')+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) 
upset_mat[is.na(upset_mat)]=0
upset_mat%>%filter(count_gt_zero==11)

filtered_mat = upset_mat%>%filter(count_gt_zero>=9|count_gt_zero<2)%>%
  dplyr::select(2:12)%>%make_comb_mat()%>%
  t()
pdf('/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/upset_deg_summary.pdf',width = 9,height = 3.5)
UpSet(t(filtered_mat),
      top_annotation = upset_top_annotation(t(filtered_mat), add_numbers = TRUE),
      right_annotation = upset_right_annotation(t(filtered_mat), add_numbers = TRUE,),
) 
dev.off()

df_all= read.csv('/scratch/nsnyderm/u01/intermediate_files/nebula/age/all_mashr_summary_filter_missing.csv')[-1]
share_genes = upset_mat%>%filter(count_gt_zero>=10)%>%dplyr::select(gene,count_gt_zero)
cor_sub = df_all%>%filter(gene %in% share_genes$gene)
#plot some genes: YEATS2,P2RY12,NR3C1,FKBP5
#YEATS2 
pdf('/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/boxplot_YEATS2.pdf',width = 6,height = 4)
df_all%>%filter(external_gene_name=='YEATS2')%>%
  mutate(sig = case_when(age_mashr_lfsr<0.05 ~ 'LFSR<0.05',
                         age_mashr_lfsr>=0.05 ~ 'LFSR>=0.05'))%>%
  ggplot(.,aes(x=region,age_mashr_pm/age_mashr_psd,color = sig))+
  geom_boxplot(outlier.size = 0.5)+facet_wrap(~celltype,ncol=4)+geom_hline(yintercept = 0,linetype = 'dotted')+ggtitle('YEATS2')
dev.off()

unique_genes = upset_mat%>%filter(count_gt_zero==1)%>%dplyr::select(gene,count_gt_zero)
temp = df_all%>%filter(gene %in% unique_genes$gene)
#P2RY12
pdf('/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/boxplot_P2RY12.pdf',width = 6,height = 4)
df_all%>%
  filter(external_gene_name=='P2RY12')%>%
  mutate(sig = case_when(age_mashr_lfsr<0.05 ~ 'LFSR<0.05',
                         age_mashr_lfsr>=0.05 ~ 'LFSR>=0.05'))%>%
  ggplot(.,aes(x=region,age_mashr_pm/age_mashr_psd,color=sig))+
  geom_boxplot(outlier.size = 0.5)+facet_wrap(~celltype,ncol=4)+geom_hline(yintercept = 0,linetype = 'dotted')+ggtitle('P2RY12')
dev.off()
#FKBP5
pdf('/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/boxplot_FKBP5.pdf',width = 6,height = 4)
df_all%>%
  filter(external_gene_name=='FKBP5')%>%
  mutate(sig = case_when(age_mashr_lfsr<0.05 ~ 'LFSR<0.05',
                         age_mashr_lfsr>=0.05 ~ 'LFSR>=0.05'))%>%
  ggplot(.,aes(x=region,age_mashr_pm/age_mashr_psd,color=sig))+
  geom_boxplot(outlier.size = 0.5)+facet_wrap(~celltype,ncol=4)+geom_hline(yintercept = 0,linetype = 'dotted')+ggtitle('FKBP5')
dev.off()
#NR3C1
pdf('/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/boxplot_NR3C1_1d.pdf',width = 11,height = 2)
df_all%>%
  filter(external_gene_name=='NR3C1')%>%
  mutate(sig = case_when(age_mashr_lfsr<0.05 ~ 'LFSR<0.05',
                         age_mashr_lfsr>=0.05 ~ 'LFSR>=0.05'))%>%
  ggplot(.,aes(x=region,age_mashr_pm/age_mashr_psd,color=sig))+
  geom_boxplot(outlier.size = 0.5)+facet_wrap(~celltype,ncol=11)+geom_hline(yintercept = 0,linetype = 'dotted')+ggtitle('NR3C1')
dev.off()







##############################
#subtype specific program - DEGs expressed in 3+ regions for a subtype
#Calculate the % of regions that a DEG is found for a subtype (total number of region is expected regions for a celltype)
#We want to look at region-specific or region-broadly changes in a cell type
#load celltype x num_region dataframe (select only the clusters with a nebula test, 205 subtypes total)
df= read.csv('/scratch/nsnyderm/u01/intermediate_files/nebula/age/all_mashr_summary_lfsr<0.05_filter_missing.csv')[-1]
cell_obs = read.csv('/scratch/nsnyderm/u01/intermediate_files/cell-class_h5ad_update/all_regions_obs.csv') ##need to update the cell count after infants
cell_obs = cell_obs%>%filter(age>3)

result = cell_obs%>%mutate(id = paste(cell_class_annotation,paste0('louvain',louvain),region,sep='_'))%>%
  filter(id %in% unique(df$id))%>%mutate(id=NULL)%>%
  ungroup%>%
  group_by(cell_class_annotation,region)%>%summarise(count=n())%>%
  pivot_wider(names_from = 'region',values_from = 'count')%>%ungroup()%>%
  mutate(across(2:12,function(x){ifelse(x>0,1,0)}))%>%
  mutate(total_region=rowSums(select(.,2:12),na.rm=TRUE))

age_up = df%>%filter(direction=='up')
temp = age_up%>%
  pivot_wider(id_cols = c(gene,celltype,louvain),names_from =region,values_from = age_mashr_pm)%>%
  replace(is.na(.), 0)
up.region = cbind(temp[1:3],ifelse(temp[4:14]>0,1,0))%>%mutate(num_region=rowSums(select(.,4:14)))%>%mutate(direction='up')
up.region = left_join(up.region,result%>%select(cell_class_annotation,total_region),
                      by=c('celltype'='cell_class_annotation'))
up.region$perc_region = up.region$num_region/up.region$total_region

up.region%>%group_by(celltype,perc_region)%>%summarise(count=n())%>%
  ggplot(., aes(x = perc_region, y = count)) +
  geom_bar(stat = "identity", color = "black") +
  facet_wrap(~ celltype) +
  labs(title = "region specificity of upregulated DEGs by celltype", x = "% region with change", y = "num DEGs") +
  theme_minimal() + scale_y_log10()+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

age_down = df%>%filter(direction=='down')
temp = age_down%>%
  pivot_wider(id_cols = c(gene,celltype,louvain),names_from =region,values_from = age_mashr_pm)%>%
  replace(is.na(.), 0)
down.region = cbind(temp[1:3],ifelse(temp[4:14]<0,1,0))%>%mutate(num_region=rowSums(select(.,4:14)))%>%
  mutate(direction='down')
down.region = left_join(down.region,result%>%select(cell_class_annotation,total_region),
                        by=c('celltype'='cell_class_annotation'))
down.region$perc_region = down.region$num_region/down.region$total_region
down.region%>%group_by(celltype,perc_region)%>%summarise(count=n())%>%
  ggplot(., aes(x = perc_region, y = count)) +
  geom_bar(stat = "identity", color = "black") +
  facet_wrap(~ celltype) +
  labs(title = "region specificity of downregulated DEGs by celltype", x = "% region with change", y = "num DEGs") +
  theme_minimal() + scale_y_log10()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#GO terms - do we have enough genes?
#how many share genes do we have for each cell class
up.share=up.region%>%filter(perc_region>0.5)%>%
  filter(!celltype %in% c('basket cells','cerebellar neurons','medium spiny neurons','midbrain neurons'))
down.share=down.region%>%filter(perc_region>0.5)%>%
  filter(!celltype %in% c('basket cells','cerebellar neurons','medium spiny neurons','midbrain neurons'))
rbind(up.share,down.share)%>%  
  group_by(celltype,direction)%>%summarise(count=n_distinct(gene))%>%
  ggplot(.,aes(x=celltype,y=count,fill=direction))+geom_bar(stat='identity',position = position_dodge(width = 0.9))+
  ggtitle('nDEG in >50% region')

up.unique=up.region%>%filter(num_region==1)
down.unique=down.region%>%filter(num_region==1)
rbind(up.unique,down.unique)%>%  
  group_by(celltype,direction)%>%summarise(count=n_distinct(gene))%>%
  ggplot(.,aes(x=celltype,y=count,fill=direction))+geom_bar(stat='identity',position = position_dodge(width = 0.9))+
  ggtitle('nDEG in 1 region')

###########################################
#GO term by cell type
#Do we aggregate by genes? (e.g. aggregate across different subtypes)
df_nebula=read.csv('/scratch/nsnyderm/u01/intermediate_files/nebula/age/nebula_NBGMM_filtered_rminfants.csv')[-1]
df_nebula$celltype = str_split_fixed(df_nebula$id,'_',3)[,1] #id: celltype_louvainX_region
df_nebula$louvain=gsub('louvain','',str_split_fixed(df_nebula$id,'_',3)[,2])%>%as.numeric()
df_nebula$region=str_split_fixed(df_nebula$id,'_',3)[,3]

ct = 'astrocytes'
temp =rbind(up.unique,down.unique)%>%filter(celltype==ct)
gene_group = temp$gene%>%unique
universe = df_nebula%>%
  #filter(celltype==ct)%>%
  pull(gene)%>%unique() #use all genes tested for now
write_lines(universe,'/home/wyang114/output/transcriptional/all_genes_input_mashr.txt',sep = '\n')

gene.df <- AnnotationDbi::select(org.Mmu.eg.db,
                  keys = gene_group,
                  columns = c("ENTREZID", "ENSEMBL"),
                  keytype = "ENSEMBL")
df_gene = left_join(temp,gene.df,by=c('gene'='ENSEMBL'))%>%drop_na()

all.gene.df = AnnotationDbi::select(org.Mmu.eg.db,
       keys = universe,
       columns = c("ENTREZID", "ENSEMBL"),
       keytype = "ENSEMBL")

#KEGG (fun=enrichKEGG,organism=mcc)
df_gene$celltype = factor(df_gene$celltype, levels = c("basket cells","cerebellar neurons",
                                                       "midbrain neurons","medium spiny neurons",
                                                       "glutamatergic neurons","GABAergic neurons",
                                                       "astrocytes","microglia","ependymal cells",
                                                       "oligodendrocytes","vascular cells"
                                                       ))
compKEGG = compareCluster(geneClusters=ENTREZID~direction, data = df_gene, fun = "enrichKEGG",
                          organism='mcc',universe = all.gene.df$ENTREZID,
                          pvalueCutoff=0.05)
dotplot(compKEGG, title = "KEGG analysis",font.size = 8)+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


###########################################
#Visualize cell type and region specificity of marker genes in heatmap
df= read.csv('/scratch/nsnyderm/u01/intermediate_files/nebula/age/all_mashr_summary_lfsr<0.05_filter_missing.csv')[-1]
#cell type specific DEGs
mat = df%>%group_by(celltype,gene)%>%summarise(mean_pm = mean(age_mashr_pm))%>%
  pivot_wider(names_from=celltype,values_from=mean_pm)%>%column_to_rownames('gene')
mat[is.na(mat)]=0
color_palette <- colorRampPalette(c("blue", "white", "red"))(50)
breaks <- seq(-0.03, 0.03, length.out = 51)
#group them into clusters
ph <- pheatmap(mat,clustering_method = "ward.D",show_rownames = FALSE)
hc <- ph$tree_row
row_clusters <- cutree(hc, k = 9)
row_group_colors <- data.frame(Group = factor(row_clusters))
row_group_colors$Group%>%table
rownames(row_group_colors) <- rownames(mat)

pheatmap(mat,cluster_rows = TRUE,cluster_cols = TRUE,show_rownames = FALSE,
         clustering_method = 'ward.D',
         breaks=breaks,color = color_palette,
         annotation_row = row_group_colors,
         treeheight_row = 20,
         treeheight_col = 20)
#GO term for each cluster
row_group_colors = row_group_colors%>%rownames_to_column('gene')
gene_group = row_group_colors$gene%>%unique
universe = read_lines('/home/wyang114/output/transcriptional/all_genes_input_mashr.txt')
gene.df <- AnnotationDbi::select(org.Mmu.eg.db,keys = gene_group,columns = c("ENTREZID", "ENSEMBL"),keytype = "ENSEMBL")
df_gene = left_join(row_group_colors,gene.df,by=c('gene'='ENSEMBL'))%>%drop_na()
all.gene.df = AnnotationDbi::select(org.Mmu.eg.db,keys = universe,columns = c("ENTREZID", "ENSEMBL"),keytype = "ENSEMBL")
compKEGG = compareCluster(geneClusters=ENTREZID~Group, data = df_gene, fun = "enrichKEGG",
                          organism='mcc',universe = all.gene.df$ENTREZID,
                          pvalueCutoff=0.05)
df_KEGG = data.frame(compKEGG)
df_KEGG$ONTOLOGY = 'KEGG'
df_KEGG$Description = gsub(" - Macaca mulatta \\(rhesus monkey\\)", "", df_KEGG$Description)
compGO =compareCluster(geneClusters=ENTREZID~Group, data = df_gene, fun = "enrichGO",OrgDb='org.Mmu.eg.db',
                       ont='ALL',universe = all.gene.df$ENTREZID,pvalueCutoff=0.05)
df_GO = data.frame(compGO)
df_combine = rbind(df_KEGG,df_GO)
df_combine <- df_combine %>%
  mutate(GeneRatio = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))%>%
  mutate(Description=paste0(ONTOLOGY,'-',Description))
###plot
top_terms <- df_combine %>%
  filter(p.adjust<0.05)%>%
  group_by(Group) %>%
  slice_min(order_by = p.adjust, n = 1)%>%  
  rbind(df_combine)

# Create dotplot
ggplot(top_terms, aes(x = Group, y = Description)) +
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


mat = read.csv('/home/wyang114/output/transcriptional/DEG_by_celltype.csv')[-1]
var_info = read.table('/scratch/wyang114/data/var_info.txt',sep=',',header=TRUE,row.names=1)
#cell-type specific genes
temp = mat%>%filter(count_gt_zero==1)%>%column_to_rownames('gene')%>%dplyr::select(count_gt_zero)
temp$external_gene_name = var_info[row.names(temp),'external_gene_name']
temp = temp%>%rownames_to_column('gene')
temp2 = mat%>%filter(count_gt_zero==1)%>%dplyr::select(1:12)%>%
  pivot_longer(cols=2:12,names_to = 'cell_type',values_to = 'mean_pm')%>%
  filter(mean_pm!=0) #2513
all_unique = left_join(temp,temp2,by='gene')
all_unique%>%group_by(cell_type)%>%summarise(n_unique=n)%>%
  ggplot(.,aes(x=cell_type,n_unique))+geom_bar(stat='identity')

write.csv(all_unique,'/home/wyang114/output/transcriptional/celltype_specific_DEG.csv')

#celltype specific GO
all_unique = all_unique%>%rownames_to_column('gene')
gene_group = all_unique$gene%>%unique
universe = read_lines('/home/wyang114/output/transcriptional/all_genes_input_mashr.txt')
gene.df <- AnnotationDbi::select(org.Mmu.eg.db,keys = gene_group,columns = c("ENTREZID", "ENSEMBL"),keytype = "ENSEMBL")
df_gene = left_join(all_unique,gene.df,by=c('gene'='ENSEMBL'))%>%drop_na()
all.gene.df = AnnotationDbi::select(org.Mmu.eg.db,keys = universe,columns = c("ENTREZID", "ENSEMBL"),keytype = "ENSEMBL")
compKEGG = compareCluster(geneClusters=ENTREZID~cell_type, data = df_gene, fun = "enrichKEGG",
                          organism='mcc',universe = all.gene.df$ENTREZID,
                          pvalueCutoff=0.05)
dotplot(compKEGG, title = "celltype-specific DEG",font.size = 6)+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
data = as.data.frame(compKEGG)

#celltype share DEG
temp = mat%>%filter(count_gt_zero>=9)%>%column_to_rownames('gene')
temp$external_gene_name = var_info[row.names(temp),'external_gene_name']
temp = temp%>%rownames_to_column('gene')
all_share = temp%>%dplyr::select(gene,external_gene_name,count_gt_zero)
color_palette <- colorRampPalette(c("blue", "white", "red"))(50)
breaks <- seq(-0.03, 0.03, length.out = 51)
#group them into clusters
pheatmap(temp[2:12],clustering_method = "ward.D2",show_rownames = FALSE,breaks=breaks,color = color_palette,
         treeheight_row = 20,treeheight_col = 20,)

#celltype share GO
gene_group = temp$gene%>%unique
universe = read_lines('/home/wyang114/output/transcriptional/all_genes_input_mashr.txt')
gene.df <- AnnotationDbi::select(org.Mmu.eg.db,keys = gene_group,columns = c("ENTREZID", "ENSEMBL"),keytype = "ENSEMBL")
df_gene = left_join(temp,gene.df,by=c('gene'='ENSEMBL'))%>%drop_na()
all.gene.df = AnnotationDbi::select(org.Mmu.eg.db,keys = universe,columns = c("ENTREZID", "ENSEMBL"),keytype = "ENSEMBL")
kegg_result <- enrichKEGG(df_gene$ENTREZID,organism = 'mcc',universe = all.gene.df$ENTREZID,
                          pvalueCutoff = 0.05)
dotplot(kegg_result, title = "celltype share DEG",font.size = 6)+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
data = as.data.frame(kegg_result)

write.csv(all_share,'/home/wyang114/output/transcriptional/celltype_share_DEG.csv')


###########################################
#region specific DEGs
mat = df%>%group_by(region,gene)%>%summarise(mean_pm = mean(age_mashr_pm))%>%
  pivot_wider(names_from=region,values_from=mean_pm)%>%column_to_rownames('gene')
mat[is.na(mat)]=0
color_palette <- colorRampPalette(c("blue", "white", "red"))(50)
breaks <- seq(-0.03, 0.03, length.out = 51)
ph <- pheatmap(mat,clustering_method = "ward.D",show_rownames = FALSE)
hc <- ph$tree_row
row_clusters <- cutree(hc, k = 9)
row_group_colors <- data.frame(Group = factor(row_clusters))
row_group_colors$Group%>%table
rownames(row_group_colors) <- rownames(mat)

pheatmap(mat,cluster_rows = TRUE,cluster_cols = TRUE,show_rownames = FALSE,
         clustering_method = 'ward.D',
         breaks=breaks,color = color_palette,
         annotation_row = row_group_colors,
         treeheight_row = 20,
         treeheight_col = 20)

#GO term for each cluster
row_group_colors = row_group_colors%>%rownames_to_column('gene')
gene_group = row_group_colors$gene%>%unique
universe = read_lines('/home/wyang114/output/transcriptional/all_genes_input_mashr.txt')
gene.df <- AnnotationDbi::select(org.Mmu.eg.db,keys = gene_group,columns = c("ENTREZID", "ENSEMBL"),keytype = "ENSEMBL")
df_gene = left_join(row_group_colors,gene.df,by=c('gene'='ENSEMBL'))%>%drop_na()
all.gene.df = AnnotationDbi::select(org.Mmu.eg.db,keys = universe,columns = c("ENTREZID", "ENSEMBL"),keytype = "ENSEMBL")
compKEGG = compareCluster(geneClusters=ENTREZID~Group, data = df_gene, fun = "enrichKEGG",
                          organism='mcc',universe = all.gene.df$ENTREZID,
                          pvalueCutoff=0.05)
dotplot(compKEGG, title = "KEGG analysis",font.size = 8)+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))# 650x1000

mat <- mat %>%rownames_to_column('gene')%>%
  rowwise() %>%
  mutate(count_gt_zero = sum(c_across(2:12) != 0)) %>%
  ungroup()
mat%>%ggplot(.,aes(x=count_gt_zero))+geom_histogram(binwidth = 0.5)+
  stat_bin(binwidth = 0.5, geom = "text", aes(label = ifelse(..count.. > 0, ..count.., "")), vjust = -0.5)+
  scale_x_continuous(breaks = seq(1, 11, 1))+xlab('num_regions')+ylab('n_DEG')
write.csv(mat,'/home/wyang114/output/transcriptional/DEG_by_region.csv')

#find Bergmann specific programs
df= read.csv('/scratch/nsnyderm/u01/intermediate_files/nebula/age/all_mashr_summary_lfsr<0.05_filter_missing.csv')[-1]
astrocytes = df%>%filter(celltype=='astrocytes')
unique_genes = astrocytes%>%group_by(gene)%>%summarise(n=n_distinct(region))%>%filter(n==1)%>%pull(gene)
lCb_gene = astrocytes%>%filter(subtype=='astrocytes_9',region=='lCb')%>%pull(gene)%>%unique()
unique_astr_genes = intersect(unique_genes,lCb_gene)
unique_cereb_genes = df%>%filter(celltype=='cerebellar neurons')%>%group_by(gene)%>%summarise(n=n_distinct(region))%>%filter(n==1)%>%pull(gene)
unique_basket_genes = df%>%filter(celltype=='basket cells')%>%group_by(gene)%>%summarise(n=n_distinct(region))%>%filter(n==1)%>%pull(gene)
unique_lCb_genes = intersect(intersect(unique_astr_genes,unique_cereb_genes),unique_basket_genes) #92 genes


df= read.csv('/scratch/nsnyderm/u01/intermediate_files/nebula/age/all_mashr_summary_filter_missing.csv')[-1]
lCb_subset = df%>%filter(gene %in% unique_lCb_genes)
write.csv(lCb_subset,'/home/wyang114/output/transcriptional/lCb_specific_DEG.csv')
lCb_subset = read.csv('/home/wyang114/output/transcriptional/lCb_specific_DEG.csv')[-1]
lCb_subset = lCb_subset%>%filter(celltype %in% c('astrocytes','cerebellar neurons','basket cells'))
mat = lCb_subset%>%dplyr::select(id,region,gene,age_mashr_pm)%>%
  pivot_wider(id_cols=id,names_from = gene,values_from = age_mashr_pm)%>%column_to_rownames('id')

annotation_df = data.frame(id=row.names(mat))
annotation_df$region=str_split_fixed(annotation_df$id,'_',3)[,3]
annotation_df$cell_class=str_split_fixed(annotation_df$id,'_',3)[,1]
annotation_df$louvain=str_split_fixed(annotation_df$id,'_',3)[,2]
annotation_df[annotation_df$id=='astrocytes_louvain9_lCb','cell_class'] ='Bergmann' 
annotation_df = annotation_df%>%column_to_rownames('id')
annotation_df$louvain=NULL

df_region = region_color%>%filter(identity %in% annotation_df$region)
region_colors <- setNames(df_region$color, df_region$identity)
df_ct<- color_map%>%filter(identity %in% unique(annotation_df$cell_class))%>%
  rbind(data.frame(color='#FFFF00',identity='Bergmann'))  
ct_colors <- setNames(df_ct$color, df_ct$identity) 

annotation_colors <- list(
  region = region_colors,
  cell_class = ct_colors
)
pheatmap(mat,
         show_colnames = FALSE,
         show_rownames = FALSE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_colors = annotation_colors,
         na_col = "grey",
         clustering_method = 'ward.D',
         annotation_row = annotation_df,
         breaks = seq(-0.02, 0.02, length.out = 100),
         fontsize = 8,
         border_color = NA,
         treeheight_row = 10,
         treeheight_col = 10)

