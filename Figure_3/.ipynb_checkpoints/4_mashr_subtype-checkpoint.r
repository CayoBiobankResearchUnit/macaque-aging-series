source('/home/wyang114/utils/plotting_theme.R')
library(dplyr)
library(tidyverse)
library(mashr)
library(ashr)
library(pheatmap)
set.seed(2024)

df = read.csv('/home/wyang114/output/hooke_result_final/cell-subtype_full_linear_ctrl-sex_age>3.csv',row.names = 1)
df = df%>%filter(!(region=='MB'&celltype=='basket-cells')) #filter out basket cells in MB because it's from dissection error

#For missing data, set coefficient to 0 and stderr to 1000
df$cell_subtype = paste0(df$celltype,'_',df$louvain)
coef_mat = df%>%mutate(cellsubtype = paste0(celltype,'_',louvain))%>%
  pivot_wider(id_cols = cell_subtype, names_from = region, values_from = age_coef)%>%
  column_to_rownames('cell_subtype') #226 x 11
coef_mat[is.na(coef_mat)] = 0
se_mat = df%>%mutate(cellsubtype = paste0(celltype,'_',louvain))%>%
  pivot_wider(id_cols = cell_subtype, names_from = region, values_from = age_se)%>%
  column_to_rownames('cell_subtype')
se_mat[is.na(se_mat)] = 1000

#Since we only have 226 x 11 conditions, we will use the data-driven vignettes
Bhat = coef_mat%>%as.matrix()
Shat = se_mat%>%as.matrix()
data   = mash_set_data(Bhat, Shat)
#select strong signals
m.1by1 = mash_1by1(data)
strong = get_significant_results(m.1by1,0.05) #29 rows from 53 results with lfsr<0.05
#compute covariance matrix by pca from strong signals (covariance small in general)
U.pca = cov_pca(data,5,subset=strong)
print(names(U.pca))
cov2cor(U.pca$PCA_1) #correlation of PCA_1 all as 1,potentially too few variance and datapoints
cov2cor(U.pca$tPCA)#correlation of tPCA looks ok 
dist = as.dist(1 - cov2cor(U.pca$tPCA))
hc <- hclust(dist, method = "ward.D2")
plot(hc) #check hierarchical structures
#Apply Extreme Deconvolution
U.ed = cov_ed(data, U.pca, subset=strong)
#fit mash: we must fit mash to all tests
m.ed = mash(data, U.ed)
print(get_loglik(m.ed),digits = 10) #-5911

#Running with both cannonical and data-driven
U.c = cov_canonical(data)
m   = mash(data, c(U.c,U.ed))
print(get_loglik(m),digits = 10) #-5824 canonical

#Use data-driven approach
res.lfsr = get_lfsr(m.ed)
res.pm = get_pm(m.ed)
res.psd = get_psd(m.ed)

df_lfsr = res.lfsr%>%data.frame()%>%rownames_to_column('cell_subtype')%>%
  pivot_longer(cols=2:12,names_to = 'region', values_to = 'age_mashr_lfsr')
df_pm = res.pm%>%data.frame()%>%rownames_to_column('cell_subtype')%>%
  pivot_longer(cols=2:12,names_to = 'region', values_to = 'age_mashr_pm')
df_psd = res.psd%>%data.frame()%>%rownames_to_column('cell_subtype')%>%
  pivot_longer(cols=2:12,names_to = 'region', values_to = 'age_mashr_psd')

df_all <- Reduce(function(x, y) merge(x, y, by = c("region","cell_subtype")), list(df, df_lfsr, df_pm, df_psd))
write.csv(df_all,'/home/wyang114/output/hooke_result_final/cell-subtype_full_linear_ctrl-sex_mash_age>3.csv')

df_all = read.csv('/home/wyang114/output/hooke_result_final/cell-subtype_full_linear_ctrl-sex_mash_age>3.csv')[-1]
region_level=c('dlPFC','M1','IPP','ACC','CN','NAc','EC','HIP','mdTN','MB','lCb')
#res.pm[res.lfsr>0.05]=0
mat = res.pm
#Set missing values to NA
coef_mat = df%>%mutate(cellsubtype = paste0(celltype,'_',louvain))%>%
  pivot_wider(id_cols = cell_subtype, names_from = region, values_from = age_coef)%>%
  column_to_rownames('cell_subtype')
mat[is.na(coef_mat)]=NA
mat = mat[rowSums(mat!=0,na.rm = T)>0,]
mat = mat[select.row,]
myBreaks <- c(seq(-0.1, 0, length.out=ceiling(50/2) + 1), 
              seq(0.1/50, 0.1, length.out=floor(50/2)))
pheatmap(mat[sort(row.names(mat)),region_level],show_rownames = TRUE,cluster_cols = FALSE,cluster_rows=FALSE,
         color=colorRampPalette(c("#219ebc","white","#900603"))(50),
         border_color=NA,breaks=myBreaks,fontsize=18,
         filename = '~/output/hooke_result_final/all_subtype_mashr_data_all_age>3.pdf',width = 12,height = 20)

#Look at prior coefficients
mat = coef_mat
padjust_mat = df%>%mutate(cellsubtype = paste0(celltype,'_',louvain))%>%
  pivot_wider(id_cols = cellsubtype, names_from = region, values_from = age_padjust)%>%
  column_to_rownames('cellsubtype')
mat[padjust_mat>0.05]=0
mat = mat[select.row,]
myBreaks <- c(seq(-0.1, 0, length.out=ceiling(50/2) + 1), 
              seq(0.1/50, 0.1, length.out=floor(50/2)))
library(pheatmap)
pheatmap(mat[,region_level],show_rownames = TRUE,cluster_cols = FALSE,cluster_rows=FALSE,
         color=colorRampPalette(c("#219ebc","white","#900603"))(50),
         border_color=NA,breaks=myBreaks,
         filename = '~/output/hooke_result_final/all_subtype_data_prior.pdf',width = 10,height = 20)

#Compare to a more relaxed threshold
padjust_mat = df%>%mutate(cellsubtype = paste0(celltype,'_',louvain))%>%
  pivot_wider(id_cols = cellsubtype, names_from = region, values_from = padjust)%>%
  column_to_rownames('cellsubtype')
pval_mat = df%>%mutate(cellsubtype = paste0(celltype,'_',louvain))%>%
  pivot_wider(id_cols = cellsubtype, names_from = region, values_from = pval)%>%
  column_to_rownames('cellsubtype')
coef_mat[pval_mat>0.05]=0
mat = coef_mat
myBreaks <- c(seq(-0.1, 0, length.out=ceiling(50/2) + 1), 
              seq(0.1/50, 0.1, length.out=floor(50/2)))
library(pheatmap)
pheatmap(mat,show_rownames = TRUE,cluster_cols = FALSE,cluster_rows=FALSE,
         color=colorRampPalette(c("#219ebc","white","#900603"))(50),
         border_color=NA,breaks=myBreaks,
         filename = '~/output/hooke_result/all_subtype_relaxedthreshold.pdf',width = 10,height = 20)

#Plot counts
df_all%>%filter(count>=50)%>%ggplot(.,aes(x=region,y=count))+geom_boxplot()
count_mat = df_all%>%mutate(cellsubtype = paste0(celltype,'_',louvain))%>%
  pivot_wider(id_cols=cellsubtype,names_from = region,values_from = count)%>%column_to_rownames('cellsubtype')
