source('/home/wyang114/utils/plotting_theme.R')
source('/home/wyang114/utils/annotation.R')
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

#df_all = read.csv('/scratch/nsnyderm/u01/intermediate_files/nebula/age_olig_split//all_mashr_summary_filter_missing.csv')[-1]
df=read.csv('/scratch/nsnyderm/u01/intermediate_files/nebula/age_olig_split/all_mashr_summary_lfsr<0.05_filter_missing_update.csv')[-1]
all_ct = unique(df$celltype)[c(1,4,5,6,7,9,10,11)]
regions = c('ACC','dlPFC','IPP','M1','EC','HIP','CN','NAc','mdTN','MB','lCb')
df$region = factor(df$region, levels=regions)

select_region = list(c('ACC','dlPFC','IPP','M1'),'EC','HIP',c('CN','NAc'),'MB','mdTN','lCb')
set_name = c('cortical','entorhinal','hippocampal','basal ganglia','midbrain','thalamas','cerebellum')

df_unique=NULL
for(i in 1:length(all_ct)){
    ct=all_ct[i]
    mat = df%>%filter(celltype==ct)%>%group_by(region,gene)%>%summarise(mean_pm = mean(age_mashr_pm))%>%
      pivot_wider(names_from=region,values_from=mean_pm)%>%column_to_rownames('gene')
    mat[is.na(mat)]=0
    for(j in 1:length(select_region)){
      set = select_region[[j]]
      if(sum(set %in% colnames(mat))==length(set)){
        set_mat <-mat %>%rownames_to_column('gene')%>%
          rowwise() %>%
          mutate(target_gt_zero = sum(c_across(set) != 0))%>%
          mutate(nontarget_gt_zero = sum(c_across(colnames(mat)[!(colnames(mat) %in% c(set,'gene','count_gt_zero','target_gt_zero'))]) != 0))

        if(intersect(set,c('ACC','dlPFC','IPP','M1'))%>%length!=0){
          rg_share= set_mat%>%filter(nontarget_gt_zero==0 & target_gt_zero >= 3)%>%pull(gene)%>%unique()
        }else{
          rg_share= set_mat%>%filter(nontarget_gt_zero==0 & target_gt_zero == length(set))%>%pull(gene)%>%unique()
        }
        if(length(rg_share)!=0){
          temp = data.frame('gene'=rg_share,'celltype'=ct,'set'=set_name[j])
          df_unique = rbind(df_unique,temp)
        }
      }
    }
}

df_unique$set = factor(df_unique$set,levels = set_name)
df_unique%>%filter(celltype!='medium spiny neurons')%>%group_by(celltype,set)%>%summarise(nDEGs=n())%>%
  ggplot(.,aes(x=set,y=nDEGs))+
  geom_bar(stat='identity')+facet_wrap(~celltype,scales = 'free_y')
write.csv(df_unique,'/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/region-specific_summary.csv')

df_unique = read.csv('/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/region-specific_summary.csv',row.names = 1)
df_unique$set = factor(df_unique$set,levels = set_name)

temp = df%>%filter(gene %in% df_unique$gene)%>%
  filter(celltype %in% all_ct)%>%
  filter(celltype!='medium spiny neurons')%>%
  left_join(df_unique)

#number of region-specific DEGs in each celltype
pdf("/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/region-spec_deg_by-celltype.pdf", width = 4, height = 6)
temp%>%group_by(celltype,subtype,set)%>%
  filter(!is.na(set))%>%
  summarise(nDEGs=n())%>%
  ggplot(.,aes(x=set,y=nDEGs,color=set))+
  ylab('number of region-specific DEGs')+
  geom_boxplot(outlier.size = 0.5)+facet_wrap(~celltype,scales = 'free_y',ncol = 1)+
  scale_y_log10()
dev.off()

# #number of subtypes showing region-specific DEGs
# temp%>%group_by(cell_type,subtype,set)%>%
#   filter(!is.na(set))%>%
#   summarise(n_subtype = n_distinct(subtype))%>%
#   ggplot(.,aes(x=set,y=n_subtype,fill=set))+
#   ylab('number of subtypes with region-specific DEGs')+
#   geom_bar(stat='identity')+facet_wrap(~cell_type,scales = 'free_y',ncol = 1)
# 
# #correlations
# res = NULL
# for(set in set_name){
#   select_set = set
#   temp2 = df_unique%>%filter(set==select_set)
#   n_con=temp2$cell_type%>%unique%>%length
#   conditions = temp2$cell_type%>%unique
#   
#   jc_mat = matrix(NA, nrow = n_con,ncol = n_con)
# 
#   for(i in 1:n_con){
#     print(i)
#     for(j in 1:n_con){
#       print(j)
#       a = temp2%>%filter(cell_type==conditions[i])%>%pull(gene)%>%unique() #using unique genes
#       b = temp2%>%filter(cell_type==conditions[j])%>%pull(gene)%>%unique()
#       if(length(a)>0 & length(b)>0){
#         jc_mat[i,j] = jaccard(a, b)
#         out = data.frame('ct1'=conditions[i],'ct1_size'=length(a),
#                           'ct2'=conditions[j], 'ct2_size'=length(b),
#                          'overlap'=intersect(a,b)%>%length(),'jaccard'=jaccard(a, b),
#                          'set'=set)
#         res = rbind(res,out)
#       }
#     }
#   }
#   row.names(jc_mat)=conditions
#   colnames(jc_mat)=conditions
#   pheatmap(jc_mat,cluster_rows = F,cluster_cols = F,show_colnames = F,
#            breaks=seq(0, 0.2, length.out = 51),
#            cellwidth = 30, cellheight = 30,
#            main=sprintf('%s-specific DEGs',select_set))
# }
########change background
df=read.csv('/scratch/nsnyderm/u01/intermediate_files/nebula/age_olig_split/all_mashr_summary_filter_missing_update.csv')[-1]

ensembl_ids = df_unique%>%filter(set=='cortical',celltype=='glutamatergic neurons')%>%pull(gene)
# ensembl_ids = intersect(df_unique%>%filter(set=='cortical',celltype=='GABAergic neurons')%>%pull(gene),
#                 df_unique%>%filter(set=='cortical',celltype=='glutamatergic neurons')%>%pull(gene))
gene.df <- AnnotationDbi::select(org.Mmu.eg.db,keys = ensembl_ids,columns = c("ENTREZID", "ENSEMBL"),
                                 keytype = "ENSEMBL")

#universe = df%>%filter(celltype=='glutamatergic neurons')%>%pull(external_gene_name)%>%unique()
#all.gene.df = AnnotationDbi::select(org.Mmu.eg.db,keys = universe,columns = c("ENTREZID", "SYMBOL"),keytype = "SYMBOL")


ego <- enrichGO(
  gene         = gene.df$ENTREZID,
  OrgDb        = org.Mmu.eg.db,
  keyType      = "ENTREZID",
  ont          = "ALL",        # use "MF" or "CC" for other domains
  pAdjustMethod= "BH",
  pvalueCutoff = 0.2,
  qvalueCutoff = 0.2,
  #universe = all.gene.df$ENTREZID,
  readable     = TRUE         # converts ENTREZ back to gene symbols
)
barplot(ego, showCategory = 10, title = "cortical-specific DEGs (glut)",font.size = 10)



#Visualize specific genes in particular cell type
ct = 'microglia'
subset = df_unique%>%filter(cell_type==ct & set=='entorhinal')%>%
  left_join(.,df_all%>%filter(cell_type==ct),by=c('gene','cell_type')) 
mat = subset%>%dplyr::select(gene,id,age_mashr_pm)%>%
  pivot_wider(id_cols=gene,names_from=id,values_from = age_mashr_pm)%>%
  column_to_rownames('gene')
p_mat = subset%>%dplyr::select(gene,id,age_mashr_lfsr)%>%
  pivot_wider(id_cols=gene,names_from=id,values_from = age_mashr_lfsr)%>%
  column_to_rownames('gene')
sig_annot = ifelse(p_mat<0.05,'*','')
sig_annot[is.na(sig_annot)]=''

row_annot = df_unique%>%filter(cell_type==ct & set=='entorhinal')%>%
  dplyr::select(gene,set)%>%column_to_rownames('gene')
col_annot = data.frame('id'=colnames(mat),'region'=str_split_fixed(colnames(mat),'_',3)[,3],
                       'subtype'=str_split_fixed(colnames(mat),'_',3)[,2])%>%column_to_rownames('id')
sorted_cols <- order(col_annot$region)
#set_palette <- c(brewer.pal(length(unique(row_annot$set)), "Set1"))
set_colors <- setNames(set_palette, unique(row_annot$set))
regions = c('ACC','dlPFC','IPP','M1','EC','HIP','CN','NAc','mdTN','MB','lCb')
df_region = region_color%>%filter(identity %in% regions)
region_colors <- setNames(df_region$color, df_region$identity)
annotation_colors <- list(
  #set=set_colors,
  region=region_colors
)

#mat[is.na(mat)]=0
#mat = mat[,colSums(mat)!=0]
df_gene = subset%>%group_by(gene,external_gene_name)%>%summarise(n=n())%>%column_to_rownames('gene')
row.names(mat) = df_gene[row.names(mat),'external_gene_name']
row.names(row_annot) = df_gene[row.names(row_annot),'external_gene_name']

color_palette <- colorRampPalette(c("blue", "white", "red"))(50)
breaks <- seq(-0.02, 0.02, length.out = 51)
pheatmap(mat[,sorted_cols],show_rownames = T,show_colnames = F,cluster_cols = F,cluster_rows = F,fontsize = 8,
         breaks=breaks,color = color_palette,
         annotation_col = col_annot,display_numbers = sig_annot[,sorted_cols],
         #annotation_row= row_annot,
         border_color = NA,cellwidth = 5,cellheight = 15,fontsize_number = 10,number_color = 'black',
         annotation_colors = annotation_colors,clustering_method = 'ward.D')

#GO term/KEGG
temp = subset
universe = df_all%>%filter(celltype=='microglia')%>%pull(gene)%>%unique
gene_group = temp%>%pull(gene)%>%unique
gene.df <- AnnotationDbi::select(org.Mmu.eg.db,
                                 keys = gene_group,
                                 columns = c("ENTREZID", "ENSEMBL"),
                                 keytype = "ENSEMBL")
df_gene = left_join(temp,gene.df,by=c('gene'='ENSEMBL'))%>%drop_na()
all.gene.df = AnnotationDbi::select(org.Mmu.eg.db,
                                    keys = universe,
                                    columns = c("ENTREZID", "ENSEMBL"),
                                    keytype = "ENSEMBL")
compKEGG = compareCluster(df_gene$ENTREZID~set, data = df_gene,
                          fun = "enrichKEGG",organism='mcc',universe = all.gene.df$ENTREZID,pvalueCutoff=0.05)
df_KEGG = data.frame(compKEGG)
print(dim(df_KEGG))
if(dim(df_KEGG)[1]!=0){
  df_KEGG$ONTOLOGY = 'KEGG'
  df_KEGG$Description = gsub(" - Macaca mulatta \\(rhesus monkey\\)", "", df_KEGG$Description)
}
compGO = compareCluster(df_gene$ENTREZID~set, data = df_gene,
                        fun = "enrichGO", OrgDb='org.Mmu.eg.db', ont='ALL',universe = all.gene.df$ENTREZID,pvalueCutoff=0.05)
df_GO = data.frame(compGO)
print(dim(df_GO))
df_combine = rbind(df_KEGG,df_GO)
