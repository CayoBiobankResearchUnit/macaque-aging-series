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
library(RColorBrewer)
library(msigdbr)


set.seed(2025)

data_dir = '/scratch/wyang114/outputs/transcriptional/GO_terms/subtype_region_gsea_C5/'
files = list.files(data_dir)
all_df = NULL
for(file in files[grep('.csv',files)]){
  temp = read.csv(file.path(data_dir,file))
  if('X' %in% names(temp)){
    temp$X=NULL
  }
  all_df = bind_rows(all_df,temp)
}
write.csv(all_df,'/scratch/wyang114/outputs/transcriptional/GO_terms/gsea_C5_summary_opc-split.csv',row.names = F)


all_df = read.csv('/scratch/wyang114/outputs/transcriptional/GO_terms/gsea_C5_summary_opc-split.csv')
ct_list = c('astrocytes','microglia','oligodendrocytes','ependymal cells','vascular cells',
            'glutamatergic neurons','GABAergic neurons',
            'medium spiny neurons','cerebellar neurons','midbrain neurons','basket cells')
region_level=c('ACC','dlPFC','IPP','M1','EC','HIP','CN','NAc','mdTN','MB','lCb')

#Filter unrelated pathways
#select_terms = data_all$Description[-grep('PANCREAS|BILE|SPERMATO|ESTROGEN|ALLOGRAFT|ANDROGEN',data_all$Description)]
#all_df = all_df%>%filter(Description %in% select_terms)

#Evaluate consistency among cell types
all_df = all_df%>%
  mutate(
    tags = str_extract(leading_edge, "(?<=tags=)\\d+"),
    list = str_extract(leading_edge, "(?<=list=)\\d+"),
    signal = str_extract(leading_edge, "(?<=signal=)\\d+")
  )%>%
  mutate(direction = case_when(
    NES>0 ~ 'up',
    NES<0 ~ 'down'
  ))%>%
  mutate(sig_score = case_when(
    direction=='up' ~ -log10(pvalue),
    direction=='down' ~ log10(pvalue)
  ))%>%
  # group_by(region,cell_type,direction)%>%
  # top_n(abs(sig_score),n = 50)%>%
  #ungroup()%>%
  mutate(id = paste(region,cell_type,sep='_'))%>%
  group_by(Description,id)%>%
  mutate(n_dir=n_distinct(direction))%>%
  ungroup()

#plot consistency
all_df%>%mutate(consistency=case_when(n_dir>1 ~ 'inconsistent',n_dir==1 ~ 'consistent'))%>%
  ggplot(.,aes(x=interaction(region,cell_type),fill=consistency))+geom_bar()+
  ylab('# of terms with subtype consistency')+coord_flip()

# #################anova
# all_df = 
#   select_all%>%
#   filter(p.adjust<0.05)%>%
#   mutate(direction = case_when(
#     NES>0 ~ 'up',
#     NES<0 ~ 'down'
#   ))%>%
#   mutate(sig_score = case_when(
#     direction=='up' ~ -log10(pvalue),
#     direction=='down' ~ log10(pvalue)
#   ))
# 
# go_mat=all_df%>%
#   mutate(id = paste(region,cell_type,louvain,sep='_'))%>%
#   dplyr::select(c(Description, id, sig_score))%>%
#   pivot_wider(id_cols = 'Description',names_from = 'id',values_from = 'sig_score',values_fill = list(sig_score = 0))%>%
#   column_to_rownames('Description')
# 
# #ANOVA on region, cell type and subtype
# data = data.frame(t(go_mat))
# data$region = str_split_fixed(colnames(go_mat),'_',3)[,1]
# data$celltype = str_split_fixed(colnames(go_mat),'_',3)[,2]
# data$subtype = paste(str_split_fixed(colnames(go_mat),'_',2)[,2],sep='_')
# anova_results = NULL
# for(covariate in c('region','celltype','subtype')){
#   print(covariate)
#   temp_results <- sapply(colnames(data)[!(colnames(data) %in% c("region", "celltype","subtype"))], function(feature) {
#     fit <- aov(as.formula(paste(feature, sprintf("~ %s",covariate))), data = data)
#     if(covariate=='subtype'){
#       fit <- aov(as.formula(paste(feature, '~ celltype + subtype')), data = data)
#     }
#     anova_summary <- summary(fit)[[1]]
#     # Extract F-statistic and p-value for both A and B
#     c(F_value = anova_summary[covariate, "F value"],
#       P_value = anova_summary[covariate, "Pr(>F)"])
#   })
#   row.names(temp_results)=gsub('value',covariate,row.names(temp_results))
#   anova_results=rbind(anova_results,temp_results)
# }
# 
# anova_df <- as.data.frame(t(anova_results))
# anova_df$Feature <- rownames(anova_df)
# anova_df%>%ggplot(.,aes(x=F_region,y=F_celltype))+geom_point()
# anova_df%>%ggplot(.,aes(x=F_region,y=F_subtype))+geom_point()
# 
# hist(log10(anova_df$F_region),breaks=50)
# hist(log10(anova_df$F_celltype),breaks=50)
# hist(log10(anova_df$F_subtype),breaks=50) 
# 
# specific_terms = anova_df%>%filter(P_region<0.05|P_celltype<0.05|P_subtype<0.05)#%>%
# #  filter(F_region>3|F_celltype>10|F_subtype>5)%>%pull(Feature)
# specific_terms = anova_df%>%filter(F_region>5)%>%pull(Feature)
# specific_terms = gsub('\\.',' ',specific_terms)
# ####################################



go_mat = all_df%>%
  #filter(Description %in% select_term)%>%
  filter(setSize>50)%>% 
  group_by(Description,cell_type,region)%>%
  summarise(mean_sig_score=mean(sig_score),.groups = "drop")%>%
  mutate(id = paste(region,cell_type,sep='_'))%>%
  dplyr::select(c(Description, id, mean_sig_score))%>%
  pivot_wider(id_cols = 'Description',names_from = 'id',values_from = 'mean_sig_score')%>%
  column_to_rownames('Description')
go_mat[is.na(go_mat)]=0
dim(go_mat)

# #sorted
# category = region_level
# main_category <- str_extract(colnames(go_mat), paste(category, collapse = "|"))
# main_category_factor <- factor(main_category, levels = category)
# sorted_vector <- colnames(go_mat)[order(main_category_factor, colnames(go_mat))]
# temp = str_extract(sorted_vector, paste(category, collapse = "|"))%>%table()
# gaps = temp[category]%>%cumsum()

sorted_vector=colnames(go_mat)
go_group=10
dist_matrix <- dist(go_mat)
hc <- hclust(dist_matrix,method='ward.D')
col_clusters <- cutree(hc, k = go_group)
col_annotation <- data.frame(Cluster = factor(col_clusters))

row_annotation <- data.frame(region=str_split_fixed(colnames(go_mat),'_',2)[,1],
                             celltype = str_split_fixed(colnames(go_mat),'_',2)[,2])
rownames(row_annotation) <- colnames(go_mat)
regions = c('dlPFC','IPP','M1','ACC','EC','HIP','CN','NAc','mdTN','MB','lCb')
df_region = region_color%>%filter(identity %in% regions)
region_colors <- setNames(df_region$color, df_region$identity)
df_ct<- color_map%>%filter(identity %in% unique(all_df$cell_type))  
ct_colors <- setNames(df_ct$color, df_ct$identity) 
annotation_colors <- list(
  region = region_colors,
  celltype = ct_colors
)

breaks <- seq(-10, 10, length.out = 101)  # 51 breaks → 50 intervals
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
pdf("/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/go_all_heatmap_cluster.pdf", width = 10, height = 8)
pheatmap(go_mat,show_rownames =F,show_colnames = T,
         breaks=breaks,color = color_palette,
         clustering_method = 'ward.D',
         cluster_rows = hc,
         cluster_cols = T,#gaps_col = gaps,
         fontsize = 8,
         annotation_col = row_annotation,
         annotation_row = col_annotation,
         annotation_colors = annotation_colors, 
         treeheight_row = 15,treeheight_col = 15,
         #cutree_cols = go_group
         )
dev.off()



save(hc,col_annotation,go_mat,file = '/scratch/wyang114/outputs/transcriptional/GO_terms/gsea_go_cluster_new_opc-split.RData')


###############working progress############
load('/scratch/wyang114/outputs/transcriptional/GO_terms/gsea_go_cluster_new_opc-split.RData')
#########################
#try ordering them by aging hall marks
hallmark_table <- read_csv("/scratch/wyang114/outputs/transcriptional/GO_terms/C5_terms_filtered_grouping.csv")%>%
  dplyr::select(c(Description,aging_hallmark_refined))
expr_df <- as.data.frame(go_mat)
expr_df$Description <- rownames(go_mat)
merged_df <- left_join(expr_df, hallmark_table, by = "Description")

excluded_terms =c('synaptic functions','signaling/transduction','neuron signaling',
                  'stem cell differentiation/neuron differentiation','cell migration',
                  'neurogenesis','structural integrity','not able to be grouped')
merged_df = merged_df %>% filter(!(aging_hallmark_refined %in% excluded_terms) & !is.na(aging_hallmark_refined))
merged_df <- merged_df %>% arrange(aging_hallmark_refined)

#row annotation
row_anno <- merged_df$aging_hallmark_refined
row_df <- data.frame(value = merged_df$aging_hallmark_refined) %>%
  count(value, name = "count") %>%
  mutate(cumsum = cumsum(count))
names(row_anno) <- merged_df$Description
plot_matrix <- merged_df%>%
  dplyr::select(-Description, -aging_hallmark_refined)
rownames(plot_matrix) <- merged_df$Description


#column annotation
column_df = data.frame(columns = colnames(plot_matrix))
column_df$celltype = str_split_fixed(column_df$columns,'_',2)[,2]
column_df$region = str_split_fixed(column_df$columns,'_',2)[,1]
ct_list = c('astrocytes','microglia','oligodendrocytes','ependymal cells','vascular cells',
            'glutamatergic neurons','GABAergic neurons',
            'medium spiny neurons','midbrain neurons','cerebellar neurons','basket cells')
region_level=c('dlPFC','M1','IPP','ACC','EC','HIP','CN','NAc','mdTN','MB','lCb')
column_df$celltype = factor(column_df$celltype,levels=ct_list)
column_df$region= factor(column_df$region,levels=region_level)


column_df = column_df%>%arrange(celltype,region)
sorted_vector = column_df$columns


pdf('/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/hallmark_go.pdf',width=16,height=10)
pheatmap(as.matrix(plot_matrix)[,sorted_vector],
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_row = data.frame(AgingHallmark = row_anno),
         annotation_col =  column_df%>%column_to_rownames('columns'),
         show_rownames = FALSE,
         show_colnames = TRUE,
         breaks = seq(-10, 10, length.out = 101),
         color = colorRampPalette(c("blue", "white", "red"))(100),
         fontsize_row = 7,
         gaps_row = row_df$cumsum,
         annotation_colors = annotation_colors,
         gaps_col = cumsum(column_df$celltype%>%table))
dev.off()


##senescence in astrocytes 1
pdf('/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/astr_louvain_senescence.pdf',width=6,height=4)
senescence_term = hallmark_table[hallmark_table$aging_hallmark_refined=='cell senescence/death/neuron death/apoptosis/apoptotic process',]$Description
all_df%>%filter(cell_type=='astrocytes')%>%filter(Description %in% senescence_term)%>%
  ggplot(.,aes(as.character(louvain),sig_score))+geom_boxplot()
dev.off()
temp = all_df%>%filter(cell_type=='astrocytes')%>%filter(Description %in% senescence_term)

############################
# #focus analysis on specific groups
# gene_group = col_annotation%>%filter(Cluster %in% c(5,8))
# gene_group$classification = ''
# gene_group[grep('development|generation|neurogenesis|differentiation|projection',row.names(gene_group)),'classification']='neuron development'
# gene_group[grep('adhesion|junction',row.names(gene_group)),'classification']='cell-cell adhesion'
# gene_group[grep('membrane|axon',row.names(gene_group)),'classification']='structure'
# gene_group[(grepl('synap',row.names(gene_group))),'classification']='synapse structure'
# gene_group[(grepl('synaptic signaling',row.names(gene_group))),'classification']='synaptic signaling'
# gene_group[gene_group$classification=='','classification']='signaling'
gene_group = read.csv('/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/share_programs_go.csv')
gene_group = gene_group%>%column_to_rownames('go.term')
row.names(gene_group)[32] = 'anterograde trans-synaptic signaling'
row.names(gene_group)[35] = "regulation of trans-synaptic signaling"
row.names(gene_group)[37] = 'trans-synaptic signaling'
row.names(gene_group)[8] = "nervous system development"


subset_mat = go_mat[row.names(gene_group),]
#sorted
category = ct_list
main_category <- str_extract(colnames(subset_mat), paste(category, collapse = "|"))
main_category_factor <- factor(main_category, levels = category)
sorted_vector <- colnames(subset_mat)[order(main_category_factor, colnames(subset_mat))]
temp = str_extract(sorted_vector, paste(category, collapse = "|"))%>%table()
gaps = temp[category]%>%cumsum()
subset_mat = subset_mat[,sorted_vector]

library(ComplexHeatmap)
mat <- t(subset_mat)
gene_group = gene_group%>%arrange(classification)%>%dplyr::select(classification)
mat = mat[,rownames(gene_group)]
ha_col <- HeatmapAnnotation(df = as.data.frame(gene_group))
row_annotation <- data.frame(region=str_split_fixed(rownames(mat),'_',2)[,1],
                             celltype = str_split_fixed(rownames(mat),'_',2)[,2])
ha_row <- rowAnnotation(df = as.data.frame(row_annotation),
                        col = annotation_colors)

# Generate heatmap
pdf("/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/shared_heatmap_original.pdf", width = 10, height = 9)
ComplexHeatmap::pheatmap(
  mat,
  col = color_palette,
  breaks = seq(-10, 10, length.out = 51),
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  top_annotation = ha_col,
  left_annotation = ha_row,
  column_names_side = c("top"),
  angle_col = c("45"),
  gaps_col = cumsum(gene_group$classification%>%table)
)
dev.off()

##########DNA damage/Immune response/mTOR
#select terms
select_term = unique(all_df$Description)[grep('immune|inflammatory|cytokine',unique(all_df$Description))]
gene_group = as.data.frame(col_annotation)%>%rownames_to_column('Description')%>%
 filter(Description %in% select_term)%>%column_to_rownames('Description')
row_annotation <- data.frame(region=str_split_fixed(colnames(go_mat),'_',2)[,1],
                             celltype = str_split_fixed(colnames(go_mat),'_',2)[,2])
rownames(row_annotation) <- colnames(go_mat)

subset_mat = go_mat[row.names(gene_group),]
#ct = 'astrocytes'
#subset_mat = go_mat[row.names(gene_group),grep(ct,colnames(go_mat))]
subset_mat=subset_mat[rowSums(subset_mat)!=0,colSums(subset_mat)!=0]

category = ct_list
main_category <- str_extract(colnames(subset_mat), paste(category, collapse = "|"))
main_category_factor <- factor(main_category, levels = category)
sorted_vector <- colnames(subset_mat)[order(main_category_factor, colnames(subset_mat))]
temp = str_extract(sorted_vector, paste(category, collapse = "|"))%>%table()
gaps = temp[category]%>%cumsum()
subset_mat = subset_mat[,sorted_vector]

mat <- t(subset_mat)
mat = mat[,rownames(gene_group)]
ha_col <- HeatmapAnnotation(df = as.data.frame(gene_group))
row_annotation = row_annotation[rownames(mat),]
ha_row <- rowAnnotation(df = as.data.frame(row_annotation),
                        col = annotation_colors)

regions = c('dlPFC','IPP','M1','ACC','EC','HIP','CN','NAc','mdTN','MB','lCb')
df_region = region_color%>%filter(identity %in% regions)
region_colors <- setNames(df_region$color, df_region$identity)
df_ct<- color_map%>%filter(identity %in% unique(all_df$cell_type))  
ct_colors <- setNames(df_ct$color, df_ct$identity) 
annotation_colors <- list(
  region = region_colors,
  celltype = ct_colors
)

color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
pdf("/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/immune_heatmap_update.pdf", width = 6, height = 9)
ComplexHeatmap::pheatmap(
  mat,
  col = color_palette,
  breaks = seq(-10, 10, length.out = 51),
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  #col_annotation = ha_row,
  left_annotation = ha_row,
  cellwidth = 10,
  cellheight = 10,
  column_names_side = c("top"),
  #angle_col = c("45"),
  gaps_col = cumsum(gene_group$classification%>%table)
)
dev.off()




###############################
#Genes involved in pathway
load('/scratch/wyang114/outputs/transcriptional/GO_terms/gsea_go_cluster_new.RData')
gene_group = col_annotation%>%filter(Cluster %in% c(8,9))
#pathways = row.names(gene_group)[34] #28: synaptic membrane, 36: synaptic signaling, 34: synapse organization
all_terms = all_df%>%dplyr::select(ID,Description)%>%distinct(ID,Description)%>%column_to_rownames('Description')
df_all = read.csv('/scratch/nsnyderm/u01/intermediate_files/nebula/age/all_mashr_summary_filter_missing.csv')[-1]
ct = 'glutamatergic neurons'

for(pathways in row.names(gene_group)){
  print(pathways)
  go_id <- all_terms[pathways,]
  
  core_sets = all_df%>%filter(cell_type %in% ct)%>%filter(Description == pathways)%>%
    separate_rows(core_enrichment, sep = "/") %>%
    mutate(core_enrichment = trimws(core_enrichment))%>%filter(core_enrichment != "")%>%
    pull(core_enrichment)
  
  enriched_groups = all_df%>%mutate(id = paste(region,cell_type,louvain,sep='_'))%>%
    filter(cell_type=='glutamatergic neurons')%>%filter(Description == pathways)%>%
    pull(id)
  #compare coefficients of all genes in a pathway
  entrez_ids <- mapIds(org.Mmu.eg.db,keys = go_id,column =c("ENTREZID"),keytype = "GOALL",multiVals = "list")[[1]] #103
  ensembl_map <- AnnotationDbi::select(org.Mmu.eg.db,keys = entrez_ids,columns = c("ENSEMBL", "SYMBOL"),keytype = "ENTREZID") #103
  ensembl_map <- ensembl_map%>%filter(ENTREZID %in% core_sets) #69 in core sets
  df_long = df_all%>%
    filter(celltype=='glutamatergic neurons')%>%
    filter(gene %in% ensembl_map$ENSEMBL)%>%
    mutate(id=paste(region,cell_type,louvain,sep='_'))%>%
    filter(id %in% enriched_groups)
  
  mat = df_long%>%
    mutate(gene_id=paste(gene,external_gene_name,sep='_'))%>%
    dplyr::select(id,gene_id,age_mashr_pm,age_mashr_psd)%>%
    pivot_wider(id_cols = 'id', names_from = gene_id, values_from = age_mashr_pm/age_mashr_psd)%>%column_to_rownames('id')
  
  sum(colSums(!is.na(mat))>(0.8*nrow(mat)))/ncol(mat) #82% of the genes in pathways, 91% of the genes in pathway after filtering for core sets
  select_col = colnames(mat)[colSums(!is.na(mat))>(0.8*nrow(mat))]
  
  # color_palette <- colorRampPalette(c("blue", "white", "red"))(50)
  # breaks <- seq(-10, 10, length.out = 101)
  # pheatmap(mat[,select_col],cluster_cols = F, cluster_rows = F,show_colnames = F,
  #          breaks=seq(-0.02, 0.02, length.out = 51),color = color_palette)
  # 
  
  df_long%>%
    mutate(gene_id=paste(gene,external_gene_name,sep='_'))%>%
    dplyr::select(id,gene_id,age_mashr_pm,age_mashr_psd)%>%
    ggplot(.,aes(id,age_mashr_pm/age_mashr_psd))+geom_boxplot()
  
  # kruskal.test(age_mashr_pm ~ id, data = df_long) #p-value = 0.7749
  # kruskal.test(age_mashr_pm/age_mashr_psd ~ id, data = df_long) #p-value = 0.1134
  # aov(age_mashr_pm/age_mashr_psd ~ id, data = df_long)%>%summary() #0.207
  
  df_long = df_long%>%filter(region %in% c('ACC','dlPFC','M1','IPP','HIP','EC'))
  df_long$region=factor(df_long$region,levels=c('ACC','dlPFC','M1','IPP','HIP','EC'))
  df_long$scaled_effect = df_long$age_mashr_pm/df_long$age_mashr_psd
  
  #region
  label_df <- stack(glut_labels)
  colnames(label_df) <- c("louvain", "supergroup")
  temp = df_long%>%left_join(label_df,by='louvain')%>%
    group_by(id,region,supergroup)%>%summarise(mean_scaled_effect=mean(scaled_effect))
  res = aov(mean_scaled_effect~region+supergroup,data=temp)%>%summary()
  pval_r = res[[1]]["region", "Pr(>F)"]
  pval_r = format(signif(pval_r, digits = 3), scientific = TRUE)
  p=temp%>%
    ggplot(.,aes(x=region,mean_scaled_effect))+
    geom_boxplot(outliers = FALSE)+
    geom_jitter(size = 1)+
    ggtitle(str_wrap(paste0(pathways,',pvalue:',pval_r), width = 40))
  pathways = gsub(' ','_',pathways)
  pdf(file=sprintf('/home/wyang114/figures/pathway_scores/%s.pdf',pathways),width=3,height=2)
  print(p)
  dev.off()
  # temp%>%
  #   ggplot(.,aes(x=supergroup,mean_scaled_effect))+geom_boxplot()+geom_jitter()
  
  
}



# #rank average
# df_ranked <- df_long %>%
#   group_by(gene) %>%
#   mutate(
#     rank = rank(age_mashr_pm),  # 1-based ranks
#     norm_rank = (rank - 1) / (n() - 1)  # scale from 0 to 1
#   ) %>%
#   ungroup()
# 
# df_ranked$region = factor(df_ranked$region,levels=region_level)
# 
# df_ranked %>%
#   ggplot(.,aes(x=id,norm_rank))+geom_boxplot()
# 
# df_ranked %>% group_by(region,louvain)%>%
#   summarise(mean_rank = mean(norm_rank))%>%
#   ggplot(.,aes(x=region,mean_rank))+geom_boxplot()+geom_jitter()
# 
# 
# df_ranked %>% 
#   group_by(region,louvain)%>%
#   summarise(mean_rank = mean(norm_rank))%>%
#   left_join(label_df, by = "louvain")%>%
#   ggplot(.,aes(x=as.character(subtype),mean_rank))+geom_boxplot()+geom_jitter()
# 
# kruskal.test(norm_rank ~ id, data = df_ranked) #p-value = 0.0001726
# 
# df_ranked %>% group_by(region,louvain)%>%
#   summarise(mean_rank = mean(norm_rank))%>%
#   kruskal.test(mean_rank ~ region, data = .) #p-value = 0.002283

# library(FSA)
# res_dunn = dunnTest(norm_rank ~ id, data = df_ranked, method = "bh")
# res_df = res_dunn$res
# res_parsed <- res_df %>%
#   separate(Comparison, into = c("Group1", "Group2"), sep = " - ") %>%
#   mutate(
#     lower_group = ifelse(Z < 0, Group1, Group2),
#     higher_group = ifelse(Z < 0, Group2, Group1),
#     direction = ifelse(Z < 0, "Group1_lower", "Group2_lower")
#   )
# sig_lower <- res_parsed %>%
#   filter(P.adj < 0.05) %>%
#   group_by(lower_group) %>%
#   summarise(times_significantly_lower = n()) %>%
#   arrange(times_significantly_lower)




##plot genes
all_mat = temp%>%
  #filter(celltype=='microglia')%>%
  #filter(age_mashr_lfsr<0.05)%>%
  mutate(id=paste(region,celltype,sep='_'))%>%
  filter(id %in% colnames(go_mat))%>%
  group_by(id,external_gene_name,gene)%>%
  summarise(mean_pm=mean(age_mashr_pm,na.rm=TRUE))%>%
  ungroup()%>%
  dplyr::select(id,gene,external_gene_name,mean_pm)%>%
  pivot_wider(id_cols = id,names_from=c(external_gene_name,gene),values_from = mean_pm)%>%
  column_to_rownames('id')

p_mat = temp%>%
  #filter(celltype=='microglia')%>%
  mutate(id=paste(region,celltype,sep='_'))%>%
  filter(id %in% colnames(go_mat))%>%
  group_by(id,external_gene_name,gene)%>%
  summarise(min_padj=min(age_mashr_lfsr,na.rm=TRUE))%>%
  ungroup()%>%
  dplyr::select(id,gene,external_gene_name,min_padj)%>%
  pivot_wider(id_cols = id,names_from=c(external_gene_name,gene),values_from = min_padj)%>%
  column_to_rownames('id')


mat = data.frame(all_mat)
mat = mat[,apply(mat, 2, function(x) mean(x, na.rm = TRUE))< (-0.004)]
mat = mat[,colSums(mat<0,na.rm = T)>40]


sig_matrix <-  ifelse(is.na(p_mat), "",
                      ifelse(p_mat < 0.001, "***",
                             ifelse(p_mat < 0.01, "**",
                                    ifelse(p_mat < 0.05, "*", ""))))
sig_matrix = sig_matrix[,colnames(mat)]

bound = max(abs(mat),na.rm = T)

pheatmap(t(mat),cluster_rows = T,cluster_cols = T,fontsize = 5,
         show_rownames = T,
         clustering_method = 'ward.D2',
         display_numbers = t(sig_matrix),
         color=colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-bound, bound, length.out = 101))


gene='PLCL1'
temp%>%
  mutate(id=paste(region,celltype,sep='_'))%>%
  filter(id %in% colnames(go_mat))%>%
  filter(external_gene_name=='PLCL1',celltype=='microglia')%>%
  ggplot(.,aes(x=region,y=age_mashr_pm))+geom_violin()+geom_jitter(aes(color = as.character(louvain)))
  


###################################
####subtype level heatmap
select_region=c('HIP')
subset_mat = all_df%>%
  filter(cell_type=='glutamatergic neurons')%>%
  group_by(Description,region, cell_type)%>%
  mutate(n_dir=n_distinct(direction))%>%
  filter(n_dir>1)%>%
  ungroup()%>%
  mutate(id = paste0(region,'_',cell_type,'_',louvain))%>%
  dplyr::select(c(Description, id, sig_score))%>%
  pivot_wider(id_cols = 'Description',names_from = 'id',values_from = 'sig_score',,values_fill = list(sig_score = 0))%>%
  column_to_rownames('Description')

# main_category <- str_extract(colnames(subset_mat), paste(region_level, collapse = "|"))
# main_category_factor <- factor(main_category, levels = region_level)
# sorted_vector <- colnames(subset_mat)[order(main_category_factor, colnames(subset_mat))]
# temp = str_extract(sorted_vector, paste(region_level, collapse = "|"))%>%table()
# gaps = temp[select_region]%>%cumsum()
sorted_vector=colnames(subset_mat)

color_palette <- colorRampPalette(c("blue", "white", "red"))(50)
breaks <- seq(-10, 10, length.out = 51)
pheatmap(subset_mat[,sorted_vector],show_rownames = T,show_colnames = T,
         breaks=breaks,color = color_palette,
         clustering_method = 'ward.D',
         cluster_rows = T,cluster_cols = F,#gaps_col = gaps,
         #annotation_row = col_annotation,
         #annotation_col = row_annotation,
         #annotation_colors = annotation_colors,
         fontsize = 8)
