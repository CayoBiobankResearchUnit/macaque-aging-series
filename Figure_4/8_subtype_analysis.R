source('/home/wyang114/utils/plotting_theme.R')
source('/home/wyang114/utils/annotation.R')
library(dplyr)
library(tidyverse)
library(reshape2)
library(pheatmap)
library(viridis)

df_query = read.csv('/scratch/nsnyderm/u01/intermediate_files/nebula/age/all_mashr_summary_lfsr<0.05_filter_missing.csv')[-1]
df_nebula = read.csv('/scratch/nsnyderm/u01/intermediate_files/nebula/age/nebula_NBGMM_filtered_rminfants.csv')
df_nebula$louvain = gsub('louvain','',str_split_fixed(df_nebula$id,'_',3)[,2])
df_nebula$region = str_split_fixed(df_nebula$id,'_',3)[,3]
df_nebula$celltype = df_nebula$cell_type

############################## 
#Identify subtype patterns of region-specific aDEGs
############################## 
set_name = 'midbrain'
select_region = c('MB')
ct = 'oligodendrocytes'
spec.genes = read_csv(paste0('~/output/transcriptional/region_specific_genes/',set_name,'_specific_aDEGs.csv'))
#filter spec.genes found in glut 
unique.genes = spec.genes%>%group_by(gene)%>%summarise(n=n())%>%filter(n==1)%>%pull(gene)
ct.genes = spec.genes%>%filter(gene %in% unique.genes)%>%filter(cell_type==ct)%>%pull(gene)
#calculate mean_pm of glutamatergic subtypes in cortical regions
all_louvain = df_query%>%filter(celltype==ct)%>%pull(louvain)%>%unique()
temp = df_query%>%filter(region %in% select_region)%>%
  filter(gene %in% ct.genes)%>%
  filter(celltype==ct)%>%
  group_by(louvain,gene)%>%summarise(mean_pm = mean(age_mashr_pm))%>%
  ungroup%>%
  pivot_wider(names_from = louvain, values_from = mean_pm)%>%column_to_rownames('gene')
temp[is.na(temp)]=0
louvain_none = all_louvain[!(all_louvain %in% colnames(temp))]
temp_none = data.frame(matrix(0, nrow = nrow(temp), ncol = length(louvain_none)))
colnames(temp_none) <- louvain_none
mat = cbind(temp,temp_none)

color_palette <- colorRampPalette(c("blue", "white", "red"))(50)
breaks <- seq(-0.05, 0.05, length.out = 51)
pheatmap(mat, show_rownames = FALSE, breaks = breaks, color = color_palette,
         clustering_method = "ward.D2",treeheight_row = 15,treeheight_col = 15,
         fontsize = 6)

##############################  
#Compute pearson's correlation among subtypes(jaccard similarity?)
############################## 
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

for(ct in unique(df_query$celltype)){
  subset = df_query%>%filter(celltype==ct)
  subset = subset%>%group_by(subtype,gene)%>%summarise(mean_pm = mean(age_mashr_pm))
  
  n_query = subset$subtype%>%unique%>%length()
  n_ref = n_query
  
  cor_df <- matrix(NA, nrow = n_ref, ncol = n_query)%>%data.frame()
  p_df <- matrix(NA, nrow = n_ref, ncol = n_query)%>%data.frame()
  row.names(cor_df) = unique(subset$subtype%>%unique)
  colnames(cor_df) = unique(subset$subtype%>%unique)
  row.names(p_df) = unique(subset$subtype%>%unique)
  colnames(p_df) = unique(subset$subtype%>%unique)
  overlap_size = NULL
  
  for(lv1 in unique(subset$subtype)){
    query_subset = subset%>%filter(subtype==lv1)%>%column_to_rownames('gene')
    for (lv2 in unique(subset$subtype)){
      ref_subset = subset%>%filter(subtype==lv2)%>%column_to_rownames('gene')
      overlap = intersect(unique(row.names(query_subset)),unique(row.names(ref_subset)))
      overlap_size = cbind(overlap_size,length(overlap))
      if(length(overlap)>=1){
        #res = cor.test(query_subset[overlap,"mean_pm"],ref_subset[overlap,"mean_pm"])
        #cor_df[lv1, lv2] = res$estimate
        #p_df[lv1, lv2] = res$p.value
        cor_df[lv1, lv2] = jaccard(row.names(query_subset),row.names(ref_subset))
      }
    }
  }
  
  #plot
  hist(overlap_size)
  cor_value_matrix = as.matrix(cor_df)
  cor_value_matrix = cor_value_matrix[sort(row.names(cor_value_matrix)),sort(colnames(cor_value_matrix))]
  p_value_matrix = as.matrix(p_df)
  p_value_matrix = p_value_matrix[row.names(cor_value_matrix),colnames(cor_value_matrix)]
  significance_matrix <- ifelse(cor_value_matrix < 0.3, "",
                                ifelse(p_value_matrix < 0.001, "***",
                                       ifelse(p_value_matrix < 0.01, "**",
                                              ifelse(p_value_matrix < 0.05, "*", ""))))
  significance_matrix[is.na(significance_matrix)]=""
  breaks <- seq(0, max(cor_value_matrix, na.rm = TRUE), length.out = 51)
  pheatmap(cor_value_matrix, display_numbers = significance_matrix, 
           color = viridis(50), number_color = "white",breaks=breaks,
           cluster_rows = TRUE, cluster_cols = TRUE,fontsize_row = 8,cellwidth = 10,cellheight = 10,
           treeheight_row = 10, treeheight_col = 10
  )
}

############################## 
#subtype heatmap
############################## 
ct = 'GABAergic neurons'

mat = df_query%>% #filter(p_adj_age<0.05)%>% #change filter
  filter(celltype==ct)%>%group_by(louvain,gene)%>%
  summarise(mean_pm = mean(age_mashr_pm))%>% #change logFC
  pivot_wider(names_from=louvain,values_from=mean_pm)%>%column_to_rownames('gene')
mat[is.na(mat)]=0

# Use a palette (e.g., RColorBrewer)
library(RColorBrewer)
df_labels = gaba_labels
df_regions = gaba_regions
df_lineage = gaba_lineage
annotation <- data.frame(annotation = rep(names(df_labels), lengths(df_labels)),
                      louvain = unlist(df_labels))
regions <- data.frame(regions = rep(names(df_regions), lengths(df_regions)),
                      louvain = unlist(df_regions))
lineage <- data.frame(lineage = rep(names(df_lineage), lengths(df_lineage)),
                      louvain = unlist(df_lineage))
all_df = left_join(annotation, regions,by='louvain')%>%
  left_join(.,lineage, by='louvain')%>%
  column_to_rownames('louvain')
group_palette <- c(brewer.pal(12, "Set3"), brewer.pal(3, "Paired"))  # Generate 3 colors
group_colors <- setNames(group_palette, unique(all_df$annotation))  # Map palette to groups
region_palette <- c(brewer.pal(8, "Set2"))  # Generate 3 colors
region_colors <- setNames(region_palette, unique(all_df$regions))
lineage_palette <- c(brewer.pal(3, "Spectral"))  # Generate 3 colors
lineage_colors <- setNames(lineage_palette[1:2], unique(all_df$lineage))
annotation_colors <- list(
  annotation = group_colors,
  region = region_colors,
  lineage = lineage_colors
)
color_palette <- colorRampPalette(c("blue", "white", "red"))(50)
breaks <- seq(-0.03, 0.03, length.out = 51)
ph <- pheatmap(mat,clustering_method = "ward.D",show_rownames = FALSE,
               annotation_col = all_df,annotation_colors = annotation_colors,
               cluster_cols = TRUE,breaks=breaks,col=color_palette, fontsize=6)

cell_obs = read.csv('/scratch/nsnyderm/u01/intermediate_files/cell-class_h5ad_update/all_regions_obs.csv') ##need to update the cell count after infants
cell_obs = cell_obs%>%filter(age>3)