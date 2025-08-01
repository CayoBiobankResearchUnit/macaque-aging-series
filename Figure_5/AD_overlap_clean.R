source('/home/wyang114/utils/plotting_theme.R')
source('/home/wyang114/utils/annotation.R')
library(dplyr)
library(tidyverse)
library(reshape2)
library(pheatmap)
library(viridis)
library(RColorBrewer)

#Need to convert ref set of DEG to macaque names
ad.sc = read.csv('/home/wyang114/data/AD_genesets/sc_all_deg_summary.csv')[-1]
ad.sc$gene = ifelse(ad.sc$gene=='9-Sep','SEPT9',ad.sc$gene)
ad.sc[ad.sc$ct_label=='Immune_Mic','ct_label']='Mic'
ad.sc[ad.sc$ct_label=='neurons','ct_label']='EX-neurons'
ad.sc[ad.sc$ct_label=='PER.END','ct_label']='END.PER'

#add bulk
#ad.bulk = readRDS('/home/wyang114/data/AD_genesets/wan_bulk_AD_genes_mmul_orthologs.rds')

biomart = read.csv('/scratch//wyang114/data/human-macaque-orthologs/ensembl113_mmul10_macaque_human.csv')
biomart =  biomart%>%filter(Human.gene.name != '',Gene.name!='')
biomart = biomart%>%filter(Human.homology.type=='ortholog_one2one')%>%dplyr::select(Human.gene.name,Gene.name) #pick 1-to-1 ortholog at gene name level
ad.sc = left_join(ad.sc,biomart,by=c('gene'='Human.gene.name'))
colnames(ad.sc)[6] = 'macaque.gene.name' 
ad.sc = ad.sc%>%filter(macaque.gene.name!='')

####fisher's exact test

###universe of celltype
df_query = read.csv('/scratch/nsnyderm/u01/intermediate_files/nebula/age_olig_split/all_mashr_summary_filter_missing_update.csv')[-1]
df_universe = df_query %>% group_by(celltype,region,louvain,id) %>% summarise(unique_values = list(unique(gene)))
save(df_universe, file = '/scratch/nsnyderm/u01/intermediate_files/nebula/age_olig_split//nebula_universe_subtypexregion.RData')

load('/scratch/nsnyderm/u01/intermediate_files/nebula/age_olig_split//nebula_universe_subtypexregion.RData') #universe of genes in each cell type
var_info = read_delim('/data/CEM/smacklab/team_scratch/wyang114/data/var_info.txt',delim = '\t')%>%column_to_rownames('ensembl_gene_id')
df_query = read.csv('/scratch/nsnyderm/u01/intermediate_files/nebula/age_olig_split/all_mashr_summary_lfsr<0.05_filter_missing_update.csv')[-1]

#test on region vs region (celltype x region vs celltype x region)
ad.sc.all = ad.sc%>%filter((dataset=='Mathys_2024'&region!='allregions')|(dataset %in% c('Grubman_2019','Morabito_2021')))
ad.sc.all = ad.sc.all%>%filter(!(ct_label %in% c('Epithelia_End','Epithelia_Per')))
ad.sc.all$`celltype-region-dataset`= paste(ad.sc.all$ct_label,ad.sc.all$region,ad.sc.all$dataset,sep='-')

#Make some plots for the AD genes
ad.sc.all$direction = ifelse(ad.sc.all$LFC<0,'up','down')

#all cell class together
res=NULL
for(dir in c('up','down')){
  for(ct in unique(df_query$celltype)){
  df_subset = df_query%>%filter(celltype==ct)%>%
    mutate(direction=ifelse(age_mashr_pm<0,'down','up'))%>%#filter(id %in% id_list)%>%
    filter(direction==dir)
  ad.subset = ad.sc.all%>%mutate(direction=ifelse(LFC<0,'down','up'))%>%filter(direction==dir)
  universe_subset = df_universe%>%filter(celltype==ct)

  for(rg in unique(df_subset$region)){
    temp = df_subset%>%filter(region==rg)
    universe_id = universe_subset%>%filter(region==rg)%>%pull(unique_values)%>%do.call('c',.)%>%unique()
    universe = var_info[universe_id,'external_gene_name']%>%unique()
    for(rg2 in unique(ad.subset$`celltype-region-dataset`)){
      temp2 = ad.subset%>%filter(`celltype-region-dataset`==rg2)
      geneset = temp2$macaque.gene.name%>%unique() #AD genes
      deg = temp$external_gene_name #aDEGs
      a = intersect(deg, geneset)%>%length() #in deg in AD
      b = setdiff(deg, geneset)%>%length() #in deg not in AD
      c = intersect(geneset,universe)%>%length() - a #in AD not in deg
      d = setdiff(universe,union(geneset,deg))%>%length() #not AD not in deg
      overlap = intersect(deg, geneset)
      summary = fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")


      out = data.frame('con. in macaque'=paste0(ct,'_',rg),
                       'con. in AD'=rg2,
                       'direction'=dir,
                       'a'=a, 'b'=b,'c'=c,'d'=d,
                       'pvalue'=summary$p.value,
                       'OR'=summary$estimate
      )
      if(length(overlap)!=0){
        out$genelist=paste(as.character(overlap), collapse = ",")
      }else{
        out$genelist=''
      }

      res = rbind(res,out)
    }
  }
}
}
res$padj = p.adjust(p = res$pvalue,method = 'bonferroni')
write.csv(res,'/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/all_overlap.csv')
res = read.csv('/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/all_overlap.csv',row.names = 1)

regions = c('ACC','dlPFC','IPP','M1','EC','HIP','CN','NAc','mdTN','MB','lCb')
ct_list=c('astrocytes','oligodendrocytes','microglia','ependymal cells','vascular cells',
          'glutamatergic neurons','GABAergic neurons','medium spiny neurons',
          'basket cells','cerebellar neurons','midbrain neurons')


for(dir in c('up','down')){
  res_sig = res%>%mutate(OR = ifelse(padj >= 0.001, NA, log2(OR)))
  mat = res_sig%>%filter(direction==dir)%>%
    pivot_wider(id_cols =con..in.macaque,names_from = con..in.AD,values_from = OR)%>%
    column_to_rownames('con..in.macaque')
  p_mat = res_sig%>%filter(direction==dir)%>%
    pivot_wider(id_cols =con..in.macaque,names_from = con..in.AD,values_from =padj)%>%
    column_to_rownames('con..in.macaque')
  p_mat = p_mat[row.names(mat),colnames(mat)]
  significance_matrix <- ifelse(p_mat < 0.001, "***",
                                       ifelse(p_mat < 0.01, "**",
                                              ifelse(p_mat < 0.05, "*", "")))
  breaks <- seq(0, 3, length.out = 101)
  if(dir=='up'){
    color_palette = colorRampPalette(brewer.pal(9, "Blues"))(100)
  }else{
    color_palette = colorRampPalette(brewer.pal(9, "Greens"))(100)
  }
  
  df = data.frame(id = row.names(mat))
  df$celltype = str_split_fixed(df$id,'_',2)[,1]%>%factor(.,levels=ct_list)
  df$region = str_split_fixed(df$id,'_',2)[,2]%>%factor(.,levels=regions)
  df = df%>%arrange(celltype,region)
  
  df2 = data.frame(id = colnames(mat))
  df2$celltype = str_split_fixed(df2$id,'-',3)[,1]
  df2$region = str_split_fixed(df2$id,'-',3)[,2]
  df2$dataset = str_split_fixed(df2$id,'-',3)[,3]
  df2[df2$celltype %in% c('astrocytes','ASC'),'cell_type']='Ast'
  df2[df2$celltype %in% c('EX'),'cell_type']='Exc'
  df2[df2$celltype %in% c('MG','microglia'),'cell_type']='Mic'
  df2[df2$celltype %in% c('ODC','olig'),'cell_type']='Oli'
  df2[df2$celltype %in% c('neurons'),'region']='EC'
  df2 = df2%>%arrange(tolower(celltype),tolower(region),dataset)
  
  pdf(sprintf('/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/%s_overlap.pdf',dir),width=12,height=9)
  pheatmap(t(mat[df$id,df2$id]),
           #display_numbers = t(significance_matrix[df$id,df2$id]),
           breaks = breaks,main = dir,fontsize = 8,number_color = "black",
           #gaps_row = c(8,27,35,42,50,58),gaps_col = c(11,22,33,37,48,56,65,67,68,69,70),
           cluster_rows = F,cluster_cols = F,color = color_palette)
  dev.off()
}

#Split the genes in two categories: 
#1. cell type specific convergence, 2. broad program:overlap because it is broadly changing in aging and AD 
ad.sc.all[ad.sc.all$ct_label %in% c('astrocytes','ASC'),'ct_label']='Ast'
ad.sc.all[ad.sc.all$ct_label %in% c('EX','EX-neurons'),'ct_label']='Exc'
ad.sc.all[ad.sc.all$ct_label %in% c('MG','microglia'),'ct_label']='Mic'
ad.sc.all[ad.sc.all$ct_label %in% c('ODC','olig'),'ct_label']='Oli'
ad.sc.all[ad.sc.all$ct_label %in% c('endo','Epithelia_End','Epithelia_Per','Epithelia_SMC','END.PER'),'ct_label']='Epithelia'
ad.sc.all$ct_label = tolower(ad.sc.all$ct_label)
ad.sc.all$ct_label%>%table
ad.sc.all$direction = ifelse(ad.sc.all$LFC<0,'down','up')
write.csv(ad.sc.all,'/home/wyang114/data/AD_genesets/sc_all_deg_summary_region.csv')

ad.sc.all <- read.csv('/home/wyang114/data/AD_genesets/sc_all_deg_summary_region.csv')
ad_count = ad.sc.all%>%dplyr::select(macaque.gene.name,ct_label)%>%group_by(macaque.gene.name)%>%
  summarise(ad_count = n_distinct(ct_label))%>%mutate(ad_perc=ad_count/7)
hist(ad_count$ad_count)

age_count = df_query%>%dplyr::select(external_gene_name,celltype)%>%group_by(external_gene_name)%>%
  summarise(age_count=n_distinct(celltype))%>%mutate(age_perc=age_count/11)
hist(age_count$age_count)

cc_count = inner_join(ad_count,age_count,by=c('macaque.gene.name'='external_gene_name'))
unique_gene = cc_count%>%filter(ad_count==1 & age_count==1)%>%pull(macaque.gene.name)
#broad_gene = cc_count%>%filter(ad_perc>0.5 & age_perc>0.5)%>%pull(macaque.gene.name)
broad_gene = cc_count%>%filter(ad_count>1|age_count>1)%>%pull(macaque.gene.name)
perc_unique = length(unique_gene)/length(unique(cc_count$macaque.gene.name)) #0.05375769
perc_broad = length(broad_gene)/length(unique(cc_count$macaque.gene.name)) #0.1267719


res = read.csv('/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/all_overlap.csv',row.names = 1)
res$macaque_ct = str_split_fixed(res$con..in.macaque,'_',2)[,1]
res$human_ct = str_split_fixed(res$con..in.AD,'-',2)[,1]
res[res$human_ct %in% c('astrocytes','ASC','Ast'),'human_ct']='astrocytes'
res[res$human_ct %in% c('EX','EX-neurons','Exc'),'human_ct']='glutamatergic neurons'
res[res$human_ct %in% c('MG','microglia','Mic'),'human_ct']='microglia'
res[res$human_ct %in% c('Inh','INH'),'human_ct']='GABAergic neurons'
res[res$human_ct %in% c('ODC','olig','Oli','Opc','OPC'),'human_ct']= 'oligodendrocytes'
res[res$human_ct %in% c('endo','Epithelia_End','Epithelia_Per','Epithelia_SMC','END.PER'),'human_ct']='vascular cells'
res_filtered = res%>%filter(padj<0.001)  #threshold is set here

overlap_summary = res_filtered%>%
  mutate(match=ifelse(human_ct==macaque_ct,'self','others'))%>%
  dplyr::select(human_ct,macaque_ct,genelist,match)%>%
  #dplyr::select(human_ct,macaque_ct,genelist,match,direction)%>%
    separate_rows(genelist, sep = ",") %>% 
    #group_by(macaque_ct,match,direction)%>%
    group_by(macaque_ct,match)%>%
    summarise(genelist = paste(unique(genelist), collapse = ","), .groups = "drop")  # collapse back to string

overlap_summary  <- overlap_summary  %>%
  mutate(n_genes = map_int(genelist, ~ length(str_trim(str_split(.x, ",")[[1]]))))%>%
  mutate(n_unique = map_int(genelist, function(x) {
    genes <- str_split(x, ",")[[1]] %>% str_trim()%>%unique()
    sum(genes %in% unique_gene)
  }),
  per_unique=n_unique/n_genes)%>%
  mutate(n_share = map_int(genelist, function(x) {
    genes <- str_split(x, ",")[[1]] %>% str_trim()%>%unique()
    sum(genes %in% broad_gene)
  }),
  per_share=n_share/n_genes)%>%
  mutate(n_others = n_genes-n_unique-n_share,
         per_others = n_others/n_genes)


overlap_summary <- overlap_summary %>%
  rowwise() %>%
  mutate(unique_pval = binom.test(n_unique, n_genes, p = perc_unique, alternative = "greater")$p.value) %>%
  ungroup()
overlap_summary <- overlap_summary %>%
  rowwise() %>%
  mutate(share_pval = binom.test(n_share, n_genes, p = perc_broad, alternative = "greater")$p.value) %>%
  ungroup()
write.csv(overlap_summary,'/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/overlap_summary_0.001.csv')
overlap_summary <- read.csv('/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/overlap_summary_0.001.csv',row.names = 1)


#check direction of microglia genes
mgl_unique = (overlap_summary%>%filter(macaque_ct=='microglia',match=='self')%>%pull(genelist)%>%str_split(., ","))[[1]]
res_test = res%>%filter(macaque_ct==human_ct,grepl('microglia',con..in.macaque))
res_test$region = str_split_fixed(res_test$con..in.macaque,'_',2)[,2]
res_test%>%filter(direction=='up')%>%ggplot(.,aes(x=factor(region,levels=regions),y=a))%>%+geom_boxplot()+
    ylab('perc of genes overlapped')+
    ggtitle('overlap of microglia up aDEG and AD DEG')
mgl_up_unique = res_test%>%filter(direction=='up')%>%pull(genelist)%>%str_split(., ",")%>%unlist()%>%unique()
df_all= read.csv('/scratch/nsnyderm/u01/intermediate_files/nebula/age_olig_split/all_mashr_summary_filter_missing_update.csv')[-1]

# subset = df_all%>%filter(external_gene_name %in% mgl_up_unique)%>%filter(celltype=='microglia')
# mat = subset%>%group_by(region,external_gene_name)%>%summarise(mean_pm = min(age_mashr_pm))%>%
#   pivot_wider(id_cols=external_gene_name,names_from = region, values_from = mean_pm)%>%
#   column_to_rownames('external_gene_name')
# p_mat = subset%>%group_by(region,external_gene_name)%>%summarise(min_p = min(age_mashr_lfsr))%>%
#   pivot_wider(id_cols=external_gene_name,names_from = region, values_from = min_p)%>%
#   column_to_rownames('external_gene_name')
# sig_mat = ifelse(p_mat<0.05,'*','')
# 
# pdf('/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/mgl_overlap_heatmap.pdf',width = 8,height = 6)
# pheatmap(mat,breaks=seq(-0.05,0.05, length.out = 101),display_numbers = sig_mat,show_rownames = F)
# dev.off()

###########plot overlap stats
#overlap conditions with padj<0.05
#percentage overlap
df_query = read.csv('/scratch/nsnyderm/u01/intermediate_files/nebula/age_olig_split//all_mashr_summary_lfsr<0.05_filter_missing_update.csv')[-1]
n_adeg = df_query%>%group_by(celltype)%>%summarise(n_adeg=n_distinct(gene))
pdf('/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/overlap_self_others.pdf',width = 7,height = 4.5)
overlap_summary %>%
  left_join(.,n_adeg,by=c('macaque_ct'='celltype'))%>%
  complete(macaque_ct, match, fill = list(n_genes = 0))%>%
  mutate(perc_overlap=n_genes/n_adeg)%>%
  ggplot(.,aes(x=macaque_ct,y=perc_overlap,fill=match))+geom_bar(stat = 'identity',
                              position = position_dodge(width=0.8),width = 0.8)+
  scale_fill_manual(name='',values=c('#FCDC94','#78ABA8'),
                    labels = c("others" = "map to other cell class", "self" = "map to the same cell class"))+
  ylim(c(0,1))+
  xlab('')+
  ylab('percentage of aDEG overlap with AD DEGs')
dev.off()

temp = overlap_summary %>%
  left_join(.,n_adeg,by=c('macaque_ct'='celltype'))%>%
  complete(macaque_ct, match, fill = list(n_genes = 0))%>%
  mutate(perc_overlap=n_genes/n_adeg)%>%select(macaque_ct, match,perc_overlap)

#barplot
all_df = data.frame(macaque_ct=rep('all',3), match=rep('',3), 
                    category=c('per_unique','per_share','per_others'),
                    percentage=c(perc_unique,perc_broad,1-perc_unique-perc_broad))
df_long <- overlap_summary %>%select(macaque_ct,match,per_unique, per_share, per_others)%>%
  pivot_longer(cols = c(per_unique, per_share, per_others), names_to = "category", values_to = "percentage")%>%
  bind_rows(.,all_df)
pdf('/home/wyang114/output/transcriptional/transcriptional_summary/olig-split/overlap_unique_share.pdf',width = 7,height = 4.5)
df_long%>%arrange(macaque_ct)%>%
  ggplot(., aes(x = interaction(match, macaque_ct), y = percentage, fill = category)) +
  geom_bar(stat = "identity", position = 'stack',width = 0.6) +
  scale_fill_brewer(palette = 'Paired')
dev.off()

#GO
organism = "org.Mmu.eg.db"
library(organism, character.only = TRUE)
library(clusterProfiler)
library(msigdbr)
#universe = cc_count%>%pull(macaque.gene.name)%>%unique() #all overlap genes?
#universe = df_all%>%filter(celltype=='microglia')%>%pull(external_gene_name)%>%unique() 
universe = df_all%>%pull(external_gene_name)%>%unique()
geneset = broad_gene
#geneset=unique_gene[unique_gene %in% mgl_unique]
gene.df <- AnnotationDbi::select(org.Mmu.eg.db,keys = geneset,columns = c("ENTREZID", "SYMBOL"),keytype = "SYMBOL")
all.gene.df = AnnotationDbi::select(org.Mmu.eg.db,keys = universe,columns = c("ENTREZID", "SYMBOL"),keytype = "SYMBOL")
ego = enrichGO(gene.df$ENTREZID, OrgDb='org.Mmu.eg.db', ont='ALL',universe = all.gene.df$ENTREZID,pvalueCutoff=0.05)
eKEGG = enrichKEGG(gene.df$ENTREZID, organism='mcc',universe = gene.df$ENTREZID,pvalueCutoff=0.05)
barplot(ego, showCategory = 10, title = "pan-cellular gene program",font.size = 10)+
  scale_fill_gradientn(colors = colorRampPalette(brewer.pal(9, "RdYlGn"))(100))

#non-linear gene groups (EC-L2-3-IT program)
cluster_series = read.csv('/scratch/wyang114/bin/non-linear_trajectory/data/AD/EC-L2-3-IT_gene_clusters_up.csv',row.names = 1)
cluster_series$gene = str_split_fixed(row.names(cluster_series),'_',2)[,2] 
df_all = read.csv('/scratch/nsnyderm/u01/intermediate_files/nebula/age/all_mashr_summary_filter_missing_update.csv')[-1]
universe = df_all%>%filter(subtype=='glutamatergic neurons_21')%>%pull(external_gene_name)%>%unique() 
geneset=cluster_series%>%filter(cluster %in% c(2))%>%pull(gene)
gene.df <- AnnotationDbi::select(org.Mmu.eg.db,keys = geneset,columns = c("ENTREZID", "SYMBOL"),keytype = "SYMBOL")
all.gene.df = AnnotationDbi::select(org.Mmu.eg.db,keys = universe,columns = c("ENTREZID", "SYMBOL"),keytype = "SYMBOL")
ego = enrichGO(gene.df$ENTREZID, OrgDb='org.Mmu.eg.db', ont='ALL',universe = all.gene.df$ENTREZID,pvalueCutoff=0.1,
               minGSSize = 10, maxGSSize = 1000,qvalueCutoff = 0.2)
barplot(ego, showCategory = 10, title = "EC-L2-3-IT_cluster2_program",font.size = 10)

eKEGG = enrichKEGG(gene.df$ENTREZID, organism='mcc',universe = all.gene.df$ENTREZID,pvalueCutoff=0.1,
                   minGSSize = 10, maxGSSize = 1000,qvalueCutoff = 0.2)
barplot(eKEGG, showCategory = 10, title = "EC-L2-3-IT_program(cluster2,0,3,4)_kegg",font.size = 10)

c2_st = msigdbr(species = "Macaca mulatta", category = "C2")%>%filter(gs_subcollection %in% c("CP:KEGG_LEGACY","CP:REACTOME","CP:WIKIPATHWAYS"))
msigdbr_t2g = c2_st%>%dplyr::distinct(gs_name, entrez_gene) %>% as.data.frame()
msig_res <- enricher(gene          = gene.df$ENTREZID,
                     TERM2GENE     = msigdbr_t2g,
                     pvalueCutoff  = 0.05)
barplot(msig_res, showCategory = 10, title = "EC-L2-3-IT_program(cluster2)_c2",font.size = 10)


#######correlation (for those that are shared between AD dataset and our dataset, which is more correlated)
#construct a correlation dataframe with celltype x region on the x and all regions (all single cell dataset)
df_query = read.csv('/scratch/nsnyderm/u01/intermediate_files/nebula/age/all_mashr_summary_lfsr<0.05_filter_missing.csv')[-1]
df_query$`celltype-region` = paste(df_query$celltype,df_query$region,sep='-')
df = df_query%>%group_by(`celltype-region`,external_gene_name)%>%summarise(mean_pm = mean(age_mashr_pm,na.rm=TRUE))
ad.sc.all = ad.sc%>%filter(!(dataset=='Mathys_2024' & region=='allregions'))
ad.sc.all$`celltype-dataset`= paste(ad.sc.all$ct_label,ad.sc.all$region,ad.sc.all$dataset,sep='-')
ad.sc.all = ad.sc.all%>%group_by(`celltype-dataset`,macaque.gene.name)%>%summarise(mean_lfc=mean(LFC)) #for mathys, multiples of the same gene in allregions

query_set = df$external_gene_name%>%unique()
ref_set = ad.sc.all$macaque.gene.name%>%unique()
length(intersect(query_set,ref_set)) #7221/9368 genes of query in ref, 7102 after liftover, 8530/15117 genes with no lfsr filter
share_set = intersect(query_set,ref_set)
df = df%>%filter(external_gene_name %in% share_set)
ad.sc.all = ad.sc.all%>%filter(macaque.gene.name %in% share_set)

n_query = df$`celltype-region`%>%unique%>%length()
n_ref = ad.sc.all$`celltype-dataset`%>%unique%>%length()

cor_df <- matrix(NA, nrow = n_ref, ncol = n_query)%>%data.frame()
p_df <- matrix(NA, nrow = n_ref, ncol = n_query)%>%data.frame()
row.names(cor_df) = unique(ad.sc.all$`celltype-dataset`)
colnames(cor_df) = unique(df$`celltype-region`)
row.names(p_df) = unique(ad.sc.all$`celltype-dataset`)
colnames(p_df) = unique(df$`celltype-region`)
overlap_size = NULL

for(cr_query in unique(df$`celltype-region`)){
  query_subset = df%>%filter(`celltype-region`==cr_query)%>%column_to_rownames('external_gene_name')
  for(cd_ref in unique(ad.sc.all$`celltype-dataset`)){
    ref_subset = ad.sc.all%>%filter(`celltype-dataset`==cd_ref)%>%column_to_rownames('macaque.gene.name')
    overlap = intersect(unique(row.names(query_subset)),unique(row.names(ref_subset)))
    overlap_size = cbind(overlap_size,length(overlap))
    if(length(overlap)>=20){
      res = cor.test(query_subset[overlap,"mean_pm"],ref_subset[overlap,"mean_lfc"])
      cor_df[cd_ref, cr_query] = res$estimate
      p_df[cd_ref, cr_query] = res$p.value
      
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
pdf('/home/wyang114/output/transcriptional/ad_overlap/all_overlap_correlation_v2.pdf',width=14,height=14)
pheatmap(cor_value_matrix, display_numbers = significance_matrix, 
         color = viridis(50), number_color = "white",breaks=breaks,
         cluster_rows = FALSE, cluster_cols = FALSE,fontsize_row = 8,cellwidth = 10,cellheight = 10,
         gaps_row = c(8,27,35,42,50,58),gaps_col = c(11,12,13,17,26,34,36,47,48,59)
)
dev.off()

###########GWAS######
advp = read.csv('/home/wyang114/data/AD_genesets/advp.variant.records.hg38.tsv',sep='\t')
# Filter for significant AD-risk SNPs
advp <- advp %>%filter(P.value < 5e-8, Phenotype.derived == "AD")
# Get non-redundant list of genes (1 gene per locus using lead SNP) #432
ad_risk_genes <- advp %>%group_by(LocusName) %>%
  slice_min(order_by = P.value) %>%pull(nearest_gene_symb) %>%
  unique()
biomart = read.csv('/home/wyang114/data/mart_export_human_primate.txt',sep='\t')
biomart =  biomart%>%filter(Human.gene.name != '')
biomart = biomart%>%group_by(Human.gene.name,Gene.name)%>%summarise(n=n())%>%filter(n==1) #pick 1-to-1 ortholog at gene name level
ad.gwas = biomart[biomart$Human.gene.name %in% ad_risk_genes,] #408/432
colnames(ad.gwas)[2] = 'macaque.gene.name' 
ad.gwas = ad.gwas%>%filter(macaque.gene.name!='')
####fisher's exact test
##which celltype x region show significant enrichment of AD set from a specfic celltype in region
#filter regional neurons/subtypes from dissection error
load('/home/wyang114/output/transcriptional/subtypexregion_percent>1.RData')

# #universe of celltype
# df_query = read.csv('/scratch/nsnyderm/u01/intermediate_files/nebula/age/all_mashr_summary_filter_missing.csv')[-1]
# df_universe = df_query %>% group_by(celltype,region,louvain,id) %>% summarise(unique_values = list(unique(gene)))
# save(df_universe, file = '/scratch/nsnyderm/u01/intermediate_files/nebula/age/nebula_universe_subtypexregion.RData')

load('/scratch/nsnyderm/u01/intermediate_files/nebula/age/nebula_universe_subtypexregion.RData')
var_info = read_delim('/data/CEM/smacklab/team_scratch/wyang114/data/var_info.txt',delim = '\t')%>%column_to_rownames('ensembl_gene_id')
df_query = read.csv('/scratch/nsnyderm/u01/intermediate_files/nebula/age/all_mashr_summary_lfsr<0.05_filter_missing.csv')[-1]

#test for enrichment of AD risk genes
res=NULL
for(ct in unique(df_query$cell_type)){
  df_subset = df_query%>%filter(celltype==ct)%>%mutate(direction=ifelse(age_mashr_pm<0,'down','up'))#%>%filter(id %in% id_list)
  universe_subset = df_universe%>%filter(celltype==ct)%>%
    mutate(id=paste(celltype,paste0('louvain',louvain),region,sep='_'))#%>%filter(id %in% id_list)
  for(dir in c('up','down')){
    for(rg in unique(df_subset$region)){
      temp = df_subset%>%filter(region==rg)%>%filter(direction==dir)
      universe_id = universe_subset%>%filter(region==rg)%>%pull(unique_values)%>%do.call('c',.)%>%unique()
      universe = var_info[universe_id,'external_gene_name']%>%unique()
      temp2 = ad.gwas
      geneset = temp2$macaque.gene.name%>%unique()
      deg = temp$external_gene_name
      a = intersect(deg, geneset)%>%length()
      b = setdiff(deg, geneset)%>%length()
      c = intersect(geneset,universe)%>%length() - a
      d = setdiff(universe,union(geneset,deg))%>%length()
      summary = fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
        
        
        out = data.frame('region in macaque'=rg,
                         'region in AD'='all-GWAS', 
                         'direction'=dir,
                         'a'=a, 'b'=b,'c'=c,'d'=d,
                         'pvalue'=summary$p.value,
                         'OR'=summary$estimate,
                         'cell_type'=ct
        )
        res = rbind(res,out)
      }
  }
}

write.csv(res,sprintf('/home/wyang114/output/transcriptional/ad_overlap/all_overlap_gwas.csv'))
regions = c('ACC','dlPFC','IPP','M1','EC','HIP','CN','NAc','mdTN','MB','lCb')
ct_list=c('astrocytes','oligodendrocytes','microglia','ependymal cells','vascular cells',
          'glutamatergic neurons','GABAergic neurons','medium spiny neurons',
          'basket cells','cerebellar neurons','midbrain neurons')
res$padj = p.adjust(p = res$pvalue,method = 'bonferroni')

mat = res%>%mutate(deg_set = paste0(cell_type,'_',direction))%>%
  pivot_wider(id_cols = deg_set,names_from = region.in.macaque,values_from = OR)%>%
  column_to_rownames('deg_set')
p_mat = res%>%mutate(deg_set = paste0(cell_type,'_',direction))%>%
  pivot_wider(id_cols = deg_set,names_from = region.in.macaque,values_from = padj)%>%
  column_to_rownames('deg_set')
significance_matrix <- ifelse(is.na(p_mat), "",
                            ifelse(p_mat < 0.001, "***",
                              ifelse(p_mat < 0.01, "**",
                                     ifelse(p_mat < 0.05, "*", ""))))
category = ct_list
main_category <- str_extract(row.names(mat), paste(category, collapse = "|"))
main_category_factor <- factor(main_category, levels = category)
sorted_vector <- row.names(mat)[order(main_category_factor, row.names(mat))]


pdf(sprintf('/home/wyang114/output/transcriptional/ad_overlap/all_overlap_gwas_bonferonni.pdf'),width=6,heigh=4)
breaks <- seq(0, max(mat,na.rm = T), length.out = 101)
#color_palette = colorRampPalette(brewer.pal(9, "YlGn"))(100)
color_palette = viridis(100)
pheatmap(mat[sorted_vector,regions],display_numbers = significance_matrix[sorted_vector,regions],
         breaks = breaks,main='gwas_overlap_AD',
         cluster_rows = F,cluster_cols = F,gaps_row = c(10),color=color_palette,border_color = NA)
dev.off()

#glut subtype
ct = 'glutamatergic neurons'
subset_query = df_query%>%filter(celltype==ct)
res=NULL
for(lv in unique(subset_query$louvain)){
  df_subset = subset_query%>%filter(louvain==lv)%>%mutate(direction=ifelse(age_mashr_pm<0,'down','up'))%>%filter(id %in% id_list)
  universe_subset = df_universe%>%
    filter(celltype==ct,louvain==lv)%>%
    mutate(id=paste(celltype,paste0('louvain',louvain),region,sep='_'))%>%
    filter(id %in% id_list)
  for(dir in c('up','down')){
    for(rg in unique(df_subset$region)){
      temp = df_subset%>%filter(region==rg)%>%filter(direction==dir)
      universe_id = universe_subset%>%filter(region==rg)%>%pull(unique_values)%>%do.call('c',.)%>%unique()
      universe = var_info[universe_id,'external_gene_name']%>%unique()
      temp2 = ad.gwas
      geneset = temp2$macaque.gene.name%>%unique()
      deg = temp$external_gene_name
      a = intersect(deg, geneset)%>%length()
      b = setdiff(deg, geneset)%>%length()
      c = intersect(geneset,universe)%>%length() - a
      d = setdiff(universe,union(geneset,deg))%>%length()
      summary = fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
      
      
      out = data.frame('region in macaque'=rg,
                       'region in AD'='all-GWAS', 
                       'direction'=dir,
                       'a'=a, 'b'=b,'c'=c,'d'=d,
                       'pvalue'=summary$p.value,
                       'OR'=summary$estimate,
                       'cell_type'=ct,
                       'louvain'=lv
      )
      res = rbind(res,out)
    }
  }
}

write.csv(res,'/home/wyang114/output/transcriptional/ad_overlap/glut_overlap_gwas.csv')

res = read.csv('/home/wyang114/output/transcriptional/ad_overlap/glut_overlap_gwas.csv')[-1]
res$padj = p.adjust(p = res$pvalue,method = 'bonferroni')
res = res%>%filter(a>=10)

mat = res%>%mutate(deg_set = paste0(direction,'_',louvain))%>%
  pivot_wider(id_cols = deg_set,names_from = region.in.macaque,values_from = OR)%>%
  column_to_rownames('deg_set')
p_mat = res%>%mutate(deg_set = paste0(direction,'_',louvain))%>%
  pivot_wider(id_cols = deg_set,names_from = region.in.macaque,values_from = padj)%>%
  column_to_rownames('deg_set')
significance_matrix <- ifelse(is.na(p_mat), "",
                              ifelse(p_mat < 0.001, "***",
                                     ifelse(p_mat < 0.01, "**",
                                            ifelse(p_mat < 0.05, "*", ""))))
# Split prefix and numeric suffix
order = c(5,0,17,27,29,16,24,30,28,18,26,6,20,1,12,23,3,22,2,19,21,4,25,8,9,10,11,13,14,15)%>%as.character()
df <- data.frame(original = row.names(mat)) %>%
  mutate(prefix = str_extract(original, "^[^_]+"),
         suffix = as.numeric(str_extract(original, "[^_]+$")))%>%
  mutate(prefix = factor(prefix, levels = c('up','down')),
         suffix = factor(suffix,levels=order)) %>%
  arrange(prefix, suffix)
ordered_vec <- df$original

pdf(sprintf('/home/wyang114/output/transcriptional/ad_overlap/glut_overlap_gwas_bonferonni_overlap>10.pdf'),width=6,height=8)
breaks <- seq(0, max(mat,na.rm = T), length.out = 101)
#color_palette = colorRampPalette(brewer.pal(9, "YlGn"))(100)
color_palette = viridis(100)
pheatmap(mat[ordered_vec ,],display_numbers = significance_matrix[ordered_vec ,],
         breaks = breaks,main='gwas_overlap_AD',
         cluster_rows = F,cluster_cols = F,gaps_row = c(27),color=color_palette,border_color = NA)
dev.off()



#some example genes in with really high genetic risk
df_query = read.csv('/scratch/nsnyderm/u01/intermediate_files/nebula/age/all_mashr_summary_filter_missing.csv')[-1]
df_query%>%filter(external_gene_name=='APOE')%>%mutate(sig = ifelse(age_mashr_lfsr<0.05,'S','NS'))%>%
  ggplot(.,aes(x=region,y=age_mashr_pm,color=sig))+geom_boxplot()+facet_wrap(~celltype)+
    geom_hline(yintercept=0,linetype = 'dotted')+ggtitle('APOE-linear')

df_query%>%filter(external_gene_name=='APP')%>%mutate(sig = ifelse(age_mashr_lfsr<0.05,'S','NS'))%>%
  ggplot(.,aes(x=region,y=age_mashr_pm,color=sig))+geom_boxplot()+facet_wrap(~celltype)+
  geom_hline(yintercept=0,linetype = 'dotted')+ggtitle('APP-linear')

df_query%>%filter(external_gene_name=='PSEN1')%>%mutate(sig = ifelse(age_mashr_lfsr<0.05,'S','NS'))%>%
  ggplot(.,aes(x=region,y=age_mashr_pm,color=sig))+geom_boxplot()+facet_wrap(~celltype)+
  geom_hline(yintercept=0,linetype = 'dotted')+ggtitle('PSEN1-linear')

df_query%>%filter(external_gene_name=='RELN')%>%mutate(sig = ifelse(age_mashr_lfsr<0.05,'S','NS'))%>%
  ggplot(.,aes(x=region,y=age_mashr_pm,color=sig))+geom_boxplot()+facet_wrap(~celltype)+
  geom_hline(yintercept=0,linetype = 'dotted')+ggtitle('RELN-linear')

label_df <- stack(glut_labels)
colnames(label_df) <- c("louvain", "supergroup")

df_query%>%filter(external_gene_name %in% c('RELN','APOE'))%>%
  pivot_wider(id_cols = c(region,celltype,louvain),names_from = external_gene_name, values_from = age_mashr_pm)%>%
  ggplot(.,aes(APOE,RELN,colour = celltype))+geom_point()+
  geom_vline(xintercept = 0,linetype='dotted')+geom_hline(yintercept = 0,linetype='dotted')

df_query%>%filter(external_gene_name %in% c('RELN','APOE'))%>%filter(region=='EC')%>%
  filter(celltype=='glutamatergic neurons')%>%left_join(.,label_df,by='louvain')%>%
  pivot_wider(id_cols = c(region,celltype,louvain,supergroup),names_from = external_gene_name, values_from = age_mashr_pm)%>%
  ggplot(.,aes(APOE,RELN,colour = supergroup))+geom_point()+
  geom_vline(xintercept = 0,linetype='dotted')+geom_hline(yintercept = 0,linetype='dotted')+
  scale_color_manual(values = c(brewer.pal(12, "Paired"),brewer.pal(6, "Set2")))
  
