library(monocle3)
library(hooke)
library(purrr)
library(tidyverse)
library(Matrix)
library(splines)

data_dir = "/scratch/nsnyderm/u01/intermediate_files/regions_monocle/"
files = list.files(data_dir)

df_all = NULL
counts_all=NULL


for(file in files){
    region = gsub('.rds','',file)
    cds = readRDS(paste0(data_dir,file))
    print(dim(cds))
    
    #cds_subset=cds[,colData(cds)$age>3] #subset to non-infants
    cds_subset = cds
    
    colData(cds_subset)$cell_class_annotation = as.character(colData(cds_subset)$cell_class_annotation)
    colData(cds_subset)$animal_id = as.character(colData(cds_subset)$animal_id)
    colData(cds_subset)$sex <- as.integer(colData(cds_subset)$sex == "F")
    sample_size = length(unique(colData(cds_subset)$animal_id))

    
    print(region)
    print(dim(cds_subset))
    
    
    ccs = new_cell_count_set(cds_subset, 
                             sample_group = "animal_id", 
                             cell_group = "cell_class_annotation") 
    colData(ccs)$age = as.numeric(colData(ccs)$age)
    
    
    ccm = new_cell_count_model(ccs, 
                               main_model_formula_str = "~age+sex",vhat_method="bootstrap",
                               penalize_by_distance=F,
                               num_bootstraps = 100
                               )
    
    #test stat (for significant testing)
    estimate = coef(ccm@full_model_family)
    se = standard_error(ccm@full_model_family,type='bootstrap')
    t_stat = estimate/se
    pval =  2 * pt(abs(t_stat), df=sample_size-2, lower.tail=F)
    
    df_coef = t(estimate)%>%data.frame%>%rename_with(~ paste0(., "_coef"))%>%rownames_to_column('celltype')
    df_se = t(se)%>%data.frame%>%rename_with(~ paste0(., "_se"))%>%rownames_to_column('celltype')
    df_pval= t(pval)%>%data.frame%>%rename_with(~ paste0(., "_pval"))%>%rownames_to_column('celltype')
    df_size = colData(cds_subset)%>%data.frame()%>%mutate(celltype =cell_class_annotation)%>%group_by(celltype)%>%summarise(count=n())
    
    df_combined <- Reduce(function(x, y) merge(x, y, by = "celltype"), list(df_coef, df_se, df_pval, df_size))
    df_combined$age_padjust = p.adjust(df_combined$age_pval, method = 'fdr', n = dim(df_combined)[1])
    df_combined$sex_padjust = p.adjust(df_combined$sex_pval, method = 'fdr', n = dim(df_combined)[1])
    
    df_combined$region = region
    df_all = rbind(df_all, df_combined)
    
    #Save normalized counts
    sel_ccs_counts = normalized_counts(ccs, norm_method="size_only", pseudocount=0) #normalize counts by sizefactor
    sel_ccs_counts_long = tibble::rownames_to_column(as.matrix(sel_ccs_counts) %>% as.data.frame, var="cell_group") %>%
      pivot_longer(!cell_group, names_to="animal", values_to="num_cells")
    
    emb_time = colData(ccs)%>%as.data.frame%>%select('age')%>%rownames_to_column(var='animal')
    sel_ccs_counts_long = left_join(sel_ccs_counts_long,emb_time,by='animal')
    
    sel_ccs_counts_long$region = region
    counts_all = rbind(counts_all,sel_ccs_counts_long)
}    
    


write.csv(df_all,'/home/wyang114/output/hooke_result_final/cell-class_full_linear_ctrl-sex.csv')
write.csv(counts_all,'/home/wyang114/output/hooke_result_final/cell-class_full_linear_ctrl-sex_counts.csv')
