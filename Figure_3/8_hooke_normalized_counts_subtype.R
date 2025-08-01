library(monocle3)
library(hooke)
library(purrr)
library(tidyverse)
library(Matrix)

data_dir = "/scratch/nsnyderm/u01/intermediate_files/cell-class_monocle/"
files = list.files(data_dir)

for(file in files[2:12]){
  ct = gsub('.rds','',file)
  cds = readRDS(paste0(data_dir,file))
  
  counts_all = NULL
  regions = colData(cds)$region%>%unique()
  
  for(region in regions){
    print(region)
    cds_subset = cds[,colData(cds)$region==region]
    colData(cds_subset)$louvain = as.character(colData(cds_subset)$louvain)
    colData(cds_subset)$animal_id = as.character(colData(cds_subset)$animal_id)
    
    ccs = new_cell_count_set(cds_subset, 
                             sample_group = "animal_id", 
                             cell_group = "louvain") 
    colData(ccs)$age = as.numeric(colData(ccs)$age)
    
    ccm = new_cell_count_model(ccs, 
                               main_model_formula_str = "ns(age, df=3)",
                               penalize_by_distance=F,
                               vhat_method = 'bootstrap',
                               num_bootstraps = 100)
    
    sel_ccs_counts = normalized_counts(ccs, norm_method="size_only", pseudocount=0) #normalize counts by sizefactor
    sel_ccs_counts_long = tibble::rownames_to_column(as.matrix(sel_ccs_counts) %>% as.data.frame, var="cell_group") %>%
      pivot_longer(!cell_group, names_to="animal", values_to="num_cells")

    emb_time = colData(ccs)%>%as.data.frame%>%select('age')%>%rownames_to_column(var='animal')
    sel_ccs_counts_long = left_join(sel_ccs_counts_long,emb_time,by='animal')
    
    sel_ccs_counts_long$region = region
    sel_ccs_counts_long$celltype = ct
    counts_all = rbind(counts_all,sel_ccs_counts_long)
    
  }
  write.csv(counts_all,paste0('/home/wyang114/output/hooke_result/',ct,'/',ct,'_hooke_normalized_counts.csv'))
}
