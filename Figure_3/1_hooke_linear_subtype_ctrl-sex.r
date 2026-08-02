#' @file 1_hooke_linear_subtype_ctrl-sex.r
#' @title Abundance Test for Cell Subcluster
#' @author Wei Yang
#'
#' @description
#' Tests for age- and sex-associated changes in cell type abundance using 
#' \code{hooke} (Poisson log-normal model).
#'
#' @section Inputs:
#' \itemize{
#'   \item \code{region_dir}: Path to Monocle objects organized by brain region.
#'   \item \code{data_dir}: Path to Monocle objects organized by cell class.
#' }
#'
#' @section Outputs:
#' \itemize{
#'   \item \code{cell-subtype_full_linear_ctrl-sex_age.csv}: Log-fold changes (logFC), 
#'         standard errors (SE), and p-values for tested cell subtypes.
#'   \item \code{cell-subtype_full_linear_ctrl-sex_age_count.csv}: Normalized cell count 
#'         matrix across ages.
#' }
#'
#' @usage
#' Rscript 1_hooke_linear_subtype_ctrl-sex.r
NULL


library(monocle3)
library(hooke)
library(purrr)
library(tidyverse)
library(Matrix)
library(splines)

##Get sample size
region_dir = "/scratch/nsnyderm/u01/intermediate_files/regions_monocle/"
files = list.files(region_dir)
df_size = NULL
for(file in files){
  region = gsub('.rds','',file)
  cds = readRDS(paste0(region_dir,file))
  cds_subset = cds[,colData(cds)$age>3]
  sample_size = length(unique(colData(cds_subset)$animal_id))
  df_size = rbind(df_size,data.frame(region=region,sample_size=sample_size))
}
write.csv(df_size,paste0('/home/wyang114/output/hooke_result/df_sample_size_age>3.csv'))


#Testing for differential abundance with in a cell class
sample_size_df = read_csv('/home/wyang114/output/hooke_result/df_sample_size_age>3.csv')
data_dir = "/scratch/nsnyderm/u01/intermediate_files/cell-class_monocle/"
files = list.files(data_dir)
df_all = NULL
counts_all=NULL

for(file in files){
  ct = gsub('.rds','',file)
  cds = readRDS(paste0(data_dir,file))
  cds = cds[,colData(cds)$age>3]         #NOTE: comment this if running test for animals of full lifespan.
  
  regions = colData(cds)$region%>%unique()
  
  for(region in regions){
    print(region)
    cds_subset = cds[,colData(cds)$region==region]
    if(dim(cds_subset)[2]>=100){
      colData(cds_subset)$louvain = as.character(colData(cds_subset)$louvain)
      colData(cds_subset)$animal_id = as.character(colData(cds_subset)$animal_id)
      colData(cds_subset)$sex <- as.integer(colData(cds_subset)$sex == "F")
      sample_size = sample_size_df[sample_size_df$region==tolower(region),]$sample_size
      
      ccs = new_cell_count_set(cds_subset, 
                               sample_group = "animal_id", 
                               cell_group = "louvain") 
      colData(ccs)$age = as.numeric(colData(ccs)$age)
      
      
      ccm = new_cell_count_model(ccs, 
                                 #main_model_formula_str = "ns(age, df=3)",
                                 main_model_formula_str = "~age+sex", 
                                 penalize_by_distance=F,
                                 vhat_method = 'bootstrap',
                                 num_bootstraps = 100)
      
      #test stat (for significant testing)
      estimate = coef(ccm@full_model_family)
      se = standard_error(ccm@full_model_family,type='bootstrap')
      t_stat = estimate/se
      pval =  2 * pt(abs(t_stat), df=sample_size-2, lower.tail=F)
      
      df_coef = t(estimate)%>%data.frame%>%rename_with(~ paste0(., "_coef"))%>%rownames_to_column('louvain')
      df_se = t(se)%>%data.frame%>%rename_with(~ paste0(., "_se"))%>%rownames_to_column('louvain')
      df_pval= t(pval)%>%data.frame%>%rename_with(~ paste0(., "_pval"))%>%rownames_to_column('louvain')
      df_size = colData(cds_subset)%>%data.frame()%>%group_by(louvain)%>%summarise(count=n())
      
      df_combined <- Reduce(function(x, y) merge(x, y, by = "louvain"), list(df_coef, df_se, df_pval, df_size))
      df_combined$age_padjust = p.adjust(df_combined$age_pval, method = 'fdr', n = dim(df_combined)[1])
      df_combined$sex_padjust = p.adjust(df_combined$sex_pval, method = 'fdr', n = dim(df_combined)[1])
      
      df_combined$region = region
      df_combined$celltype = ct
      df_all = rbind(df_all, df_combined)
      
      
      sel_ccs_counts = normalized_counts(ccs, norm_method="size_only", pseudocount=0) #normalize counts by sizefactor
      sel_ccs_counts_long = tibble::rownames_to_column(as.matrix(sel_ccs_counts) %>% as.data.frame, var="cell_group") %>%
        pivot_longer(!cell_group, names_to="animal", values_to="num_cells")
      
      emb_time = colData(ccs)%>%as.data.frame%>%select('age')%>%rownames_to_column(var='animal')
      sel_ccs_counts_long = left_join(sel_ccs_counts_long,emb_time,by='animal')
      
      sel_ccs_counts_long$region = region
      sel_ccs_counts_long$celltype = ct
      counts_all = rbind(counts_all,sel_ccs_counts_long)
    }
  }
}

# ##results output for test of full lifep
# write.csv(df_all,'/home/wyang114/output/hooke_result_final/cell-subtype_full_linear_ctrl-sex.csv')
# write.csv(counts_all,'/home/wyang114/output/hooke_result_final/cell-subtype_full_linear_ctrl-sex_counts.csv')
##results output for test of age>3 years
write.csv(df_all,'/home/wyang114/output/hooke_result/cell-subtype_full_linear_ctrl-sex_age>3.csv')
write.csv(counts_all,'/home/wyang114/output/hooke_result/cell-subtype_full_linear_ctrl-sex_age>3_count.csv')                            
