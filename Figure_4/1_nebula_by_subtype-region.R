#' @file 1_nebula_by_subtype-region.R
#' @title Differential expression test across subclusters x region
#' @author Wei Yang
#'
#' @description
#' Tests for age- and sex-associated gene expression change with
#' \code{nebula} (negative binomial mixed model) controlling for log2umi, sex and age as covariate. 
#' Subcluster x region with <30 animals were excluded.
#'
#' @section Inputs:
#' \itemize{
#'   \item \code{data_dir}: Path to data matrix (mtx.gz), feature data frame (var.txt.gz) 
#'   and metadata (obs.txt.gz) organized by cell class.
#' }
#'
#' @section Outputs:
#' \itemize{
#'   \item \code{<cell class>_<subcluster>.txt}: data frame of nebula output with logFC, se and p-value, 
#'   adjusted p-value after FDR correction for each subcluster.
#' }
#'
#' @usage
#' #SBATCH --mem=64G
#' #SBATCH --cpus-per-task=20
#' #SBATCH --time=4:00:00
#' module load gsl-2.7.1-gcc-11.2.0
#' module load r-4.2.2-gcc-11.2.0
#' module load gcc-11.2.0-gcc-11.2.0
#' Rscript 1_nebula_by_subtype-region.R $cell_type $SLURM_CPUS_ON_NODE
NULL

arguments = commandArgs(trailing=TRUE)

ct = arguments[1]

if (length(arguments) > 1) {
  ncpus = as.integer(arguments[2])
} else {
  ncpus = 1
}


library(data.table)
library(nebula)
library(Matrix)
library(tidyverse)
library(future)
library(rngtools)
library(foreach)
library(purrr)
set.seed(2025)

data_dir = paste0('/scratch/wyang114/outputs/nebula/prepared_data/',ct)
files = list.files(data_dir)
files = files[grep('mtx',files)]
louvains = str_split_fixed(files,'_',3)[,2]
output_dir = paste0('/scratch/wyang114/outputs/nebula/by_region_NBLMM/',ct)
dir.create(file.path(output_dir,ct))

# Read in files 
for(n in louvains){
  
  mat = readMM(sprintf('%s/%s_%s_exp.mtx', data_dir,ct,n))
  var = read.table(sprintf('%s/%s_%s_var.txt.gz', data_dir,ct,n), sep = "\t", header=TRUE, row.names=1)
  obs = read.table(sprintf('%s/%s_%s_obs.txt.gz', data_dir,ct,n), sep = "\t", header=TRUE, row.names=1)
  
  all_summary=NULL
  
  for(region in unique(obs$region)){
    start.time <- Sys.time()  
    tryCatch({
      e.keep = mat[,obs$region==region]
      print(dim(e.keep))
      obs.keep = obs[obs$region==region,]
      e.keep@Dimnames = list(var$ensembl_gene_id, rownames(obs.keep))
      # Sort everything (required for nebula)
      e.keep = e.keep[,order(obs.keep$animal_id)]
      obs.keep = obs.keep[order(obs.keep$animal_id),]
      
      if (length(unique(obs.keep$animal_id))<=30) {
        print(paste("Skipping regions with <=30 animals:", ct, "_", region, "_", n))
        next
      }
      
      sample_data = all.data = list(
        count = e.keep,
        pred = model.matrix(as.formula(paste0('~log2umi + age + sex')), data=obs.keep),
        id = obs.keep$animal_id,
        offset = obs.keep$n_umi
      )
      # Run nebula 
      re = nebula(sample_data$count, sample_data$id, pred=sample_data$pred, offset=sample_data$offset,model='NBGMM', ncore = ncpus)
      
      # Clean output to only what we care about, add external gene names 
      re_summary = re$summary
      re_summary$convergence = re$convergence
      if(is.null(dim(re$overdispersion))){
        re_summary$subject_overdispersion = re$overdispersion
      }else{
        re_summary$subject_overdispersion = re$overdispersion[,1]
        re_summary$cell_overdispersion = re$overdispersion[,2]
      }
      re_summary <- left_join(re_summary, var, by = c("gene" = "ensembl_gene_id"))
      re_summary = re_summary %>% select(logFC_log2umi, logFC_age, se_log2umi, se_sexM, se_age, p_log2umi, p_sexM, p_age, gene, external_gene_name, chr)
      
      # Adjust p-values using FDR 
      re_summary$p_adj_age = p.adjust(re_summary$p_age, method='fdr')
      re_summary$p_adj_sex = p.adjust(re_summary$p_sex, method='fdr')
      re_summary$p_adj_log2umi = p.adjust(re_summary$p_log2umi, method='fdr')
      
      #add region info
      re_summary$region=region
      re_summary$cell_type=ct
      re_summary$louvain=n
      re_summary = re_summary%>%mutate(id = paste(celltype, louvain, region, sep = "_"))
      all_summary = rbind(all_summary,re_summary)
    }, error = function(e) {
    message(paste("Skipping:", ct, "_", region, "_", n))
  })
    
    end.time <- Sys.time()
    print(end.time-start.time)
  }
  # Write output
  write_csv(all_summary, sprintf('%s/%s/%s_%s.txt', output_dir,ct,ct,n))
    
}





