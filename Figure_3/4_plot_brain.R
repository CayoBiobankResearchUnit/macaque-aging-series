#' @file 3_plot_brain.R
#' @title Plot Abundance Coefficients Across Brain Regions
#' @author Wei Yang
#'
#' @description
#' Generates brain maps colored by cell abundance posterior means (pm) from \code{mashr}.
#'
#' @section Inputs:
#' \itemize{
#'   \item \code{df}: Data frame containing \code{mashr} outputs (\code{lfsr}, \code{pm}, \code{psd}).
#' }
#'
#' @section Outputs:
#' \itemize{
#'   \item Brain map visualizations colored by cell abundance.
#' }
#'
#' @usage
#' source("3_plot_brain.R")  # Or run line-by-line in RStudio
NULL

## Packages to install
# https://github.com/coolbutuseless/svgparse

devtools::install_github("kchiou/macaquesimplebrain")
library(dplyr)
library(macaquesimplebrain)
library(RColorBrewer)

# EXAMPLES
all.regions = c('ACC','CN','dlPFC','EC','HIP','IPP','lCb','M1','MB','mdTN','NAc')
manual.colors = c(ACC = "#5555b0", NAc = "#148e34", CN = "#2caa4d", dlPFC = "#345296",  HIP = "#84db55", lCb = "#e6cd69", IPP = "#6c9cc0", M1 = "#467aa1",  EC = "#532b70", MB = "#fc7462", mdTN = "#117867")
anatomical.colors = c(ACC = "")
# Map colors manually


# subtype coefficient
df = readr::read_csv('~/output/hooke_result_final/cell-subtype_full_linear_ctrl-sex_mash_age>3.csv')
df$logfoldchange = log(exp(df$age_mashr_pm))
df = df%>%filter(age_mashr_lfsr<0.05)
id='opc-olig_7'
subtype_df = df%>%filter(cell_subtype==id)
#threshold = max(abs(min(subtype_df$logfoldchange)),abs(max(subtype_df$logfoldchange)))
threshold=0.1
plot_brain(values=subtype_df$logfoldchange,
           color.family='brewer',color.theme='RdBu',#reverse=TRUE,
           regions=subtype_df$region,legend=T,labels=T,color.limits = c(-threshold,threshold),
           legend.title='',
           display.options=list(stroke_width='1px')
)
