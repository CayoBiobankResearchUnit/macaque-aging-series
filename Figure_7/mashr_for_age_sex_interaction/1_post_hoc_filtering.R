library(tidyverse)

add = read_csv('/scratch/klwatki5/R_analysis/bbi/mashr/mash_new/results/master_merged_update.csv')
add = add %>%  mutate(gene.id = paste0(ensembl_gene_id,"_",id))

int = read_csv('/scratch/klwatki5/R_analysis/bbi/nebula/interaction/raw_combined_results.csv')
int = int %>% mutate(cell_type = case_when(
  (cell_type == "oligodendrocytes" & louvain %in% c(12,13)) ~ "OPC",
  TRUE ~ cell_type)) %>%
  mutate(id = paste0(cell_type,"_louvain",louvain,"_", region)) %>%
  mutate(gene.id = paste0(gene,"_",id))

add_subset = add %>% filter(age_lfsr < 0.05 | sexM_lfsr < 0.05)
int_subset = int %>% filter(gene.id %in% add_subset$gene.id)

################################################################################
# Size filtering with OPC and oligs separately 
################################################################################
####Use 30 donors and 300 cells as cutoff. (filter all basket cells in MB)
cell_count = read_csv("/scratch/klwatki5/R_analysis/bbi/nebula/infants_removed/cell_counts_updated.csv")
cell_count = cell_count %>% mutate(cell_type = case_when(
  (cell_type == "oligodendrocytes" & louvain %in% c(12,13)) ~ "OPC",
  TRUE ~ cell_type)) %>% 
  mutate(identity = paste0(cell_type,"_louvain",louvain,"_", region)) 

all_identity = cell_count$identity%>%unique()

e.filter.others = cell_count%>% filter(unique_animal_count<30 | cell_count < 300)%>%pull(identity)

e.filter.basket = cell_count%>%filter(cell_type=='basket cells')%>%
  filter(region=='MB'| unique_animal_count<30 | cell_count < 300)%>%pull(identity) 
e.filter = c(e.filter.basket,e.filter.others)

###############
#filtering
df_filter = int_subset %>% 
  filter(!id %in% e.filter) %>%         # filter out region x subtype lower than threshold
  filter(!is.na(p_sexM))                # filter out genes with NA in p_age or p_sex
df_filter = df_filter %>% filter(!(cell_type == "medium spiny neurons" & region == "HIP")) # 5,146,072

table(is.na(int_subset))

write_csv(df_filter, file = "/scratch/klwatki5/R_analysis/bbi/mashr/mash_new/interaction/interaction_nebula_results_filtered_for_input.csv")
