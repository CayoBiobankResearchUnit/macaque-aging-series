#Script to annotate cell class and cell subcluster with tacco
#author: wei yang
import os
import tacco as tc
import sys
import scanpy as sc
import numpy as np
import pandas as pd
import random

#Params
region = 'mdtn'
ref_region = ['mdTN']
ref_dir = '/scratch/wyang114/processed_data/reference/anndata_rna_biccn-all_filtered.h5ad'
#query_dir = f'/scratch/wyang114/processed_data/anndata_rna_{region.lower()}_processed.h5ad'
query_dir='/scratch/wyang114/processed_data/filtered_data_Feb/filtered_data_wSiletti/mdtn_processed_w_siletti.h5ad'

#script
print('processing '+region)
print('loading ref...')
reference = sc.read_h5ad(ref_dir) 
select_genes = reference.var.index[np.where(~reference.var['chromosome'].isin(['X','Y']))]

#filter cells
df = reference.obs[reference.obs['region'].isin(ref_region)]
temp = df['cell_class'].value_counts(normalize=True)
select_cell_class = temp[temp>0].index
select_cells = reference.obs[reference.obs['cell_class'].isin(select_cell_class)].index
reference_subset = reference[select_cells,select_genes].copy()
reference_subset = reference_subset[reference_subset.obs['region'].isin(ref_region)]
print(reference_subset.shape)


print('loading query...')
query_data = sc.read_h5ad(query_dir) #enter query data directory
print(query_data.shape)

#Cell class
annotation_key = 'cell_class' 
result_key = 'tacco_cell_class' 
level = 'cell_class' 
print('Running tacco on cell class..')
query_data.X = query_data.raw.X
#query_data.X = query_data.layers['raw'] #(optional)run this line if the count data is stored in raw
tc.tl.annotate(query_data, reference_subset, annotation_key= annotation_key,result_key= result_key,counts_location='X',normalize_to='reference') #specify annotation_key to the relevant column in reference
print('Outputing..check if the dimension match(1)')
sum(query_data.obsm[result_key].index == query_data.obs.index)/query_data.shape[0] 
query_data.obs[result_key]=query_data.obsm[result_key].idxmax(axis=1) #extract tacco assignment information
query_data.obs[result_key+'_score']=query_data.obsm[result_key].max(axis=1) #extract tacco score information
sc.pl.umap(query_data,color= result_key, save='_%s_tacco_%s_umap_adata_filtered.png'%(region,level))
sc.pl.umap(query_data,color= result_key+'_score', save='_%s_tacco_%s_score_umap_adata_filtered.png'%(region,level))
df_ref = reference_subset.obs[annotation_key].value_counts(normalize=True)
df_query = query_data.obs[result_key].value_counts(normalize=True)
df = pd.DataFrame({'reference':df_ref,f'query_tacco_{region}_ref':df_query})
df.fillna(0, inplace=True)
#df.to_csv(f'/scratch/wyang114/outputs/tacco_cell_fraction/{ref_region}_{level}_ref_query_tacco.csv',float_format='%.6f')

# #Cell subcluster
# annotation_key = 'cell_subcluster' 
# result_key = 'tacco_cell_subcluster' 
# level = 'cell_subcluster' 
# print('Running tacco on cell class..')
# query_data.X = query_data.layers['raw'] #(optional)run this line if the count data is stored in raw
# tc.tl.annotate(query_data, reference_subset, annotation_key= annotation_key,result_key= result_key,counts_location='X',normalize_to='reference') #specify annotation_key to the relevant column in reference
# print('Outputing..check if the dimension match(1)')
# sum(query_data.obsm[result_key].index == query_data.obs.index)/query_data.shape[0] 
# query_data.obs[result_key]=query_data.obsm[result_key].idxmax(axis=1) #extract tacco assignment information
# query_data.obs[result_key+'_score']=query_data.obsm[result_key].max(axis=1) #extract tacco score information
# sc.pl.umap(query_data,color= result_key, save='_%s_tacco_%s_umap_adata.png'%(region,level))
# sc.pl.umap(query_data,color= result_key+'_score', save='_%s_tacco_%s_score_umap_adata.png'%(region,level))
# df_ref = reference_subset.obs[annotation_key].value_counts(normalize=True)
# df_query = query_data.obs[result_key].value_counts(normalize=True)
# df = pd.DataFrame({'reference':df_ref,f'query_tacco_{region}_ref':df_query})
# df.fillna(0, inplace=True)
# df.to_csv(f'/scratch/wyang114/outputs/tacco_cell_fraction/{ref_region}_{level}_ref_query_tacco.csv',float_format='%.6f')


#print('Saving files..')
#query_data.write_h5ad('/scratch/wyang114/processed_data/anndata_rna_%s_processed.h5ad'%(region.lower())) #save anndata object with tacco information as h5ad