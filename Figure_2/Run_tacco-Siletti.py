#Script to annotate cell class and cell subcluster with tacco
#author: wei yang
import os
import tacco as tc
import sys
import scanpy as sc
import numpy as np
import pandas as pd

#Params
region = 'MB'
ref_dir = '/scratch/klchiou/bican/Siletti_subsampled_AIT21_scos.h5ad'
query_dir = f'/scratch/wyang114/processed_data/filtered_data_Feb/anndata_rna_{region.lower()}_filtered.h5ad'

#script
print('processing '+region)
print('loading ref...')
reference = sc.read_h5ad(ref_dir) 
select_genes = reference.var.index


#filter cells
df = reference.obs
ref_region = ['Human CaB']
non_neurons = ['Astrocyte','Bergmann glia','Committed oligodendrocyte precursor','Ependymal','Fibroblast','Microglia','Oligodendrocyte',
              'Oligodendrocyte precursor','Vascular',
              'Miscellaneous'
              ]
df_subset = df[(df['roi'].isin(ref_region)) | (df['supercluster_term'].isin(non_neurons))]
#df_subset = df #use full dataset
select_cells = df_subset.index
reference_subset = reference[select_cells,:].copy()
print(reference_subset.shape)


print('loading query...')
query_data = sc.read_h5ad(query_dir) #enter query data directory
#query_data = query_data[:, query_data.var.highly_variable]
print(query_data.raw.shape)

#Cell class
annotation_key = 'AIT21_all.class' 
result_key = 'siletti_tacco_cell_class' 
level = 'supercluster_term' 
print('Running tacco on cell class..')
#query_data.X = query_data.layers['raw'] #(optional)run this line if the count data is stored in raw
tc.tl.annotate(query_data, reference_subset, annotation_key= annotation_key,result_key= result_key,counts_location='X',normalize_to='data') #specify annotation_key to the relevant column in reference
print('Outputing..check if the dimension match(1)')
sum(query_data.obsm[result_key].index == query_data.obs.index)/query_data.shape[0] 
query_data.obs[result_key]=query_data.obsm[result_key].idxmax(axis=1) #extract tacco assignment information
query_data.obs[result_key+'_score']=query_data.obsm[result_key].max(axis=1) #extract tacco score information
sc.pl.umap(query_data,color= result_key, save='_%s_tacco_%s_umap_adata.png'%(region,level))
sc.pl.umap(query_data,color= result_key+'_score', save='_%s_tacco_%s_score_umap_adata.png'%(region,level))
df_ref = reference_subset.obs[annotation_key].value_counts(normalize=True)
df_query = query_data.obs[result_key].value_counts(normalize=True)
df = pd.DataFrame({'reference':df_ref,f'query_tacco_{region}_ref':df_query})
df.fillna(0, inplace=True)
#df.to_csv(f'/scratch/wyang114/siletti/{region}_{level}_ref_query_tacco.csv',float_format='%.6f')

# #Cell subcluster
annotation_key = 'AIT21_all.subclass' 
result_key = 'siletti_tacco_cell_subcluster' 
level = 'AIT21_all.subclass' 
print('Running tacco on cell subclass..')
query_data.X = query_data.layers['raw'] #(optional)run this line if the count data is stored in raw
#query_data.X = query_data.raw.X
tc.tl.annotate(query_data, reference_subset, annotation_key= annotation_key,result_key= result_key,counts_location='X',normalize_to='data') #specify annotation_key to the relevant column in reference
print('Outputing..check if the dimension match(1)')
sum(query_data.obsm[result_key].index == query_data.obs.index)/query_data.shape[0] 
query_data.obs[result_key]=query_data.obsm[result_key].idxmax(axis=1) #extract tacco assignment information
query_data.obs[result_key+'_score']=query_data.obsm[result_key].max(axis=1) #extract tacco score information
sc.pl.umap(query_data,color= result_key, save='_%s_tacco_%s_umap_adata.png'%(region,level))
sc.pl.umap(query_data,color= result_key+'_score', save='_%s_tacco_%s_score_umap_adata.png'%(region,level))
df_ref = reference_subset.obs[annotation_key].value_counts(normalize=True)
df_query = query_data.obs[result_key].value_counts(normalize=True)
df = pd.DataFrame({'reference':df_ref,f'query_tacco_{region}_ref':df_query})
df.fillna(0, inplace=True)
df.to_csv(f'/scratch/wyang114/siletti/{region}_{level}_ref_query_tacco.csv',float_format='%.6f')


print('Saving files..')
query_data.write_h5ad(f'/scratch/wyang114/processed_data/filtered_data_Feb/anndata_rna_{region.lower()}_filtered_wfrank.h5ad') #save anndata object with tacco information as h5ad