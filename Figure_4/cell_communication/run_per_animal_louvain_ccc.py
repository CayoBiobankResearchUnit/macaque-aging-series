#!/usr/bin/env python3
"""
Run LIANA per individual animal with Louvain cluster resolution - WITHIN ONE REGION
v2: OPC and Oligo split by per-cell cell_class_annotation instead of louvain majority
Uses CellPhoneDB resource with neurotransmitter signaling, expr_prop=0.05, FDR on cellphone_pvals
"""
import scanpy as sc
import pandas as pd
import numpy as np
import liana as li
import pandas as _pd
from statsmodels.stats.multitest import multipletests
import os
import h5py
import sys

# Load and validate CellPhoneDB resource
print("Loading CellPhoneDB resource...")
cpdb_resource = _pd.read_csv("/scratch/easmit31/cell_cell/cellphonedb_interactions_liana_format.csv")
assert len(cpdb_resource) == 2909, f"ERROR: expected 2716 LR pairs, got {len(cpdb_resource)}"
neuro_check = cpdb_resource[cpdb_resource.apply(
    lambda r: r.astype(str).str.contains('GABR|GRIN|DRD|HTR').any(), axis=1)]
assert len(neuro_check) > 50, f"ERROR: neurotransmitter pairs missing, only {len(neuro_check)} found"
print(f"✓ CellPhoneDB resource loaded: {len(cpdb_resource)} LR pairs, {len(neuro_check)} neurotransmitter-related")

base_path = "/scratch/nsnyderm/u01/intermediate_files/cell-class_h5ad_update/"
cell_type_files = {
    'Astrocyte': 'Res1_astrocytes_update.h5ad',
    'Basket': 'Res1_basket-cells_update.h5ad',
    'Cerebellar': 'Res1_cerebellar-neurons_subset.h5ad',
    'Ependymal': 'Res1_ependymal-cells_new.h5ad',
    'GABA': 'Res1_GABAergic-neurons_subset.h5ad',
    'Glutamatergic': 'Res1_glutamatergic-neurons_update.h5ad',
    'MSN': 'Res1_medium-spiny-neurons_subset.h5ad',
    'Microglia': 'Res1_microglia_new.h5ad',
    'Midbrain': 'Res1_midbrain-neurons_update.h5ad',
    'OPC_Oligo': 'Res1_opc-olig_subset.h5ad',
    'Vascular': 'Res1_vascular-cells_subset.h5ad'
}

OPC_OLIGO_LABEL_MAP = {
    'oligodendrocyte precursor cells': 'OPC',
    'oligodendrocytes': 'Oligo'
}

def load_h5ad_animal_with_louvain(h5ad_path, animal_id, region, cell_type):
    with h5py.File(h5ad_path, 'r') as f:
        animal_id_group = f['obs']['animal_id']
        animal_cats = [c.decode('utf-8') if isinstance(c, bytes) else str(c)
                      for c in animal_id_group['categories'][:]]
        animal_codes = animal_id_group['codes'][:]
        animal_ids = np.array([animal_cats[c] for c in animal_codes])

        region_group = f['obs']['region']
        region_cats = [c.decode('utf-8') if isinstance(c, bytes) else str(c)
                      for c in region_group['categories'][:]]
        region_codes = region_group['codes'][:]
        regions = np.array([region_cats[c] for c in region_codes])

        mask = (animal_ids == animal_id) & (regions == region)
        n_cells = mask.sum()

        if n_cells == 0:
            return None

        louvain_group = f['obs']['louvain']
        if 'categories' in louvain_group:
            louvain_cats = [c.decode('utf-8') if isinstance(c, bytes) else str(c)
                           for c in louvain_group['categories'][:]]
            louvain_codes = louvain_group['codes'][:][mask]
            louvains = np.array([louvain_cats[c] for c in louvain_codes])
        else:
            louvains = np.array(['0'] * n_cells)

        age = f['obs']['age'][:][mask][0]

        sex_group = f['obs']['sex']
        sex_cats = [c.decode('utf-8') if isinstance(c, bytes) else str(c)
                   for c in sex_group['categories'][:]]
        sex_codes = sex_group['codes'][:][mask]
        sex = sex_cats[sex_codes[0]]

        cell_class_labels = None
        if cell_type == 'OPC_Oligo' and 'cell_class_annotation' in f['obs']:
            cca_group = f['obs']['cell_class_annotation']
            if 'categories' in cca_group:
                cca_cats = [c.decode('utf-8') if isinstance(c, bytes) else str(c)
                           for c in cca_group['categories'][:]]
                cca_codes = cca_group['codes'][:][mask]
                cell_class_labels = np.array([cca_cats[c] for c in cca_codes])

    adata = sc.read_h5ad(h5ad_path)
    adata = adata[mask].copy()

    if 'external_gene_name' in adata.var.columns:
        gene_names = adata.var['external_gene_name'].values
        if hasattr(gene_names, 'categories'):
            gene_names = gene_names.astype(str)
        adata.var_names = gene_names
        adata.var_names = adata.var_names.astype(str)
        adata.var_names_make_unique()

    if cell_type == 'OPC_Oligo' and cell_class_labels is not None:
        ct_labels = []
        for cca, lou in zip(cell_class_labels, louvains):
            prefix = OPC_OLIGO_LABEL_MAP.get(cca, 'OPC_Oligo')
            ct_labels.append(f"{prefix}_{lou}")
        adata.obs['celltype_louvain'] = ct_labels
    else:
        adata.obs['celltype_louvain'] = [f"{cell_type}_{lou}" for lou in louvains]

    adata.obs['age'] = age
    adata.obs['sex'] = sex
    adata.obs['animal_id'] = animal_id
    adata.obs['region'] = region

    return adata, age, sex


def run_animal_analysis(animal_id, region, output_dir=None):
    print(f"\n{'='*70}")
    print(f"Processing: {animal_id} - {region}")
    print(f"{'='*70}")

    if output_dir is None:
        output_dir = "/scratch/easmit31/cell_cell/results/per_animal_louvain_corrected"
    os.makedirs(output_dir, exist_ok=True)

    print(f"Output dir: {output_dir}")

    print(f"[Step 1] Loading cell types...")
    adatas = []
    age = None
    sex = None

    for cell_type, filename in cell_type_files.items():
        h5ad_path = os.path.join(base_path, filename)
        result = load_h5ad_animal_with_louvain(h5ad_path, animal_id, region, cell_type)
        if result is not None:
            adata_ct, age, sex = result
            print(f"  ✓ {cell_type}: {adata_ct.n_obs:,} cells")
            if cell_type == 'OPC_Oligo':
                print(f"    celltype_louvain groups: {adata_ct.obs['celltype_louvain'].value_counts().to_dict()}")
            adatas.append(adata_ct)

    if len(adatas) == 0:
        print(f"⚠️  No cells found for {animal_id} in {region}!")
        return

    print(f"\n[Step 2] Combining {len(adatas)} cell types...")
    adata_combined = sc.concat(adatas, join='outer')
    print(f"  Total cells: {adata_combined.n_obs:,}")
    print(f"  Age: {age}, Sex: {sex}")

    cluster_counts = adata_combined.obs['celltype_louvain'].value_counts()
    print(f"\n  Unique celltype-louvain groups: {len(cluster_counts)}")
    print(f"  Groups with >10 cells: {(cluster_counts > 10).sum()}")

    print(f"\n[Step 3] Normalizing and log-transforming...")
    sc.pp.normalize_total(adata_combined, target_sum=1e4)
    sc.pp.log1p(adata_combined)

    print(f"\n[Step 4] Running LIANA cell-cell communication...")
    li.mt.rank_aggregate(
        adata_combined,
        groupby='celltype_louvain',
        resource=cpdb_resource,
        expr_prop=0.05,
        n_perms=100,
        verbose=True,
        use_raw=False
    )

    print(f"\n[Step 5] Extracting and processing results...")
    results_df = adata_combined.uns['liana_res'].copy()

    results_df['animal_id'] = animal_id
    results_df['age'] = age
    results_df['sex'] = sex
    results_df['region'] = region

    _, results_df['qval'], _, _ = multipletests(results_df['cellphone_pvals'], method='fdr_bh')

    output_file = os.path.join(output_dir, f"{animal_id}_{region}_louvain_results.csv")
    results_df.to_csv(output_file, index=False)

    print(f"\n✓ Results saved to: {output_file}")
    print(f"  Total interactions tested: {len(results_df):,}")
    print(f"  Significant (qval < 0.05): {(results_df['qval'] < 0.05).sum():,}")

    return results_df


if __name__ == "__main__":
    if len(sys.argv) not in (3, 4):
        print("Usage: python run_per_animal_louvain_corrected.py <ANIMAL_ID> <REGION> [OUTPUT_DIR]")
        print("Example: python run_per_animal_louvain_corrected.py 0B9 HIP")
        print("Example: python run_per_animal_louvain_corrected.py 0B9 dlPFC /scratch/easmit31/cell_cell/results/per_animal_louvain_dlpfc_attempt2")
        sys.exit(1)

    animal_id = sys.argv[1]
    region = sys.argv[2]
    output_dir = sys.argv[3] if len(sys.argv) == 4 else None
    run_animal_analysis(animal_id, region, output_dir)
    print(f"\n{'='*70}")
    print(f"✓ Completed: {animal_id} - {region}")
    print(f"{'='*70}")
