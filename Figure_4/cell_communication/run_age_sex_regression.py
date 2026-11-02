#!/usr/bin/env python3
"""
Run age + sex regression for cell communication.
Model: lr_means ~ age + sex (OLS) for each unique interaction.
Reads from pre-built no-threshold lr_matrix CSV.
Filters to interactions present in >= min_animals animals.
FDR correction (Benjamini-Hochberg) applied separately for age and sex.
Usage:
    python run_age_sex_regression.py --region HIP
    python run_age_sex_regression.py --region HIP --min_animals 10
    python run_age_sex_regression.py --region ACC --exclude_animals 8H2
    python run_age_sex_regression.py --region dlPFC --matrix_file /path/to/matrix.csv --output_suffix _attempt2
"""

import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

parser = argparse.ArgumentParser()
parser.add_argument('--region', required=True)
parser.add_argument('--min_animals', type=int, default=10)
parser.add_argument('--min_age', type=float, default=1.0)
parser.add_argument('--exclude_animals', type=str, default=None,
                    help='Comma-separated animal IDs to exclude (e.g. 8H2,0B9)')
parser.add_argument('--matrix_file', type=str, default=None,
                    help='Override default matrix file path')
parser.add_argument('--output_suffix', type=str, default='',
                    help='Suffix to append to output filename (e.g. _attempt2)')
args = parser.parse_args()

region = args.region
min_animals = args.min_animals
min_age = args.min_age

print("="*80)
print(f"AGE + SEX REGRESSION: {region} [whole region]")
print("="*80)
print(f"  Min animals: {min_animals}")
print(f"  Min age:     {min_age}")
if args.exclude_animals:
    print(f"  Excluding:   {args.exclude_animals}")
if args.output_suffix:
    print(f"  Output suffix: {args.output_suffix}")

if args.matrix_file:
    matrix_file = args.matrix_file
else:
    matrix_file = f'/scratch/easmit31/cell_cell/results/lr_matrices_corrected/{region}_nothresh_minage{str(min_age).replace(".","p")}_matrix.csv'
print(f"\nLoading matrix: {matrix_file}")

mat = pd.read_csv(matrix_file, index_col=0, low_memory=False)
age_row = mat.loc['age'].astype(float)
sex_row = mat.loc['sex']
mat = mat.drop(index=['age', 'sex'])

if args.exclude_animals:
    exclude = [a.strip() for a in args.exclude_animals.split(',')]
    before = mat.shape[1]
    mat = mat.drop(columns=[a for a in exclude if a in mat.columns], errors='ignore')
    age_row = age_row.drop([a for a in exclude if a in age_row.index], errors='ignore')
    sex_row = sex_row.drop([a for a in exclude if a in sex_row.index], errors='ignore')
    print(f"Excluded {before - mat.shape[1]} animals: {exclude} ({before} -> {mat.shape[1]} remaining)")

print(f"Matrix: {mat.shape[0]:,} interactions x {mat.shape[1]} animals")

print(f"\nFiltering to interactions in >={min_animals} animals...")
n_valid = mat.notna().sum(axis=1)
mat = mat[n_valid >= min_animals]
print(f"Interactions to model: {len(mat):,}")

print("\nRunning OLS: lr_means ~ age + sex...")
results = []

for i, (interaction, row) in enumerate(mat.iterrows()):
    if i % 10000 == 0:
        print(f"  Model {i+1:,}/{len(mat):,}...")

    valid = row.dropna()
    ages = age_row[valid.index].values.astype(float)
    sexes = sex_row[valid.index].values
    lr_vals = valid.values.astype(float)

    sex_binary = (sexes == 'M').astype(int)
    if len(np.unique(sex_binary)) < 2:
        continue

    try:
        X = sm.add_constant(np.column_stack([ages, sex_binary]))
        model = sm.OLS(lr_vals, X).fit()
        results.append({
            'interaction': interaction,
            'n_animals': len(valid),
            'mean_lr_means': lr_vals.mean(),
            'age_coef': model.params[1],
            'age_stderr': model.bse[1],
            'age_pval': model.pvalues[1],
            'sex_coef': model.params[2],
            'sex_stderr': model.bse[2],
            'sex_pval': model.pvalues[2],
            'r_squared': model.rsquared
        })
    except Exception:
        continue

print(f"\nModeled {len(results):,} interactions")
results_df = pd.DataFrame(results)

_, age_qvals, _, _ = multipletests(results_df['age_pval'], method='fdr_bh')
_, sex_qvals, _, _ = multipletests(results_df['sex_pval'], method='fdr_bh')
results_df['age_qval'] = age_qvals
results_df['sex_qval'] = sex_qvals

output_dir = Path(f'/scratch/easmit31/cell_cell/results/within_region_analysis_corrected/regression_{region}')
output_dir.mkdir(parents=True, exist_ok=True)

excl_tag = '_' + '_'.join(args.exclude_animals.split(',')) + '_excluded' if args.exclude_animals else ''
output_file = output_dir / f'whole_{region.lower()}_age_sex_regression{excl_tag}{args.output_suffix}.csv'
results_df.to_csv(output_file, index=False)
print(f"\n✓ Saved: {output_file}")

print(f"\nTotal modeled: {len(results_df):,}")
print(f"Age  p<0.05: {(results_df['age_pval'] < 0.05).sum():,}  q<0.05: {(results_df['age_qval'] < 0.05).sum():,}")
print(f"Sex  p<0.05: {(results_df['sex_pval'] < 0.05).sum():,}  q<0.05: {(results_df['sex_qval'] < 0.05).sum():,}")

if (results_df['age_qval'] < 0.05).sum() > 0:
    print(f"\nTop 10 age associations (by p-value):")
    print(results_df[results_df['age_qval'] < 0.05].nsmallest(10, 'age_pval')[
        ['interaction', 'age_coef', 'age_pval', 'age_qval', 'n_animals']].to_string(index=False))
