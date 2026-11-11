# Single-cell, multi-region profiling of the macaque brain across the lifespan

This repository contains analysis and figure-generation code associated with:

> Yang W, Watkins KL et al. **Single-cell, multi-region profiling of the macaque brain across the lifespan.** bioRxiv (2025). [https://doi.org/10.1101/2025.10.31.685880](https://www.biorxiv.org/content/10.1101/2025.10.31.685880v2)

We present a single-cell transcriptomic atlas spanning 11 brain regions and the rhesus macaque lifespan (0-21 yrs). The code in this repository covers data preprocessing, cell-type annotation, cell-abundance analysis, age- and sex- associated differential expression and pathway analysis, and comparison with Alzheimer's disease-associated transcriptional changes.

## Repository organization

```text
.
├── Figure_1/   # Data preprocessing, overview plots, PCA, and supplemental analyses
├── Figure_2/   # Cell annotation and subcluster/regulon analyses
├── Figure_3/   # Cell-abundance testing and pseudotime analyses
├── Figure_4/   # Age-associated differential expression and pathway analyses
├── Figure_5/   # Overlap with Alzheimer's disease-associated gene-expression changes
├── Figure_6/   # Sex-associated differential expression, age and sex overlap, and pathway analyses
└── Figure_7/   # Age-by-sex interaction analyses  
```

### Figure 1: atlas preprocessing and overview

- `Figure_1/data_preprocessing/`: region-specific Jupyter notebooks for preprocessing the single-cell datasets.
- `Figure_1/Figure1.ipynb`: analyses and plots used for the main atlas overview.
- `Figure_1/PCA_analysis.ipynb`: principal-component analyses of age-associated trajectories, including region-specific analyses.
- `Figure_1/supplemental_1.ipynb`: analyses and plots for Supplemental Figure 1, including infant-focused analyses.

### Figure 2: cell annotation and subcluster characterization

- `Figure_2/Run_tacco.py`: transfers cell-class and cell-subcluster labels from the BICCN macaque reference atlas to a query AnnData object using TACCO.
- `Figure_2/Run_tacco-Siletti.py`: transfers labels from the Siletti human brain atlas; non-neuronal reference cells are pooled across regions where needed.
- `Figure_2/scenic_subcluster.ipynb`: uses SCENIC/regulon activity to characterize novel cell subclusters and explore age-associated regulons.

### Figure 3: age-associated cell abundance and pseudotime

- `Figure_3/4_plot_brain.R`: maps `mashr` posterior cell-abundance coefficients across brain regions.
- `Figure_3/5_plotting_cell_abundance.R`: generates Figure 3 abundance plots, heatmaps, and correlation summaries.
- `Figure_3/7_clustering_hooke_spline.R`: clusters nonlinear age-associated cell-abundance trajectories for Figure 3I.
- `Figure_3/8_hooke_normalized_counts_subtype.R`: calculates normalized cell counts with `hooke` from cell-class and subtype models.
- `Figure_3/9_normalized_counts_heatmap.R`: visualizes normalized-count results as heatmaps.
- `Figure_3/pseudotime/`: microglial pseudotime analysis and age-associated gene visualization.

These scripts are primarily interactive analysis documents. Open them in RStudio and run them sequentially after updating their input and output paths.

### Figure 4: age-associated differential expression

The main documented workflow is:

1. `Figure_4/1_nebula_by_subtype-region.R` fits negative-binomial mixed models with `nebula` for each subcluster-by-region combination. The models test age- and sex-associated expression while controlling for library size; combinations represented by fewer than 30 animals are excluded.
2. `Figure_4/2_nebula_mashr_prepare_files.R` converts NEBULA log-fold changes and standard errors into matrices for `mashr`.
3. `Figure_4/3_nebula_mashr_by-celltype_fullset.R` applies multivariate adaptive shrinkage by cell class and produces posterior means, posterior standard deviations, and local false-sign rates (LFSR).
4. `Figure_4/4_mashr_summary_clean.R` combines NEBULA and `mashr` results, filters significant gene-by-condition results, and performs downstream summaries.
5. `Figure_4/5_GSEA_regionxsubtype_v2.R` performs gene-set enrichment analysis using GO, KEGG, Hallmark, SynGO, and SASP-related gene sets.

Additional scripts examine region-specific and cell-type-specific genes, summarize enriched ontology terms, identify gene modules, compare subtype patterns, and calculate pathway gene scores.

For a next flow pipeline of age-associated differential expression: https://github.com/vinettey/macaque-aging-de

### Figure 5: comparison with Alzheimer's disease

- `Figure_5/AD_overlap_clean.R`: tests overlap between macaque age-associated differentially expressed genes and AD-associated genes from three independent age-matched case-control single-cell studies.
- `Figure_5/plot_AD_overlap.ipynb`: generates downstream visualizations of shared age- and AD-associated expression patterns.

## Figure 6: Sex-associated differential expression, overlap of age and sex on expression

- `Figure_6/mashr_for_sex_degs/`

The documented workflow is:

1. `1_prep_files.R` converts NEBULA log-fold changes and standard error output from the sex covariate in `Figure_4/1_nebula_by_subtype-region.R` into matrices for `mashr`. Genes on autosomes and X/Y chromosomes are analyzed separately to avoid artificial inflation. 
2. `2_run_mashr_sex.R` applies multivariate adaptive shrinkage by cell class and produces posterior means, posterior standard deviations, and local false-sign rates (LFSR). 
3. `3_create_results_dataframes.R` combines NEBULA and `mashr` results, filters significant gene-by-condition results, and performs downstream summaries.

- `Fig_6_plots.R`: performs statistical analyses and generates plots shown in Figure 6.
- `Fig_6_ORA_analysis.R`: performs statistical over-representation analysis of age + sex DEGs among gene ontology datasets from `msigdbr` and SynGO. 
- `Fig_6_AD_DEG.R`: tests overlap between macaque a+sDEGs and AD-associated genes from Mathys et al., 2024.

These scripts are primarily interactive analysis documents. Open them in RStudio and run them sequentially after updating their input and output paths.

## Figure 7: age-by-sex interactions on expression 

- `Figure_7/mashr_for_age_sex_interaction/`

The documented workflow is: 

1. `1_post_hoc_filtering.R` filters NEBULA interaction model output for genes that showed either a significant age or sex interaction in the additive model. 
2. `2_prepare_interaction_input_files.R` converts NEBULA log-fold changes and standard error output from the age, sex, interaction, and intercept terms into matrices for `mashr`. Genes on autosomes and X/Y chromosomes are analyzed separately to avoid artificial inflation.
3. `3_run_mashr_interaction.R` applies multivariate adaptive shrinkage by cell class and produces posterior means, posterior standard deviations, and local false-sign rates (LFSR) for each term. 
4. `4_create_results_dataframes.R` combines NEBULA and `mashr` results, filters significant gene-by-condition results, and performs downstream summaries.

- `Fig_7_plots.R`: performs statistical analyses and generates plots shown in Figure 7.

These scripts are primarily interactive analysis documents. Open them in RStudio and run them sequentially after updating their input and output paths.

## Software requirements

The repository combines Python/Jupyter and R workflows. Core software used across the scripts includes:

- Python: `anndata`, `scanpy`, `tacco`, `pandas`, `numpy`, and Jupyter.
- R: `tidyverse`, `nebula`, `mashr`, `hooke`, `monocle3`, `Matrix`, `ggplot2`, `pheatmap`, `clusterProfiler`, `msigdbr`, and `org.Mmu.eg.db`.
- Selected analyses also use SCENIC/pySCENIC, SynGO annotations, MSigDB gene sets, and the `macaquesimplebrain` package.

Python package versions are documented in the corresponding jupyter notebooks. R pacakge versions are to be updated.

## Brain-region abbreviations

The analyses include 11 regions: anterior cingulate cortex (ACC), caudate nucleus (CN), dorsolateral prefrontal cortex (dlPFC), entorhinal cortex (EC), hippocampus (HIP), inferior posterior parietal cortex (IPP), lateral cerebellum (lCb), primary motor cortex (M1), midbrain (MB), mediodorsal thalamic nucleus (mdTN), and nucleus accumbens (NAc).

## Citation

If you use this code, please cite the associated preprint:

```bibtex
@article{yang2025macaqueaging,
  title   = {Single-cell, multi-region profiling of the macaque brain across the lifespan},
  author  = {Yang, Wei and Watkins, Kelsi L. and DeCasien, Alex R. and others},
  journal = {bioRxiv},
  year    = {2025},
  doi     = {10.1101/2025.10.31.685880},
  url     = {https://www.biorxiv.org/content/10.1101/2025.10.31.685880v2}
}
```

## License and contact

No license file is currently included. Please contact the repository authors before reusing or redistributing the code beyond the permissions provided by applicable law. Questions about individual analyses should be directed to the corresponding authors listed in the preprint.
