#!/usr/bin/env python3
"""
plot_chord_category.py

Circular "chord" plots of age-associated cell-cell communication (CCC) for one
functional LR category, one plot per brain region (or all 11 regions in a grid).

What is drawn
-------------
- Nodes   = broad cell classes (AST, INH, EXC, MGL, OPC, OLIG, VASC, BC, CER,
            EPEN, MSN, MBN), placed alphabetically by class name, evenly spaced
            counter-clockwise starting at 3 o'clock. Louvain subtypes are
            collapsed to their parent class (trailing "_<n>" stripped).
- Arrows  = sender -> receiver class pairs with at least one LR pair in the
            chosen category that is significant for age (age_qval < 0.05).
- Color   = mean age coefficient (age_coef) across those significant
            interactions; RdBu_r centered at 0 (red = stronger with age,
            blue = weaker with age).
- Width   = number of significant interactions for that sender -> receiver
            pair, scaled 1-6 relative to the max within each plot.
- Self-loops (a class signalling to itself) are drawn as arcs outside the
  node ring.

Colorbar modes (--colorbar)
---------------------------
   1         one color scale shared across all loaded regions
             (grid: single colorbar on the right)
   0         shared color scale, no colorbar drawn
   omitted   color scale set separately per region, colorbar per plot

Inputs (read only)
------------------
1. Regression results, one file per region:
     /scratch/easmit31/cell_cell/results/within_region_analysis_corrected/
       regression_results/regression_{region}/whole_{region}_age_sex_regression.csv
   Columns used:
     interaction   "sender|receiver|ligand|receptor"
     age_qval      filtered to < 0.05
     age_coef      age effect, plotted as color
   Produced by run_age_sex_regression.py (OLS, age + sex), which runs on the
   per-animal LIANA / CellPhoneDB v5 output from run_per_animal_louvain_ccc.py
   (per_animal_louvain_corrected/).

2. LR functional annotations:
     /scratch/easmit31/cell_cell/cpdb_lr_annotations.csv
   Columns used: ligand, receptor, broad_category.
   This is the corrected annotation file (2,909 LR pairs), renamed from
   cpdb_lr_annotations_corrected_v2.csv; older versions are archived.
   Originally built by build_cpdb_annotations.py from the CellPhoneDB v5
   interaction table. --category must match a broad_category value exactly;
   if not, the script prints the available categories and exits.

No other project modules are imported (does not use plot_style.py).

Outputs
-------
  {BASE_DIR}/chord_plots_{Category}/{region}_chord_{Category}.png
      one per region (default mode)
  {BASE_DIR}/chord_plots_{Category}/grid_chord_{Category}_11regions{suffix}.png
      --grid mode; suffix = _colorbar1, _colorbar0 or _percbar
  Existing files with the same name are overwritten.

Usage
-----
    python plot_chord_category.py --category "Synaptic adhesion" --colorbar 1
    python plot_chord_category.py --category "Glutamate signaling" --grid --colorbar 1
    python plot_chord_category.py --category "GABA signaling" --region dlPFC --colorbar 1
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from matplotlib.colors import TwoSlopeNorm
from matplotlib.cm import ScalarMappable
from pathlib import Path
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--category', type=str, required=True,
                    help='Functional category to plot (must match broad_category in annotations)')
parser.add_argument('--colorbar', type=int, default=-1,
                    help='1=shared colorbar, 0=no colorbar, omit=per-plot colorbar')
parser.add_argument('--grid', action='store_true',
                    help='3x4 grid of all 11 regions')
parser.add_argument('--region', type=str, default=None,
                    help='Single region to plot (e.g. dlPFC)')
args = parser.parse_args()

REGIONS      = ["ACC","CN","dlPFC","EC","HIP","IPP","lCb","M1","MB","mdTN","NAc"]
GRID_REGIONS = ["ACC","CN","dlPFC","EC","HIP","IPP","lCb","M1","MB","mdTN","NAc"]
BASE_DIR     = Path("/scratch/easmit31/cell_cell/results/within_region_analysis_corrected/regression_results")
ANN_FILE     = Path("/scratch/easmit31/cell_cell/cpdb_lr_annotations.csv")

# output dir named after category (spaces replaced with underscores)
cat_slug = args.category.replace(' ', '_').replace('/', '_')
OUT_DIR  = BASE_DIR / f"chord_plots_{cat_slug}"
OUT_DIR.mkdir(exist_ok=True)

CELLTYPE_COLORS = {
    'Astrocyte':  '#d53a50', 'GABA':       '#4d6aa8',
    'Glut':       '#3d9995', 'Microglia':  '#bb471c',
    'OPC':        '#fbae64', 'Oligo':      '#fb8d64',
    'Vascular':   '#d1538f', 'Basket':     '#4e379a',
    'Cerebellar': '#44839f', 'Ependymal':  '#fc6238',
    'MSN':        '#5e3c99', 'Midbrain':   '#65bfb0',
}

CELLTYPE_LABELS = {
    'Astrocyte':  'AST',  'GABA':       'INH',
    'Glut':       'EXC',  'Microglia':  'MGL',
    'OPC':        'OPC',  'Oligo':      'OLIG',
    'Vascular':   'VASC', 'Basket':     'BC',
    'Cerebellar': 'CER',  'Ependymal':  'EPEN',
    'MSN':        'MSN',  'Midbrain':   'MBN',
}

def get_class(s):
    return next((k for k in CELLTYPE_COLORS if s.startswith(k)), 'Unknown')

def get_color(c):
    return CELLTYPE_COLORS.get(c, '#aaaaaa')

ann = pd.read_csv(ANN_FILE)
ann['lr_pair'] = ann['ligand'] + '|' + ann['receptor']
cat_pairs = set(ann[ann['broad_category'] == args.category]['lr_pair'])
print(f"Category: {args.category}")
print(f"LR pairs: {len(cat_pairs)}")

if len(cat_pairs) == 0:
    print("ERROR: No LR pairs found for this category. Check category name.")
    print("Available categories:", sorted(ann['broad_category'].unique()))
    exit(1)

regions_to_load = [args.region] if args.region else REGIONS

all_edges = {}
for region in regions_to_load:
    csv = BASE_DIR / f"regression_{region}" / f"whole_{region.lower()}_age_sex_regression.csv"
    if not csv.exists():
        print(f"  ⚠ Not found: {csv}")
        continue
    print(f"Loading {region}...")
    chunks = []
    for chunk in pd.read_csv(csv, chunksize=500000):
        parts = chunk['interaction'].str.split('|', expand=True)
        chunk['sender_class']   = parts[0].str.replace(r'_\d+$','',regex=True).apply(get_class)
        chunk['receiver_class'] = parts[1].str.replace(r'_\d+$','',regex=True).apply(get_class)
        chunk['lr_pair']        = parts[2] + '|' + parts[3]
        sub = chunk[(chunk['age_qval'] < 0.05) & (chunk['lr_pair'].isin(cat_pairs))]
        if len(sub): chunks.append(sub)
    if not chunks:
        print(f"  No significant interactions for {region}")
        continue
    sig = pd.concat(chunks)
    edges = (sig.groupby(['sender_class','receiver_class'])
               .agg(n_sig=('age_coef','count'), mean_coef=('age_coef','mean'))
               .reset_index())
    all_edges[region] = edges

if not all_edges:
    print("No data loaded for any region.")
    exit(1)

all_coefs = pd.concat(all_edges.values())['mean_coef']
global_vmax = max(all_coefs.abs().max(), 0.001)
global_norm = TwoSlopeNorm(vmin=-global_vmax, vcenter=0, vmax=global_vmax)
cmap = plt.cm.RdBu_r

def draw_on_ax(ax, region, edges, use_norm, show_cbar=False):
    celltypes = sorted(set(edges['sender_class']) | set(edges['receiver_class']))
    n = len(celltypes)
    angles = {ct: 2*np.pi*i/n for i, ct in enumerate(celltypes)}
    pos    = {ct: (np.cos(angles[ct]), np.sin(angles[ct])) for ct in celltypes}

    NODE_RADIUS = 0.12
    max_sig = edges['n_sig'].max()
    min_lw, max_lw = 1.0, 6.0

    for _, row in edges.iterrows():
        src, tgt = row['sender_class'], row['receiver_class']
        x0, y0 = pos[src]; x1, y1 = pos[tgt]
        lw = min_lw + (max_lw - min_lw) * (row['n_sig'] / max_sig)
        color = cmap(use_norm(row['mean_coef']))

        if src == tgt:
            angle = angles[src]
            ax.add_patch(mpatches.Arc((1.35*np.cos(angle), 1.35*np.sin(angle)),
                                      0.25, 0.25, angle=0, theta1=0, theta2=300,
                                      color=color, lw=lw))
        else:
            dx, dy = x1-x0, y1-y0
            dist = np.sqrt(dx**2+dy**2)
            ux, uy = dx/dist, dy/dist
            sx, sy = x0+ux*NODE_RADIUS, y0+uy*NODE_RADIUS
            ex, ey = x1-ux*NODE_RADIUS, y1-uy*NODE_RADIUS
            mid_x = (sx+ex)/2 - 0.25*uy
            mid_y = (sy+ey)/2 + 0.25*ux
            t = np.linspace(0, 1, 100)
            bx = (1-t)**2*sx + 2*(1-t)*t*mid_x + t**2*ex
            by = (1-t)**2*sy + 2*(1-t)*t*mid_y + t**2*ey
            ax.plot(bx, by, color=color, lw=lw, alpha=0.8, zorder=2)
            t2 = 0.92
            ax_pt = (1-t2)**2*sx + 2*(1-t2)*t2*mid_x + t2**2*ex
            ay_pt = (1-t2)**2*sy + 2*(1-t2)*t2*mid_y + t2**2*ey
            ax.annotate('', xy=(ex,ey), xytext=(ax_pt,ay_pt),
                        arrowprops=dict(arrowstyle='->', color=color,
                                        lw=lw, mutation_scale=12+lw*2), zorder=3)

    for ct in celltypes:
        x, y = pos[ct]
        ax.add_patch(mpatches.Circle((x,y), NODE_RADIUS, color=get_color(ct),
                                      zorder=5, ec='white', lw=2))
        lx = 1.22*np.cos(angles[ct]); ly = 1.22*np.sin(angles[ct])
        ha = 'left' if lx>0.1 else ('right' if lx<-0.1 else 'center')
        va = 'bottom' if ly>0.1 else ('top' if ly<-0.1 else 'center')
        label = CELLTYPE_LABELS.get(ct, ct)
        ax.text(lx, ly, label, ha=ha, va=va, fontsize=14, fontweight='bold',
                color=get_color(ct),
                path_effects=[pe.withStroke(linewidth=2, foreground='white')])

    ax.set_xlim(-1.6, 1.6); ax.set_ylim(-1.6, 1.6)
    ax.set_aspect('equal'); ax.axis('off')
    ax.set_title(region, fontsize=13, fontweight='bold', pad=4)

    if show_cbar:
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        cax = inset_axes(ax, width="4%", height="40%", loc='lower right',
                         bbox_to_anchor=(0.12, 0.02, 1, 1),
                         bbox_transform=ax.transAxes, borderpad=0)
        sm = ScalarMappable(cmap=cmap, norm=use_norm); sm.set_array([])
        cbar = plt.colorbar(sm, cax=cax)
        cbar.set_label('Age Effect', fontsize=7)
        cbar.ax.tick_params(labelsize=6)

# ── GRID MODE ──────────────────────────────────────────────────────────────
if args.grid:
    ncols = 3; nrows = 4
    n_regions = len(GRID_REGIONS)

    if args.colorbar == 1:
        fig = plt.figure(figsize=(ncols*7 + 1.2, nrows*6))
        gs  = fig.add_gridspec(nrows, ncols+1, width_ratios=[7]*ncols + [0.6],
                                hspace=0.05, wspace=0.02)
        axes = [fig.add_subplot(gs[r, c]) for r in range(nrows) for c in range(ncols)]
        cbar_ax = fig.add_subplot(gs[:, ncols])
    else:
        fig = plt.figure(figsize=(ncols*7, nrows*6))
        gs  = fig.add_gridspec(nrows, ncols, hspace=0.05, wspace=0.02)
        axes = [fig.add_subplot(gs[r, c]) for r in range(nrows) for c in range(ncols)]
        cbar_ax = None

    for i, ax in enumerate(axes):
        if i >= n_regions:
            ax.axis('off'); continue
        region = GRID_REGIONS[i]
        if region not in all_edges:
            ax.axis('off'); continue
        if args.colorbar == 1:
            use_norm = global_norm; show_cbar = False
        elif args.colorbar == 0:
            use_norm = global_norm; show_cbar = False
        else:
            coef_abs_max = max(all_edges[region]['mean_coef'].abs().max(), 0.001)
            use_norm = TwoSlopeNorm(vmin=-coef_abs_max, vcenter=0, vmax=coef_abs_max)
            show_cbar = True
        draw_on_ax(ax, region, all_edges[region], use_norm, show_cbar=show_cbar)
    if args.colorbar == 1 and cbar_ax is not None:
        sm = ScalarMappable(cmap=cmap, norm=global_norm); sm.set_array([])
        cbar = plt.colorbar(sm, cax=cbar_ax)
        cbar.set_label('Age Effect', fontsize=11)
        cbar.ax.tick_params(labelsize=9)
    fig.suptitle(f'{args.category} Age Effects', fontsize=14, fontweight='bold', y=1.002)
    suffix = f"_colorbar{args.colorbar}" if args.colorbar in [0,1] else "_percbar"
    fname = OUT_DIR / f"grid_chord_{cat_slug}_11regions{suffix}.png"
    fig.savefig(fname, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"✓ {fname.name}")
# ── INDIVIDUAL PLOTS ───────────────────────────────────────────────────────
else:
    for region, edges in all_edges.items():
        if args.colorbar == 1:
            use_norm = global_norm
        elif args.colorbar == 0:
            use_norm = global_norm
        else:
            coef_abs_max = max(edges['mean_coef'].abs().max(), 0.001)
            use_norm = TwoSlopeNorm(vmin=-coef_abs_max, vcenter=0, vmax=coef_abs_max)
        fig, ax = plt.subplots(figsize=(9, 9))
        draw_on_ax(ax, region, edges, use_norm, show_cbar=False)
        if args.colorbar != 0:
            sm = ScalarMappable(cmap=cmap, norm=use_norm); sm.set_array([])
            plt.colorbar(sm, ax=ax, shrink=0.35, pad=0.02).set_label('Age Effect', fontsize=9)
        ax.set_title(f"{region} — {args.category}\n"
                     f"Color=mean age effect | Width=n sig LR pairs (q<0.05)",
                     fontsize=12, fontweight='bold', pad=20)
        plt.tight_layout()
        fname = OUT_DIR / f"{region.lower()}_chord_{cat_slug}.png"
        fig.savefig(fname, dpi=200, bbox_inches='tight')
        plt.close()
        print(f"✓ {fname.name}")
print("Done.")
