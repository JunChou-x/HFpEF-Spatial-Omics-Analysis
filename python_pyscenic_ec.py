import os
import glob
import numpy as np
import scanpy as sc
import pandas as pd
import loompy as lp
import anndata as ad 
import matplotlib.pyplot as plt
import scipy.sparse as sp
import scipy.stats as stats
import seaborn as sns

sns.set(style="whitegrid", font_scale=1.0)

h5ad_path = './pysenic_ec.h5ad'
loom_path = './ec_pyscenic_output.loom'

adata = sc.read_h5ad(h5ad_path)

with lp.connect(loom_path, mode='r', validate=False) as lf:
    auc_matrix = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)

common_cells = adata.obs_names.intersection(auc_matrix.index)
adata_sub = adata[common_cells].copy()
auc_matrix_sub = auc_matrix.loc[common_cells]
adata_sub.obsm['X_aucell'] = auc_matrix_sub

adata_sub.obs = pd.concat([adata_sub.obs, auc_matrix_sub], axis=1).copy()
adata.obs["condition"].unique()

adata_auc = ad.AnnData(X=auc_matrix_sub, obs=adata_sub.obs)

key_markers = {
    "Vascular Endothelial Cells": [
        "Ppara(+)", "Rxrg(+)",
        "Nfe2l1(+)", "Cebpb(+)",
        "Atf6(+)"
    ],
    "Endocardial Endothelial Cells": [
        "Cebpb(+)",
        "Sox10(+)", "Gli3(+)"
    ],
    "Lymphatic Endothelial Cells": [
        "Irf8(+)", "Rara(+)", "Srebf2(+)",
        "Nfe2l1(+)", "Cebpb(+)",
        "Foxc1(+)"
    ]
}

key_markers = {
    "Vascular Endothelial Cells": [
        "Ppara(+)", "Rxrg(+)",
        "Nfe2l1(+)", "Cebpb(+)",
        "Atf6(+)"
    ],
    "Endocardial Endothelial Cells": [
        "Cebpb(+)",
        "Sox10(+)", "Gli3(+)"
    ],
    "Lymphatic Endothelial Cells": [
        "Irf8(+)", "Rara(+)", "Srebf2(+)",
        "Nfe2l1(+)", "Cebpb(+)",
        "Foxc1(+)"
    ]
}

def plot_volcano(adata, cluster_name, highlight_genes, ax):
    subset = adata[adata.obs['L5'] == cluster_name].copy()
    
    try:
        sc.tl.rank_genes_groups(subset, groupby='condition', reference='Control', method='wilcoxon')
        df = sc.get.rank_genes_groups_df(subset, group='HFpEF')
    except:
        print(f"skip {cluster_name}: ")
        return


    min_pval = df[df['pvals_adj'] > 0]['pvals_adj'].min()
    df['pvals_adj'] = df['pvals_adj'].replace(0, min_pval / 10)
    df['-log10(padj)'] = -np.log10(df['pvals_adj'])

    df['color'] = 'ns'
    df.loc[(df['logfoldchanges'] > 0.1) & (df['pvals_adj'] < 0.05), 'color'] = 'up'
    df.loc[(df['logfoldchanges'] < -0.1) & (df['pvals_adj'] < 0.05), 'color'] = 'down'

    sns.scatterplot(
        data=df, x='logfoldchanges', y='-log10(padj)',
        hue='color',
        palette={'ns': 'lightgrey', 'up': '#E63946', 'down': '#4E84C4'},
        alpha=0.6, s=40, legend=False, ax=ax
    )

    texts = []
    to_label = df[df['names'].isin(highlight_genes) & (df['color'] != 'ns')]
    for _, row in to_label.iterrows():
        ax.scatter(row['logfoldchanges'], row['-log10(padj)'], color='black', s=5, zorder=10)
        texts.append(ax.text(row['logfoldchanges'], row['-log10(padj)'], row['names'],
                             fontsize=10, fontweight='bold', color='black'))
    if texts:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='black', lw=0.5))

    ax.set_title(cluster_name, fontsize=14, fontweight='bold')
    ax.set_xlabel("log2 Fold Change (HFpEF/Control)")
    ax.set_ylabel("-log10 Adj. P-value")
    ax.axvline(x=0, color='grey', linestyle='--', linewidth=1)
    ax.axhline(y=-np.log10(0.05), color='grey', linestyle='--', linewidth=1)

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

plt.tight_layout()
plt.show()

pdf_path = "./HFpEF_volcano_TF_three_clusters.pdf"
fig.savefig(pdf_path, bbox_inches="tight")
