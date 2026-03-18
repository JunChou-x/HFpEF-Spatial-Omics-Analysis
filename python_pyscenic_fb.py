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
from adjustText import adjust_text

sns.set(style="whitegrid", font_scale=1.0)

h5ad_path = './pysenic_fb.h5ad'
loom_path = './fb_pyscenic_output.loom'

adata = sc.read_h5ad(h5ad_path)

with lp.connect(loom_path, mode='r', validate=False) as lf:
    auc_matrix = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)

common_cells = adata.obs_names.intersection(auc_matrix.index)
adata_sub = adata[common_cells].copy()
auc_matrix_sub = auc_matrix.loc[common_cells]
adata_sub.obsm['X_aucell'] = auc_matrix_sub

adata_sub.obs = pd.concat([adata_sub.obs, auc_matrix_sub], axis=1).copy()
adata_auc = ad.AnnData(X=auc_matrix_sub, obs=adata_sub.obs)

target_groups = {
    'Ventricular Fibroblasts': ['Nfe2(+)', 'Atf5(+)', 'Nr1d1(+)', 'Bhlhe40(+)', 'Lmx1a(+)', 'Egr2(+)','Rxrg(+)','Ppara(+)'],
    'Myofibroblasts': ['Esrrb(+)', 'Tal2(+)', 'Tfap2a(+)', 'Wt1(+)', 'Lhx8(+)', 'Foxf1(+)'],
    'Fibro-Immune Interface': ['Esrrb(+)', 'Atf5(+)', 'Nr1d1(+)', 'Nfe2(+)', 'Nfe2l3(+)', 'Wt1(+)'],
    'Stress-Response Fibroblasts': ['Esrrb(+)', 'Foxf1(+)', 'Hmga2(+)', 'Nfe2l1(+)', 'Atf5(+)', 'Tfap2a(+)'],
    'Angiogenic Fibroblasts': ['Ppara(+)', 'Rxrg(+)', 'Vdr(+)','Esrra(+)'], 
    'Atrial Fibroblasts': ['Bhlhe40(+)', 'Gli1(+)', 'Atf5(+)', 'Sox10(+)', 'Foxc2(+)']
}

fig, axes = plt.subplots(2, 3, figsize=(20, 12)) 
axes = axes.flatten() 

logfc_threshold = 0.1
pval_threshold = 0.05


for idx, (cell_type, labels_raw) in enumerate(target_groups.items()):
    ax = axes[idx]
    print(f"[{idx+1}/6] : {cell_type} ...")
    
    labels_to_show = set(labels_raw)

    current_cells = adata_sub.obs[adata_sub.obs['L5'] == cell_type].index
    
    current_auc = auc_matrix_sub.loc[current_cells, :]
    current_obs = adata_sub.obs.loc[current_cells, :]
    
    adata_calc = sc.AnnData(X=current_auc, obs=current_obs)
    
    if 'HFpEF' not in adata_calc.obs['condition'].unique() or 'Control' not in adata_calc.obs['condition'].unique():
        print(f"  Warning: A Control or HFpEF sample is missing in {cell_type}. Skip the analysis.")
        continue

    sc.tl.rank_genes_groups(adata_calc, groupby='condition', reference='Control', method='wilcoxon')
    df = sc.get.rank_genes_groups_df(adata_calc, group='HFpEF')
    
    min_nonzero_pval = df.loc[df['pvals_adj'] > 0, 'pvals_adj'].min()
    if pd.isna(min_nonzero_pval): min_nonzero_pval = 1e-300
    df['pvals_adj_clean'] = df['pvals_adj'].replace(0, min_nonzero_pval * 0.1)
    df['log10p'] = -np.log10(df['pvals_adj_clean'])

    df['group'] = 'NS'
    df.loc[(df['logfoldchanges'] > logfc_threshold) & (df['pvals_adj'] < pval_threshold), 'group'] = 'Up'
    df.loc[(df['logfoldchanges'] < -logfc_threshold) & (df['pvals_adj'] < pval_threshold), 'group'] = 'Down'

    sns.scatterplot(data=df[df['group'] == 'NS'], x='logfoldchanges', y='log10p',
                    color='lightgrey', alpha=0.3, s=15, ax=ax, rasterized=True)
    
    sns.scatterplot(data=df[df['group'] == 'Up'], x='logfoldchanges', y='log10p',
                    color='#E64B35', alpha=0.8, s=40, label='Up in HFpEF', ax=ax)
    
    sns.scatterplot(data=df[df['group'] == 'Down'], x='logfoldchanges', y='log10p',
                    color='#3C5488', alpha=0.8, s=40, label='Down in HFpEF', ax=ax)

    texts = []
    for _, row in df.iterrows():
        if row['names'] in labels_to_show:
            texts.append(ax.text(row['logfoldchanges'], row['log10p'], row['names'],
                                 fontsize=10, fontweight='bold', color='black', fontstyle='italic'))
    
    if texts:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='grey', lw=0.5))

    ax.set_title(cell_type, fontsize=12, fontweight='bold')

    ax.axvline(x=logfc_threshold, linestyle='--', color='grey', alpha=0.5, lw=0.8)
    ax.axvline(x=-logfc_threshold, linestyle='--', color='grey', alpha=0.5, lw=0.8)
    ax.axhline(y=-np.log10(pval_threshold), linestyle='--', color='grey', alpha=0.5, lw=0.8)
    
    ax.grid(False)
    
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color('black')
        spine.set_linewidth(1.0)
    
    if idx >= 3: ax.set_xlabel('log2 Fold Change (HFpEF vs Control)')
    else: ax.set_xlabel('')
    
    if idx % 3 == 0: ax.set_ylabel('-log10 (adj. P-value)')
    else: ax.set_ylabel('')

    ax.legend(loc='best', frameon=False, prop={'size': 8}, markerscale=0.8)

plt.tight_layout()
save_path = './Fibroblasts_Subtypes_Volcano_Panel.pdf'
plt.savefig(save_path, format='pdf', bbox_inches='tight', dpi=300)
plt.show()