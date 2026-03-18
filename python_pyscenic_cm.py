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
h5ad_path = './pysenic_cm.h5ad'
loom_path = './cm_pyscenic_output.loom'
adata = sc.read_h5ad(h5ad_path)

with lp.connect(loom_path, mode='r', validate=False) as lf:
    auc_matrix = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)

common_cells = adata.obs_names.intersection(auc_matrix.index)
adata_sub = adata[common_cells].copy()
auc_matrix_sub = auc_matrix.loc[common_cells]
adata_sub.obsm['X_aucell'] = auc_matrix_sub   
adata_sub.obs = pd.concat([adata_sub.obs, auc_matrix_sub], axis=1).copy()

target_cell = 'Stressed Cardiomyocytes'
adata_subset = adata_auc[adata_auc.obs['L5'] == target_cell].copy()

sc.tl.rank_genes_groups(adata_subset, groupby='condition', reference='Control', method='wilcoxon')
diff_results = sc.get.rank_genes_groups_df(adata_subset, group='HFpEF')

up_tfs = diff_results[
    (diff_results['pvals_adj'] < 0.05) & 
    (diff_results['logfoldchanges'] > 0.1)
].sort_values('logfoldchanges', ascending=False)

down_tfs = diff_results[
    (diff_results['pvals_adj'] < 0.05) & 
    (diff_results['logfoldchanges'] < -0.1)
].sort_values('logfoldchanges', ascending=True)

print("\n" + "="*60)
print(f"  Top Upregulated (LogFC > 0.1)")
print("="*60)
print(up_tfs[['names', 'logfoldchanges', 'pvals_adj']].head(10))

print("\n" + "="*60)
print(f"  Top Downregulated (LogFC < -0.1)")
print("="*60)
print(down_tfs[['names', 'logfoldchanges', 'pvals_adj']].head(10))

logfc_threshold = 0.1
break_start = -2.5  
break_end = -1.5    

df = diff_results.copy()

manual_genes = ['Nfkb1', 'Nr1d2', 'Bhlhe40', 'Wt1', 'Irf6(+)', 
                'Alx4', 'Arid5b', 'Tbx15', 'Cebpg', 'Hinfp', 'Fos', 'Nr4a1']
labels_to_show = set()
for g in manual_genes:
    if g.endswith('(+)'): labels_to_show.add(g)
    else: labels_to_show.add(g + '(+)')

min_nonzero_pval = df.loc[df['pvals_adj'] > 0, 'pvals_adj'].min()
if pd.isna(min_nonzero_pval): min_nonzero_pval = 1e-300
df['pvals_adj'] = df['pvals_adj'].replace(0, min_nonzero_pval * 0.1)
df['log10p'] = -np.log10(df['pvals_adj'])

df['group'] = 'NS'
df.loc[(df['logfoldchanges'] > logfc_threshold) & (df['pvals_adj'] < 0.05), 'group'] = 'Up'
df.loc[(df['logfoldchanges'] < -logfc_threshold) & (df['pvals_adj'] < 0.05), 'group'] = 'Down'

fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(8, 8), 
                               gridspec_kw={'width_ratios': [1, 4]})

def plot_on_ax(target_ax):

    sns.scatterplot(data=df[df['group'] == 'NS'], x='logfoldchanges', y='log10p', 
                    color='lightgrey', alpha=0.5, s=20, ax=target_ax, label='NS')
    sns.scatterplot(data=df[df['group'] == 'Up'], x='logfoldchanges', y='log10p', 
                    color='#E64B35', alpha=0.8, s=50, ax=target_ax, label='Upregulated')
    sns.scatterplot(data=df[df['group'] == 'Down'], x='logfoldchanges', y='log10p', 
                    color='#3C5488', alpha=0.8, s=50, ax=target_ax, label='Downregulated')
    
    target_ax.axvline(x=0, linestyle='--', color='grey', linewidth=0.8)
    target_ax.axvline(x=logfc_threshold, linestyle='--', color='grey', linewidth=0.8, alpha=0.6)
    target_ax.axvline(x=-logfc_threshold, linestyle='--', color='grey', linewidth=0.8, alpha=0.6)
    target_ax.axhline(y=-np.log10(0.05), linestyle='--', color='grey', linewidth=0.8)
    target_ax.grid(False)

plot_on_ax(ax1)
plot_on_ax(ax2)

x_min_total = df['logfoldchanges'].min() - 0.5
x_max_total = df['logfoldchanges'].max() + 0.5
ax1.set_xlim(x_min_total, break_start)
ax2.set_xlim(break_end, x_max_total)

ax1.spines['top'].set_visible(True)
ax1.spines['bottom'].set_visible(True)
ax1.spines['left'].set_visible(True)
ax1.spines['right'].set_visible(False) 

ax2.spines['top'].set_visible(True)
ax2.spines['bottom'].set_visible(True)
ax2.spines['right'].set_visible(True)
ax2.spines['left'].set_visible(False) 

for ax in [ax1, ax2]:
    for spine in ax.spines.values():
        spine.set_color('black')
        spine.set_linewidth(1.0)
    ax.tick_params(axis='both', colors='black')

ax1.tick_params(axis='y', which='both', right=False)
ax2.tick_params(axis='y', which='both', left=False, labelleft=False)

d = 0.015 
kwargs = dict(transform=ax1.transAxes, color='black', clip_on=False)
ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs) 
ax1.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs) 
kwargs.update(transform=ax2.transAxes)
ax2.plot((-d, +d), (-d, +d), **kwargs) 
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs) 

texts1, texts2 = [], []
for i, row in df.iterrows():
    if row['names'] in labels_to_show:
        x_val, y_val = row['logfoldchanges'], row['log10p']
        if x_val <= break_start:
            texts1.append(ax1.text(x_val, y_val, row['names'], fontsize=10, fontstyle='italic'))
        elif x_val >= break_end:
            texts2.append(ax2.text(x_val, y_val, row['names'], fontsize=10, fontstyle='italic'))

if texts1: adjust_text(texts1, ax=ax1, arrowprops=dict(arrowstyle='-', color='grey', lw=0.5), expand_points=(1.5, 1.5))
if texts2: adjust_text(texts2, ax=ax2, arrowprops=dict(arrowstyle='-', color='grey', lw=0.5), expand_points=(1.5, 1.5))

fig.suptitle('Differential Regulon Activity: Stressed CM (HFpEF vs Control)', fontsize=14, y=0.95)
fig.text(0.5, 0.02, 'log2 Fold Change', ha='center', fontsize=12)
ax1.set_ylabel('-log10 (adj. P-value)', fontsize=12)
ax1.set_xlabel('')
ax2.set_xlabel('')

if ax1.get_legend(): ax1.get_legend().remove()
handles, labels = ax2.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax2.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False, borderaxespad=0.)

plt.subplots_adjust(wspace=0.05)
plt.savefig('./Stressed_CM_Volcano_BrokenAxis_0.1.pdf', format='pdf', bbox_inches='tight', dpi=300)
plt.show()



adata_control = adata_auc[adata_auc.obs['condition'] == 'Control'].copy()

target_cells = ['Stressed Cardiomyocytes', 'Ventricular Cardiomyocytes']
adata_control = adata_control[adata_control.obs['L5'].isin(target_cells)].copy()


sc.tl.rank_genes_groups(adata_control, groupby='L5', reference='Ventricular Cardiomyocytes', method='wilcoxon')

diff_results_identity = sc.get.rank_genes_groups_df(adata_control, group='Stressed Cardiomyocytes')

up_tfs_identity = diff_results_identity[
    (diff_results_identity['pvals_adj'] < 0.05) & 
    (diff_results_identity['logfoldchanges'] > 0)
].sort_values('logfoldchanges', ascending=False)

down_tfs_identity = diff_results_identity[
    (diff_results_identity['pvals_adj'] < 0.05) & 
    (diff_results_identity['logfoldchanges'] < 0)
].sort_values('logfoldchanges', ascending=True)

if not up_tfs_identity.empty:
    print(up_tfs_identity[['names', 'logfoldchanges', 'pvals', 'pvals_adj']].head(20))

if not down_tfs_identity.empty:
    print(down_tfs_identity[['names', 'logfoldchanges', 'pvals', 'pvals_adj']].head(20))

df = diff_results_identity.copy()

raw_up_genes = ['Alx4', 'Arid5b', 'Tbx15', 'Cebpg', 'Hinfp', 'Nfatc4', 'Myc']
raw_down_genes = ['Zfp182(+)', 'Zfp212(+)', 'Ebf1(+)', 'Ilf2(+)', 'Zfhx2(+)']

genes_to_label_up = [g + "(+)" for g in raw_up_genes]
genes_to_label_down = raw_down_genes
labels_to_show = set(genes_to_label_up + genes_to_label_down)

min_nonzero_pval = df.loc[df['pvals_adj'] > 0, 'pvals_adj'].min()
if pd.isna(min_nonzero_pval):
    min_nonzero_pval = 1e-300
df['pvals_adj_clean'] = df['pvals_adj'].replace(0, min_nonzero_pval * 0.1)
df['log10p'] = -np.log10(df['pvals_adj_clean'])

logfc_threshold = 0.25
pval_threshold = 0.05

df['group'] = 'NS' 
df.loc[(df['logfoldchanges'] > logfc_threshold) & (df['pvals_adj'] < pval_threshold), 'group'] = 'Stressed Identity'
df.loc[(df['logfoldchanges'] < -logfc_threshold) & (df['pvals_adj'] < pval_threshold), 'group'] = 'Normal Identity'

plt.figure(figsize=(8, 8))


sns.scatterplot(
    data=df[df['group'] == 'NS'],
    x='logfoldchanges', y='log10p',
    color='lightgrey', alpha=0.3, s=15, label='NS'
)


sns.scatterplot(
    data=df[df['group'] == 'Stressed Identity'],
    x='logfoldchanges', y='log10p',
    color='#E64B35', alpha=0.8, s=45, label='Stressed CM Identity'
)


sns.scatterplot(
    data=df[df['group'] == 'Normal Identity'],
    x='logfoldchanges', y='log10p',
    color='#3C5488', alpha=0.8, s=45, label='Normal Ventricular Identity'
)

texts = []
for i, row in df.iterrows():
    if row['names'] in labels_to_show:
        color = 'black' 
        fontweight = 'normal'
        if 'Prrx1' in row['names'] or 'Nfatc4' in row['names']:
             fontweight = 'bold'

        texts.append(plt.text(
            row['logfoldchanges'],
            row['log10p'],
            row['names'],
            fontsize=10,
            fontweight=fontweight,
            color=color,
            fontstyle='italic'
        ))

adjust_text(
    texts,
    arrowprops=dict(arrowstyle='-', color='grey', lw=0.5),
    expand_points=(1.5, 1.5)
)

plt.title('Baseline Identity: Stressed CM vs. Normal Ventricular CM (Control Group)', fontsize=14)
plt.xlabel('log2 Fold Change (Enriched in Stressed <---> Enriched in Normal)', fontsize=12)
plt.ylabel('-log10 (adj. P-value)', fontsize=12)

plt.axvline(x=logfc_threshold, linestyle='--', color='grey', alpha=0.5, lw=0.8)
plt.axvline(x=-logfc_threshold, linestyle='--', color='grey', alpha=0.5, lw=0.8)
plt.axhline(y=-np.log10(pval_threshold), linestyle='--', color='grey', alpha=0.5, lw=0.8)

plt.grid(False)

ax = plt.gca() 
ax.spines['top'].set_visible(True)
ax.spines['right'].set_visible(True)
ax.spines['bottom'].set_visible(True)
ax.spines['left'].set_visible(True)
for spine in ax.spines.values():
    spine.set_color('black')
    spine.set_linewidth(1.0) 


plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., frameon=False)

plt.savefig('./Volcano_Stressed_Identity_Boxed.pdf', format='pdf', bbox_inches='tight')

plt.tight_layout()
plt.show()