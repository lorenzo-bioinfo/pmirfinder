from sys import argv
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from scipy.cluster import hierarchy
from matplotlib import pyplot as plt
import seaborn as sns

sns.set_context('paper')

basedir = argv[1]

counts_df = pd.read_csv(f'{basedir}/report/counts_matrix.tsv', sep = '\t', header = 0, index_col = 0, na_values = 0).dropna(how='all', axis = 0).fillna(0)
pa_df = pd.read_csv(f'{basedir}/report/pres_abs_matrix.tsv', sep = '\t', header = 0, index_col = 0, na_values = 0).dropna(how='all', axis = 0).fillna(0)

#normalizing, clustering and plotting counts data
fitter = MinMaxScaler().fit(counts_df)
norm_counts_df = pd.DataFrame(fitter.transform(counts_df))
norm_counts_df.index = counts_df.index
norm_counts_df.columns = counts_df.columns
cluster_col = hierarchy.linkage(norm_counts_df.T, method="ward", metric="euclidean")
cluster_row = hierarchy.linkage(norm_counts_df, method="ward", metric="euclidean")
clusterfig = sns.clustermap(norm_counts_df, row_linkage = cluster_row, col_linkage = cluster_col, cmap = 'magma')
index_col = clusterfig.dendrogram_col.reordered_ind #sequences
index_row = clusterfig.dendrogram_row.reordered_ind #srrids
#exporting ordered lists
sequences = list(norm_counts_df.columns)
ids = list(norm_counts_df.index)
with open(f'{basedir}/results/sequences_counts_ordered.txt', 'w') as f:
    for i in index_col:
        f.write(f'{sequences[i]}\n')
with open(f'{basedir}/results/srrids_counts_ordered.txt', 'w') as f:
    for i in index_row:
        f.write(f'{ids[i]}\n')
plt.savefig(f'{basedir}/report/graphs/counts_clustering.png', dpi = 300)
plt.clf()
clusterfig = sns.clustermap(norm_counts_df, row_linkage = cluster_row, col_linkage = cluster_col, cmap = 'magma', figsize = (len(norm_counts_df) / 2, len(norm_counts_df.columns) / 2))
plt.savefig(f'{basedir}/report/graphs/big_counts_clustering.png', dpi = 300)
plt.clf()

#clustering and plotting presence/absence data

cluster_col = hierarchy.linkage(pa_df.T, method="ward", metric="euclidean")
cluster_row = hierarchy.linkage(pa_df, method="ward", metric="euclidean")
clusterfig = sns.clustermap(pa_df, row_linkage = cluster_row, col_linkage = cluster_col, cmap = 'magma')
index_col = clusterfig.dendrogram_col.reordered_ind #sequences
index_row = clusterfig.dendrogram_row.reordered_ind #srrids
#exporting ordered lists
sequences = list(pa_df.columns)
ids = list(pa_df.index)
with open(f'{basedir}/results/sequences_presabs_ordered.txt', 'w') as f:
    for i in index_col:
        f.write(f'{sequences[i]}\n')
with open(f'{basedir}/results/srrids_presabs_ordered.txt', 'w') as f:
    for i in index_row:
        f.write(f'{ids[i]}\n')
plt.savefig(f'{basedir}/report/graphs/pres_abs_clustering.png', dpi = 300)
plt.clf()
clusterfig = sns.clustermap(pa_df, row_linkage = cluster_row, col_linkage = cluster_col, cmap = 'magma', figsize = (len(pa_df) / 2, len(pa_df.columns) / 2))
plt.savefig(f'{basedir}/report/graphs/big_pres_abs_clustering.png', dpi = 300)
plt.clf()
