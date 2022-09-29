from sys import argv
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

basedir = argv[1]

sns.set_context('paper')
sns.set(style = 'whitegrid')
fig, ax = plt.subplots(figsize = (12, 8))

df = pd.read_csv(f'{basedir}/results/reported_alignments.tsv', sep = '\t')
order = list(set(df['Common name']))
order.sort()
sns.barplot(data = df, x = 'Common name', y = "Alignments", order = order, errorbar = 'se')
plt.xticks(rotation = 45)
plt.savefig(f'{basedir}/report/graphs/alignments_barplot_.png', dpi = 600)
plt.clf()

sns.barplot(data = df, x = 'Common name', y = "Normcounts", order = order, errorbar = 'se')
plt.xticks(rotation = 45)
plt.ylabel('Normalized alignments count')
plt.savefig(f'{basedir}/report/graphs/alignments_barplot_norm.png', dpi = 600)
plt.clf()
