from sys import argv
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

basedir = argv[1]

sns.set_context('paper')
df = pd.read_csv(f'{basedir}/report/final_report.tsv', sep = '\t')
df.columns = ['N', 'Sequence', 'RPM', 'SeqID', 'Best P-match', 'Plant miRs',
       'Best HSa-match', 'Hsa miRs', 'Best O-match', 'Other miRs',
       'Best-HRefSeq', 'Hsa RefSeq', 'RefSeq-Ident', 'RefSeq-Gaps', 'Individuals']
print(df)
pbit = list(df['Plant miRs'].fillna(0.0))
hbit = list(df['Hsa miRs'].fillna(0.0))
obit = list(df['Other miRs'].fillna(0.0))
rsbit = list(df['Hsa RefSeq'].fillna(0.0))
etiquettes = []
etiquettes.extend(['Plant miRs' for x in pbit])
etiquettes.extend(['Hsa miRs' for x in hbit])
etiquettes.extend(['Other miRs' for x in obit])
etiquettes.extend(['Hsa RefSeq' for x in rsbit])
scores = []
scores.extend(pbit)
scores.extend(hbit)
scores.extend(obit)
scores.extend(rsbit)
dataf = pd.DataFrame(zip(scores, etiquettes), columns = ['Bitscore', 'Database'])
sns.set(style = 'whitegrid')
sns.kdeplot(data = dataf, x = 'Bitscore', hue = 'Database', fill = True, cumulative = False)
plt.ylabel('')
plt.savefig(f'{basedir}/report/graphs/clean_bitscores_distribution.png', dpi = 600)
plt.clf()
