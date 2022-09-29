import pandas as pd
from sys import argv

basedir = argv[1]


df = pd.read_csv(f'{basedir}/results/risultati_full.tsv', header = 0, index_col = 0, sep = '\t')
print(df)
print(df.columns)
df_norefseq = df[(df['RefSeq-Ident'] < 100) | (df['RefSeq-Ident'].isna())]
print(df_norefseq.columns)
df_clean = df_norefseq[((df_norefseq['HRefSeq-BitScore'] < df_norefseq['P-Bitscore']) & (df_norefseq['P-Bitscore'] > df_norefseq['Hsa-Bitscore']) & (df_norefseq['P-Bitscore'] > df_norefseq['O-Bitscore']) | df_norefseq['HRefSeq-BitScore'].isna())]
df_clean = df_clean[(df_clean['P-Bitscore'] > df_clean['Hsa-Bitscore']) | (df_clean['Hsa-Bitscore'].isna())]
print(df_clean.columns)
df_clean.to_csv(f'{basedir}/report/final_report.tsv', sep = '\t', index = False)
seqs = list(df_clean['Sequence'])
ids = list(df_clean.index)
rpms = list(df_clean['RPM'])
with open(f'{basedir}/report/pmirnas.fasta', 'w') as f:
    for i in range(len(ids)):
        f.write(f'>{ids[i]} | RPM({rpms[i]:.3f})\n{seqs[i]}\n')
