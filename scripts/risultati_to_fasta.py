import pandas as pd
df = pd.read_csv('risultati_pmirfinder_puliti.tsv', sep = '\t')
print(df.columns)
ids = list(df['SeqID'])
seqs = list(df['Sequence'])
with open('risultati_definitivi.fasta', 'w') as f:
    for i in range(len(ids)):
        f.write(f'>{ids[i]}\n{seqs[i]}\n')