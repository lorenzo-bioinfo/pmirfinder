import pandas as pd
import sys
import os

basedir = sys.argv[1]
summaryfile = os.path.join(basedir, 'results/final_summary.tsv')
df = pd.read_csv(summaryfile, header = 0, index_col = 0, sep = '\t')
seqs = list(df['Sequence'])
outfile = f'{basedir}/results/final_summary.fasta'
print(outfile)
with open(outfile, 'w') as f:
	for seq in seqs:
		data = list(df[df['Sequence'] == seq].values[0])
		f.write(f'>{data[2]}\n{seq}\n')
print(data)
