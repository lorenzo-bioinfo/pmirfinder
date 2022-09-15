import pandas as pd
import sys
basedir = sys.argv[1]

srrids = []
with open(f'{basedir}/data/allids.txt', 'r') as f:
	for line in f:
		srrids.append(line.strip())

resultsdf = pd.read_csv(f'{basedir}/results/risultati_finali.tsv', sep = '\t')
sequences = list(resultsdf['Sequence'])

ind_seqs = {}
for srrid in srrids:
	s = []
	with open(f'{basedir}/results/individuals/{srrid}.txt', 'r') as f:
		for line in f:
			s.append(line.strip())
		ind_seqs[srrid] = set(s)

counts = []
for seq in sequences:
	c = 0
	for ind in ind_seqs:
		if seq in ind_seqs[ind]:
			c += 1
	counts.append(c)

resultsdf['Individuals'] = counts
resultsdf.to_csv(f'{basedir}/results/risultati_full.tsv', sep = '\t')
