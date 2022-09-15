import pandas as pd
import sys
import os

basedir = sys.argv[1]

#reading results
resfile = os.path.join(basedir, 'results/final_summary.tsv')
results = pd.read_csv(resfile, sep = '\t', header = 0, index_col = 0)
pmir_ids = list(results['SeqID'])

#obtaining mirbase ids for plants and other organisms

pids = []
with open(f'{basedir}/data/pids.txt') as f:
	for line in f:
		pids.append(line.strip().lower())
oids = []
with open(f'{basedir}/data/oids.txt') as f:
	for line in f:
		oids.append(line.strip().lower())
blast_out_header = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'.split(' ')
print(blast_out_header)

#reading final blast results
blastres = pd.read_csv(f'{basedir}/results/final_blastres.tsv', sep = '\t')
blastres.columns = blast_out_header
ids_to_clean = list(blastres['sseqid'])
blast_ids = []
for i in ids_to_clean:
	blast_ids.append(i.split('-')[0])
blastres['spid'] = blast_ids
best_pmatch = {}
best_npmatch = {}
best_hmatch = {}
for mir in pmir_ids:
	mirdf = blastres[blastres['qseqid'] == mir]
	pmirdf = mirdf[mirdf['spid'].isin(pids)].sort_values(by = 'bitscore', ascending = False)
	hdf = mirdf[mirdf['spid'] == 'hsa'].sort_values(by = 'bitscore', ascending = False)
	odf = mirdf[mirdf['spid'].isin(oids)].sort_values(by = 'bitscore', ascending = False)
	try:
		best_pmatch[mir] = [pmirdf.iloc[[0]]['sseqid'].values[0], pmirdf.iloc[[0]]['bitscore'].values[0]]
	except IndexError:
		best_pmatch[mir] = [None, None]
	try:
		best_npmatch[mir] = [odf.iloc[[0]]['sseqid'].values[0], odf.iloc[[0]]['bitscore'].values[0]]
	except IndexError:
		best_npmatch[mir] = [None, None]
	try:
		best_hmatch[mir] = [hdf.iloc[[0]]['sseqid'].values[0], hdf.iloc[[0]]['bitscore'].values[0]]
	except IndexError:
		best_hmatch[mir] = [None, None]
pmatch = []
pbit = []
hmatch = []
hbit = []
omatch = []
obit = []
for mir in pmir_ids:
	pmatch.append(best_pmatch[mir][0])
	pbit.append(best_pmatch[mir][1])
	hmatch.append(best_hmatch[mir][0])
	hbit.append(best_hmatch[mir][1])
	omatch.append(best_npmatch[mir][0])
	obit.append(best_npmatch[mir][1])
results['Best P-match'] = pmatch
results['P-Bitscore'] = pbit
results['Best HSa-match'] = hmatch
results['Hsa-Bitscore'] = hbit
results['Best O-match'] = omatch
results['O-Bitscore'] = obit

results.to_csv(f'{basedir}/results/final_results_ok.tsv', sep = '\t')
