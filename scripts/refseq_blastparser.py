import pandas as pd
import sys
import os

basedir = sys.argv[1]

#reading results
results = pd.read_csv(f'{basedir}/results/final_results_ok.tsv', sep = '\t', header = 0, index_col = 0)
print(results)
pmir_ids = list(results['SeqID'])

#reading blast output
blast_out_header = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'.split(' ')
print(blast_out_header)

#reading final blast results
blastres = pd.read_csv(f'{basedir}/results/refseq_blast_final.tsv', sep = '\t')
blastres.columns = blast_out_header
bestmatch = {}
for mir in pmir_ids:
	mirdf = blastres[blastres['qseqid'] == mir].sort_values(by = 'bitscore', ascending = False)
	try:
		i = 0
		searching = True
		while searching:
			sid = mirdf.iloc[[i]]['sseqid'].values[0]
			if sid.startswith('X'):
				i += 1
				print(f'Skipping predicted RNA {sid}')
			else:
				bestmatch[mir] = [mirdf.iloc[[i]]['sseqid'].values[0], mirdf.iloc[[i]]['bitscore'].values[0], mirdf.iloc[[i]]['pident'].values[0], mirdf.iloc[[i]]['gapopen'].values[0]]
				searching = False	
	except IndexError:
		bestmatch[mir] = [None, None, None, None]
		i = 0
match = []
score = []
pident = []
gapopen = []
for mir in bestmatch:
	match.append(bestmatch[mir][0])
	score.append(bestmatch[mir][1])
	pident.append(bestmatch[mir][2])
	gapopen.append(bestmatch[mir][3])
results['Best-HRefSeq'] = match
results['HRefSeq-BitScore'] = score
results['RefSeq-Ident'] = pident
results['RefSeq-Gaps'] = gapopen
results.to_csv(f'{basedir}/results/risultati_finali.tsv', sep = '\t')
