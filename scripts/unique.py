import pandas as pd
import sys
import os

basedir = sys.argv[1]

plants = [x.split('.')[0] for x in os.listdir(os.path.join(basedir, 'data/genomes/'))]
seqs = {}
for plant in plants:
	resfile = os.path.join(basedir, f'results/{plant}_results.tsv')
	df = pd.read_csv(resfile, sep = '\t')
	sequences = list(df['Seq'])
	counts = list(df['NormCount'])
	tuples = zip(sequences, counts) 
	for tup in tuples:
		if tup[0] in seqs.keys():
			seqs[tup[0]][0] += tup[1]
			seqs[tup[0]][1] += 1
		else:
			seqs[tup[0]] = [tup[1], 1]
seqs_clean = {}
for i, seq in enumerate(seqs):
	seqs_clean[i] = [seq, seqs[seq][0] / seqs[seq][1]]
final_results = pd.DataFrame.from_dict(seqs_clean, orient = 'index')
final_results.columns = ['Sequence', 'RPM']
final_ordered = final_results.sort_values(by = 'RPM', ascending = False)
id_col = []
for i in range(len(final_ordered)):
	id_col.append(f'pmir-{i}')
final_ordered['SeqID'] = id_col
exp_file = os.path.join(basedir, 'results/final_summary.tsv')
final_ordered.to_csv(exp_file, sep = '\t')
