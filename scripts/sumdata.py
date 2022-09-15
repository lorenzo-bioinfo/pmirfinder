import sys
import pandas as pd
import os

#reading cli arguments
plant = sys.argv[1]
basedir = sys.argv[2]

#Getting SRA Projects IDs
projs = [x.split('.')[0] for x in os.listdir(os.path.join(basedir, 'data/projects/'))]

plant_counts = {}
for proj in projs:
	summary_path = os.path.join(basedir, f'{proj}/{proj}.summary/{proj}.{plant}_summary.tsv')
	df = pd.read_csv(summary_path, sep = '\t')
	full = True
	try:
		seqs = list(df['Seq'])
	except KeyError:
		print(f'{proj} file for {plant} was empty')
		full = False
	if full:
		for seq in seqs:
			if seq in (plant_counts.keys()):
				plant_counts[seq][0] += df[df['Seq'] == seq]['ProjMeanCount'].values[0]
				plant_counts[seq][1] += df[df['Seq'] == seq]['ProjNormCount'].values[0]
			else:
				normcounts = df[df['Seq'] == seq]['ProjNormCount'].values[0]
				counts = df[df['Seq'] == seq]['ProjMeanCount'].values[0]
				plant_counts[seq] = [counts, normcounts]
clean_dict = {}
for i, key in enumerate(plant_counts):
	clean_dict[str(i)] = (key, plant_counts[key][0]/len(projs), plant_counts[key][1]/len(projs))
df_to_export = pd.DataFrame.from_dict(clean_dict, orient = 'index', columns = ['Seq', 'MeanCount', 'NormCount'])
if os.path.isdir(os.path.join(basedir, 'results')):
	print('Results folder found')
	pass
else:
	os.mkdir(os.path.join(basedir, 'results'))
	print('Created results folder')
exp_file = os.path.join(basedir, f'results/{plant}_results.tsv')
df_to_export.sort_values(by = 'MeanCount', ascending = False).to_csv(exp_file, sep = '\t')
