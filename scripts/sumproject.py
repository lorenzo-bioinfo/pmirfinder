import sys
import pandas as pd
import os

#reading cli arguments
plant = sys.argv[1]
proj = sys.argv[2]
basedir = sys.argv[3]
#Starting reading each tsv file and summarizing
#mean read counts for each sequence for each plant
srr_ids = []
with open(os.path.join(basedir, f'data/projects/{proj}.txt')) as f:
	for line in f:
		srr_ids.append(line.strip())

plant_counts = {}
downloaded = len(srr_ids)
for srrid in srr_ids:
	try:
		df = pd.read_csv(os.path.join(basedir, f'{proj}/{proj}.pmirs/{srrid}.{plant}.tsv'), header = 0, sep = '\t')
		full = True
		try:
			seqs = list(df['Seq'])
		except KeyError:
			print(f'{srrid} file for {plant} was empty')
			full = False
		if full:
			for seq in seqs:
				if seq in (plant_counts.keys()):
					plant_counts[seq][0] += df[df['Seq'] == seq]['Count'].values[0]
					plant_counts[seq][1] += df[df['Seq'] == seq]['Normcount'].values[0]
				else:
					normcounts = df[df['Seq'] == seq]['Normcount'].values[0]
					counts = df[df['Seq'] == seq]['Count'].values[0]
					plant_counts[seq] = [counts, normcounts]
	except FileNotFoundError:
		print(f'{srrid} not found. Proably download went wrong in the first place...')
		downloaded = downloaded - 1
clean_dict = {}
for i, key in enumerate(plant_counts):
	clean_dict[str(i)] = (key, plant_counts[key][0]/downloaded, plant_counts[key][1]/downloaded)
df_to_export = pd.DataFrame.from_dict(clean_dict, orient = 'index', columns = ['Seq', 'ProjMeanCount', 'ProjNormCount'])

exp_dir = os.path.join(basedir, f'{proj}/{proj}.summary')
try:
	os.mkdir(exp_dir)
except FileExistsError:
	pass
df_to_export.sort_values(by = 'ProjMeanCount', ascending = False).to_csv(f'{exp_dir}/{proj}.{plant}_summary.tsv', sep = '\t')
