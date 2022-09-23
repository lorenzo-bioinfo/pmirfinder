import sys
from gtfparse import read_gtf
import pandas as pd
import os

#reading cli arguments
plant = sys.argv[1]
proj = sys.argv[2]
basedir = sys.argv[3]

#reading srr_ids
srr_ids = []
with open(os.path.join(basedir, f'data/projects/{proj}.txt')) as f:
	for line in f:
		srr_ids.append(line.strip())

#formatting annotation file path
gtf_path = os.path.join(basedir, f'bt-index/{plant}/{plant}.gtf')
#getting gtf annotation file
df = read_gtf(gtf_path)

#filtering miRNAs annotations
mirannot = df[(df['gene_biotype'] == 'pre_miRNA') | (df['gene_biotype'] == 'miRNA')]
idx = list(mirannot.index)

#creating directory
os.path.join(basedir, f'{proj}/proj.pmirs')
if os.path.isdir(os.path.join(basedir, f'{proj}/{proj}.pmirs')):
	print(f'{proj}.pmirs directory found')
else:
	os.mkdir(os.path.join(basedir, f'{proj}/{proj}.pmirs'))

for srrid in srr_ids:
	try:
		#initializing collection variables for data
		annots = [] #contains start and end coordinates and gene_id for mirnas annots
		genes = {} #contains reads that align to mirnas regions and reads count

		#finalizing variables setup
		for i in idx:
			row = mirannot.loc[[i]]
			annots.append([int(row['start'].values[0]), int(row['end'].values[0]), row['gene_id'].values[0]])
			genes[str(row['gene_id'].values[0])] = [0, [], row['gene_biotype'].values[0]]
		#checking sam file for reads that align to mirs coordinates
		with open(f'{basedir}/{proj}/{proj}.{plant}/{srrid}.{plant}.sam', 'r') as f:
			for line in f:
				if line.startswith('SRR'):
					features = line.strip().split('\t')
					if int(features[3]) > 0:
						read_desc = [features[9], int(features[3]), len(features[9])]
						for annot in annots:
							if read_desc[1] >= annot[0] and read_desc[1] + read_desc[2] <= annot[1]:
								genes[annot[2]][0] += 1
								genes[annot[2]][1].append(read_desc[0])
		#saving results in fasta format
		with open(os.path.join(basedir, f'{proj}/{proj}.pmirs/{srrid}.{plant}.fa'), 'w') as f:
			for gene in genes:
				if genes[gene][0] > 0:
					sequences = list(set(genes[gene][1]))
					for seq in sequences:
						if genes[gene][1].count(seq) > 0:
							f.write(f'>{gene}|L-{len(seq)}|C-{genes[gene][1].count(seq)}|biotype-{genes[gene][2]}\n{seq}\n')
	except FileNotFoundError:
		print(f'{srrid} not found. Probably the download went wrong in the first place...')
