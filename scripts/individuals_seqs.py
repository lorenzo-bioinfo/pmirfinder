import os
from sys import argv

basedir = argv[1]

def fastaparse(filepath):
	found = False
	sequences = []
	with open(filepath) as f:
		for line in f:
			if found:
				sequences.append(line.strip())
				found = False
			else:
				found = True
	return(sequences)


ls = os.listdir(f'{basedir}/data/projects/')
projects = [x.split('.')[0] for x in ls]
for proj in projects:
	srrids = []
	with open(f'{basedir}/data/projects/{proj}.txt') as f:
		for line in f:
			srrids.append(line.strip())
	dirlist = os.listdir(f'{basedir}/{proj}/{proj}.pmirs/')
	for srrid in srrids:
		srr_seqs = []
		srrfiles = [x for x in dirlist if (x.startswith(srrid) and x.endswith('.fa'))]
		for file in srrfiles:
			srr_seqs.extend(fastaparse(f'{basedir}/{proj}/{proj}.pmirs/{file}'))
		with open(f'{basedir}/results/individuals/{srrid}.txt', 'w') as f:
			srr_seqs = set(srr_seqs)
			for seq in srr_seqs:
				f.write(f'{seq}\n')
