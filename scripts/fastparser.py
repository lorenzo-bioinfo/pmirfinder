import sys
import os
import pandas as pd

def getMillionReads(proj, plant, srrid, basedir):
	with open(os.path.join(basedir, f'{proj}/{proj}.human_mirs/{srrid}.alignment_report.txt')) as f:
		million_reads = int(f.readline().strip().split(' ')[3]) / 1000000
		return million_reads

def fastaParse(fasta_file, proj, plant, srrid, basedir):
	scaling_factor = getMillionReads(proj, plant, srrid, basedir)
	with open(fasta_file, 'r') as f:
		genes = {}
		for i, line in enumerate(f):
			if line.startswith('>'):
				feats = line.strip().split('|')
				genename = feats[0][1:]
				genes[str(i)] = [genename, int(feats[1].split('-')[1]), int(feats[2].split('-')[1]), int(feats[1].split('-')[1]) / scaling_factor]
			else:
				try:
					genes[str(i-1)].append(line.strip())
				except KeyError:
					print(f'{i} not in dict')
	df = pd.DataFrame.from_dict(genes, orient = 'index')
	try:
		df.columns = ['Name', 'Len', 'Count', 'Normcount', 'Seq']
	except ValueError:
		print(f'Skipping {plant} mirs for {srrid} as no PmiRNAs were found')
	return(df)

#reading cli arguments
plant = sys.argv[1]
proj = sys.argv[2]
basedir = sys.argv[3]


#reading srr_ids
srr_ids = []
with open(os.path.join(basedir, f'data/projects/{proj}.txt')) as f:
	for line in f:
		srr_ids.append(line.strip())

for srrid in srr_ids:
	fastafile = os.path.join(basedir, f'{proj}/{proj}.pmirs/{srrid}.{plant}.fa')
	df = fastaParse(fastafile, proj, plant, srrid, basedir)
	df.to_csv(os.path.join(basedir, f'{proj}/{proj}.pmirs/{srrid}.{plant}.tsv'), sep = '\t')
