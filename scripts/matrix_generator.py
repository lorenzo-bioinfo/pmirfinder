import os
import pandas as pd
import numpy as np
from sys import argv

#defining function for parsing fasta files
def fastparse(filename):
    diz = {}
    try:
        with open(filename, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    count = int(line.strip().split('|')[2].split('-')[1])
                else:
                    diz[line.strip()] = count
        return diz
    except FileNotFoundError:
        return diz


#getting basedir
basedir = argv[1]
#listing projects dir
projdir = os.listdir(f'{basedir}/data/projects/')
#getting projects names
projects = [x.split('.')[0] for x in projdir]
#getting plants aliases
plantsdir = os.listdir(f'{basedir}/data/genomes/')
plants = [x.split('.')[0] for x in plantsdir]
#getting list of all srrids
with open(f'{basedir}/data/allids.txt', 'r') as f:
    srrids = [line.strip() for line in f]
#getting individual sequences
df = pd.read_csv(f'{basedir}/report/final_report.tsv', sep = '\t')
pmirnas = list(df['Sequence'])

#initializing matrices
seqs_dict = dict(zip(pmirnas, np.zeros(len(pmirnas) - 1)))
counts_matrix = []
pa_matrix = [] #presence/absence

#parsing files and populating matrices
for srrid in srrids:
    sr_counts = dict(zip(pmirnas, np.zeros(len(pmirnas) - 1)))
    sr_pa = dict(zip(pmirnas, np.zeros(len(pmirnas) - 1)))
    for proj in projects:
        for plant in plants:
            seqs = fastparse(f'{basedir}/{proj}/{proj}.pmirs/{srrid}.{plant}.fa')
            for seq in seqs:
                try:
                    sr_counts[seq] += seqs[seq]
                    sr_pa[seq] = 1
                except KeyError:
                    pass
    counts_matrix.append(sr_counts)
    pa_matrix.append(sr_pa)

counts_df = pd.DataFrame.from_records(data = counts_matrix)
counts_df.index = srrids
counts_df = counts_df.applymap(lambda x: x/len(plants))
pa_df = pd.DataFrame.from_records(data = pa_matrix)
pa_df.index = srrids
counts_df.to_csv(f'{basedir}/report/counts_matrix.tsv', sep = '\t')
pa_df.to_csv(f'{basedir}/report/pres_abs_matrix.tsv', sep = '\t')
