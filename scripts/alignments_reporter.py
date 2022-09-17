from sys import argv
from os import listdir
import pandas as pd

basedir = argv[1]
#getting basic data
#projects list
ls_projects = listdir(f'{basedir}/data/projects/')
projs = [x.split('.')[0] for x in ls_projects]

#plants list
ls_plants = listdir(f'{basedir}/data/genomes/')
plants = [x.split('.')[0] for x in ls_plants]

#dictionary for plants common name
common_names = {}
with open(f'{basedir}/data/common_names.txt', 'r') as f:
	for line in f:
		key, value = line.strip().split(',')
		common_names[key] = value

#genomes size
with open(f'{basedir}/data/genome_sizes.txt', 'r') as f:
	sizes = []
	for line in f:
		size, p_alias = line.split(' ')
		p_alias_ok = p_alias.strip().split('/')[-1][0:-3]
		sizes.append((p_alias_ok, int(size)))
size_dict = dict(sizes)

#getting alignments count for each plant in every project
alignment_count = []
for proj in projs:
	for plant in plants:
		ls_reports = listdir(f'{basedir}/{proj}/{proj}.{plant}/')
		reports_files = [x for x in ls_reports if x.endswith('txt')]
		for file in reports_files:
			with open(f'{basedir}/{proj}/{proj}.{plant}/{file}', 'r') as f:
				for line in f:
					if line.startswith('Reported'):
						n_al = int(line.split(' ')[1])
						alignment_count.append((plant, n_al))
#creating dataframe
report_df = pd.DataFrame.from_records(alignment_count, columns = ['Alias', 'Alignments'])
#adding common names column
alias_order = list(report_df['Alias'])
common_column = [common_names[x] for x in alias_order]
report_df['Common name'] = common_column
#adding normalized counts column
norm_counts = []
counts = list(report_df['Alignments'])
df_dict = dict(zip(alias_order, counts))
for al in alias_order:
	norm_counts.append(df_dict[al] / (size_dict[al] / 1000000000))
report_df['Normcounts'] = norm_counts
#exporting data
report_df.to_csv(f'{basedir}/results/reported_alignments.tsv', sep = '\t')
