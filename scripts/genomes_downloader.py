from sys import argv
import os
from wget import download
from shutil import copy
from subprocess import run

basedir = argv[1]

files_list = os.listdir(os.path.join(basedir, 'data/genomes/'))
plants_aliases = [x.split('.')[0] for x in files_list]

#creating directories and downloading genomes
if os.path.isdir(os.path.join(basedir, 'bt-index/')):
	pass
else:
	os.mkdir(os.path.join(basedir, 'bt-index/'))
os.chdir(os.path.join(basedir, 'bt-index/'))

#downloading human references for filtering
#human genome
print('\nDownloading human genome')
os.mkdir('hgenome')
os.chdir(os.path.join(basedir, 'bt-index/hgenome'))
download('http://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz')
print('\nDecompressing genome file')
run(['gunzip', 'Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'])
os.rename('Homo_sapiens.GRCh38.dna.primary_assembly.fa', 'hgenome.fa')
#human cdna
os.chdir(os.path.join(basedir, 'bt-index/'))
print('\nDownloading human cDNA')
os.mkdir('hcdna')
os.chdir(os.path.join(basedir, 'bt-index/hcdna'))
download('http://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz')
print('\nDecompressing file')
run(['gunzip', 'Homo_sapiens.GRCh38.cdna.all.fa.gz'])
os.rename('Homo_sapiens.GRCh38.cdna.all.fa', 'hcdna.fa')
#human nonchromosomal
os.chdir(os.path.join(basedir, 'bt-index/'))
print('\nDownloading human non-chromosomal DNA')
os.mkdir('hnonchr')
os.chdir(os.path.join(basedir, 'bt-index/hnonchr'))
download('http://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.nonchromosomal.fa.gz')
print('\nDecompressing file')
run(['gunzip', 'Homo_sapiens.GRCh38.dna.nonchromosomal.fa.gz'])
os.rename('Homo_sapiens.GRCh38.dna.nonchromosomal.fa', 'hnonchr.fa')
os.chdir(os.path.join(basedir, 'bt-index/'))
#human refseq mrna
print('Downloading human RefSeq mRNA data')
os.mkdir(os.path.join(basedir, 'blast_db/'))
os.mkdir(os.path.join(basedir, 'blast_db/hrefseq/'))
os.chdir(os.path.join(basedir, 'blast_db/hrefseq'))
for i in range(1, 7):
	download(f'https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/human.{i}.rna.fna.gz')
print('\nDecompressing and merging data...')
ziplist = os.listdir()
for zipped in ziplist:
	run(['gunzip', zipped])
flist = os.listdir()
with open('hrefseq.fna', 'w') as f:
	for file in flist:
		with open(file) as temp:
			f.write(temp.read())
		os.remove(file)
print('Done')
os.chdir(os.path.join(basedir, 'bt-index'))
for file in files_list:
	plant_alias = file.split('.')[0]
	os.mkdir(plant_alias)
	os.chdir(os.path.join(basedir, f'bt-index/{plant_alias}/'))
	with open(os.path.join(basedir, f'data/genomes/{file}'), 'r') as f:
		print(f'\nStarted download of {plant_alias} genome')
		for line in f:
			download(line)
		print(f'\nFinished downloading {plant_alias} genome')
		print('Decompressing files...')
		gzip_files = os.listdir()
		for zipped in gzip_files:
			run(['gunzip', f'{zipped}'])
		print(f'Starting assembly of {plant_alias} genome from chromosome files...')
		chromosome_files_temp = os.listdir()
		chromosome_files = []
		for i in range(len(chromosome_files_temp)):
			chromosome_files.append('x')
		for filename in chromosome_files_temp:
			pos = filename.split('.')[-2]
			try:
				chromosome_files[int(pos) - 1] = filename
			except ValueError:
				pos = int(pos[1:]) - 1
				chromosome_files[pos] = filename
		with open(f'{plant_alias}.fa', 'w') as f:
			for chromosome_file in chromosome_files:
				with open(chromosome_file) as temp_file:
					f.write(temp_file.read())
				os.remove(chromosome_file)
		print('Done')
	os.chdir(os.path.join(basedir, 'bt-index/'))
for alias in plants_aliases:
	os.chdir(os.path.join(basedir, f'bt-index/{alias}/'))
	with open(os.path.join(basedir, f'data/gtfs/{alias}.txt')) as f:
		print(f'\nDownloading gtf annotation file for {alias}')
		dl_link = f.readline()
		download(dl_link)
		fname = dl_link.split('/')[-1][:-3]
		print('\nDecompressing and renaming gtf file...')
		run(['gunzip', dl_link.split('/')[-1]])
		os.rename(fname, f'{alias}.gtf')
		print('Done')
