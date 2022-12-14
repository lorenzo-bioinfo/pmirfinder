## <u>PmiRFinder</u>

This is a tool I developed to identify plant-derived micro-RNAs in human smallRNA-sequencing experiments. It's built as a pipeline of Python3 and BASH scripts/tools.
### Usage
To run the tool with default settings and data options, just clone the repository and run */pmirfinder.sh* on your machine.

Please note that the tool is intended to be ran on a GNU/Linux system (preferably Debian-based) and uses several bioinformatics tools and Python3 libraries (Software requirements are listed below). The installation of the required packages and libraries is left to the user.

Minimum hardware requirements to use this software with default data options are:

- 8-cores CPU;
- 16 GB of RAM;
- 500 GB of free disk space.

A good internet connection speed is also a must as a lot of data needs to be downloaded in the process.

By default, the tool performs the analysis on the sequencing experiments I selected for my study and the search is done on selected plant genomes. The instructions to expand/change/modify the data to be used are detailed below.

Please note that sequencing data must come from NCBI-SRA database. If you need to use this tool with local sequencing data, feel free to edit the script to match your needs, but no support is provided.

The steps of the analysis can be summarized as follow:

1) Quality control and trimming of raw reads files;
2) Alignment of reads to Human reference genome, transcriptome, extrachromosomal DNA and human known miRNAs to filter out human sequences;
3) Alignment of the remaining reads to different plant genomes;
4) Selection of reads that aligned in regions annotated as miRNA or pre-miRNA for every plant genome in the corresponding GTF annotation file;
5) Local BLAST search of the sequences identified as putative pmiRNAs against miRBase and Human RefSeq mRNA;
6) Parsing of the obtained results and generation of statistics and graphs.

The tools used in this script are all FOSS, and all the code you find in this repository is distributed under the GNU-GPL v2. If you use any significant chunk of this code, please consider citing the authors of 3rd party tools and my work :)

### - <u>Software requirements</u>
This is a list of the tools and libraries you need to have installed on your machine in order to run pmiRFinder. The installation of these packages is left to the user. The tools must also be in your PATH in order for the script to call them when needed.
#### <u>List of tools</u>
Please install specified version or higher

 - fasterq-dump v2.11.3 from SRA-toolkit (https://hpc.nih.gov/apps/sratoolkit.html)
- FastQC v0.11.9  (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- MultiQC v1.11 (https://multiqc.info/)
- Trim Galore! v0.6.7 (https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
- Bowtie v1.3.1 (http://bowtie-bio.sourceforge.net/index.shtml)
- Blast+ (https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- MUSCLE v3.8.1551 (http://www.drive5.com/muscle)
- Python3

### - <u>Python libraries</u>

- pandas
- numpy
- matplotlib
- seaborn v0.12.0 or higher
- requests
- wget
- gtfparse
- shutil
- subprocess
- biopython

### - <u>Managing input data</u>
The following instructions explain how to add smallRNA-sequencing data to be processed or plant genomes to perform pmiRNAs search against.

#### <u>Adding plant genomes</u>

Default plant genomes used in the study are downloaded by Ensembl Plants FTP website: (https://plants.ensembl.org/info/data/ftp/index.html).

After chosing a new plant genome to add to the analysis, you should chose an alias for the species (e.g. *Arabidopsis thaliana* = ath).

<u>**Important Note**</u>: the .gtf file must contain the attributes *miRNA* or *pre-miRNA* in the *gene_biotype* attribute. Not all Plant genomes have this kind of annotation, so please check it out befor adding one to the script.

Create a file named <plant_alias>.txt containing the links to the chromosomes fasta files, and put it in the folder */data/genomes/* (plenty of examples can be found in the folder).

Next, create another file named <plant_alias.txt> containing the link to the .gtf genome annotation file, also found on the Ensembl Plants FTP website. Put the file in the folder */data/gtfs/* (again, plenty of examples can be found in the folder).

After doing this, edit the script *pmirfinder.sh* and add your new plant alias to the array named PLANTS, which can be found in the first few lines of code.

Run the script and all should be fine.

To remove a plant genome from the process, just delete the two files indicated above and the plant alias from the PLANTS array in */pmirfinder.sh*.

#### <u>Adding SRA projects</u> (sequencing data)

NCBI-SRA data used by default in this tool is organized in projects identified by the prefix PRJNA-

To add sequencing data create a file named *PRJNA<project_number>.txt* and put it in the folder */data/projects* (plenty of examples can be found in the folder).

After doing this, edit the script */pmirfinder.sh* and add *PRJNA<project_number>* to the array named PROJS, which can be found in the first few lines of code.

Run the script and all should be fine

To remove an SRA project from the process just delete the corresponding file in */data/projects/* and remove the corresponding PRJNA- entry from the PROJS array in */pmirfinder.sh*

If you wish to just remove a single experiment (i.e. an SRR-id), edit the corresponding PRJNA- file and delete the SRR-id.

### - <u>Understanding output</u>

pmiRFinder produces several files and graphs, which can be found in the /report/ folder at the end of the analysis. Following, a brief description of the output:
#### Text files
- *final_report.tsv* contains the list of the identified putative pmiRNAs, their average abundance espressed in Reads Per Million (RPM) and the matches they produced against different databases via BLAST, with the relative BitScore;
- *counts_matrix.tsv* contains a table listing all the sequencing experiments (SRRids) and all the putative pmiRNAs with their counts in RPM. It is used for clustering data;
- *pres_abs_matrix.tsv* just like *counts_matrix.tsv* but contains 1 if the pmiRNA is present or 0 if it is absent for each SRRid. Used for clustering data;
- *pmirnas.fasta* contains (guess what) the identified putative pmiRNAs in fasta format;
- *pmirnas_consensus.fasta* contains the consensus sequence calculated on putative pmiRNAs.
#### Graphs
- *alignments_barplot.png* and *alignments_barplot_norm.png* shows the average number of reads aligned to each plants genome. The second one is just normalized for genome size;
- *unclean_bitscores_distribution.png* and *clean_bitscores_distribution.png* display the Kernel Density Estimate (KDE) for the Bitscores reported by BLAST on different databases. The unclean one refers to putative pmiRNAs before they are filtered based on Bitscores comparison across different databases;
- *counts_clustering.png* and *pres_abs_clustering.png* show heatmaps and dendrograms for clustering performed on the relative matrices (see above in Text files). A big version is also provided, which shows every experiment and every sequence (BIG pictures, handle with care).
