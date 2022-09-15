#!/usr/bin/bash

PLANTS=('ath' 'bevul' 'brapa' 'daucar' 'gmax' 'osat' 'phavul' 'soltub' 'solyc' 'vitvin')
PROJ=('PRJNA593788' 'PRJNA637898' 'PRJNA735638' 'PRJNA857648' 'PRJNA858789')
BASEDIR=$(pwd)
NUMTHREADS=$(grep -c ^processor /proc/cpuinfo)

read -p "Do you wish to download and index reference genomes (This will take a lot of time)?[y/n]" download
download=${download:-y}
if [ $download == 'y' ]
then
	python $BASEDIR/scripts/genomes_downloader.py $BASEDIR
	echo "Building indexes for human filter data"
	mkdir $BASEDIR/bt-index/hmirbase
	cp $BASEDIR/data/mirbase/hmirnas.fa $BASEDIR/bt-index/hmirbase/
	echo "Building Human Genome index (This will take some time)..."
	bowtie-build --quiet --threads $NUMTHREADS -f $BASEDIR/bt-index/hgenome/hgenome.fa $BASEDIR/bt-index/hgenome/hgenome
	echo "Building Human cDNA index (This will take some time)..."
	bowtie-build --quiet --threads $NUMTHREADS -f $BASEDIR/bt-index/hcdna/hcdna.fa $BASEDIR/bt-index/hcdna/hcdna
	echo "Building Human non-chromosomal dna index..."
	bowtie-build --quiet --threads $NUMTHREADS -f $BASEDIR/bt-index/hnonchr/hnonchr.fa $BASEDIR/bt-index/hnonchr/hnonchr
	echo "Building Human miRBase index..."
	bowtie-build --quiet --threads $NUMTHREADS -f $BASEDIR/bt-index/hmirbase/hmirnas.fa $BASEDIR/bt-index/hmirbase/hmirbase
	echo "Building indexes for plants genomes"
	for plant in ${PLANTS[@]}
	do
		echo "Building $plant genome index (likely to take some time)..."
		bowtie-build --quiet --threads $NUMTHREADS -f "$BASEDIR/bt-index/$plant/$plant.fa" "$BASEDIR/bt-index/$plant/$plant"
	done
fi
#############################################################################################################
for project in ${PROJ[@]}
do
	mkdir $BASEDIR/$project

	echo "Starting analysis of $project"
	mkdir "$BASEDIR/$project/$project.raw"
	cat "$BASEDIR/data/projects/$project.txt" | while read line
	do
		echo "Started downloading $line from $project"
		fasterq-dump --threads $NUMTHREADS -O $BASEDIR/$project/$project.raw/  --skip-technical $line
	done
###########################################################################################################
	echo "Startig quality check for raw reads of $project"
	mkdir "$BASEDIR/$project/$project.rawQC"
	fastqc -t $NUMTHREADS -o $BASEDIR/$project/$project.rawQC $BASEDIR/$project/$project.raw/*.fastq
	echo "Summarizing raw QC results from $project"
	multiqc $BASEDIR/$project/$project.rawQC -o $BASEDIR/$project/$project.rawQC
############################################################################################################
	echo "Starting trimming for raw reads of $project"
	mkdir $BASEDIR/$project/$project.trimmed
	cat $BASEDIR/data/projects/$project.txt | while read line
	do
		echo "Trimming reads from $line in $project..."
		trim_galore --phred33 --max_length 26 -j 8 -o "$BASEDIR/$project/$project.trimmed" "$BASEDIR/$project/$project.raw/$line.fastq" > "$BASEDIR/$project/$project.trimmed/$line.trigalore_output.txt" 2> /dev/null
		echo "Finished trimming for $line, removing raw reads file..."
		rm $BASEDIR/$project/$project.raw/$line.fastq
	done
	echo "Trimming completed for all files, removing raw reads folder"
	rmdir $BASEDIR/$project/$project.raw
################################################################################################################
	echo "Startig quality check for trimmed reads of $project"
	mkdir $BASEDIR/$project/$project.trimmedQC
	fastqc -t $NUMTHREADS -o $BASEDIR/$project/$project.trimmedQC $BASEDIR/$project/$project.trimmed/*.fq
	echo "Summarizing trimmed QC results from $project"
	multiqc $BASEDIR/$project/$project.trimmedQC -o $BASEDIR/$project/$project.trimmedQC
################################################################################################################
	mkdir $BASEDIR/$project/$project.non_human_mirs
	mkdir $BASEDIR/$project/$project.human_mirs
	cat $BASEDIR/data/projects/$project.txt | while read line
	do
		echo "Filtering reads on human miRBase for $line in $project"
		touch $BASEDIR/$project/$project.human_mirs/$line.alignment_report.txt
		bowtie -x $BASEDIR/bt-index/hmirbase/hmirbase -n 0 -l 20 --norc --best --strata -m 1 -p $NUMTHREADS --un $BASEDIR/$project/$project.non_human_mirs/$line.nonhmirs.fastq --al $BASEDIR/$project/$project.human_mirs/$line.hmirs.fastq -S $BASEDIR/$project/$project.trimmed/${line}_trimmed.fq > $BASEDIR/$project/$project.human_mirs/$line.sam 2>> $BASEDIR/$project/$project.human_mirs/$line.alignment_report.txt
		echo "Removing SAM file and trimmed reads for ${line}"
		rm $BASEDIR/$project/$project.human_mirs/$line.sam
		rm $BASEDIR/$project/$project.trimmed/${line}_trimmed.fq
	done
	echo "Running first alignment check..."
	python $BASEDIR/scripts/hmirs_alignment_checker.py $BASEDIR $project
########################################################################################################################
	mkdir $BASEDIR/$project/$project.non_human_chr
	mkdir $BASEDIR/$project/$project.human_chr
	cat $BASEDIR/data/projects/$project.txt | while read line
	do
		echo "Filtering reads on human genome (Grch38) for $line in $project"
		touch $BASEDIR/$project/$project.human_chr/$line.alignment_report.txt
		bowtie -x $BASEDIR/bt-index/hgenome/hgenome -n 0 -l 20 --norc --best --strata -m 1 -p $NUMTHREADS --un $BASEDIR/$project/$project.non_human_chr/$line.nonhchr.fastq --al $BASEDIR/$project/$project.human_chr/$line.hchr.fastq -S $BASEDIR/$project/$project.non_human_mirs/$line.nonhmirs.fastq > $BASEDIR/$project/$project.human_chr/$line.sam 2>> $BASEDIR/$project/$project.human_chr/$line.alignment_report.txt
		echo "Removing SAM file and trimmed reads for ${line}"
		rm $BASEDIR/$project/$project.human_chr/$line.sam
		rm $BASEDIR/$project/$project.non_human_mirs/$line.nonhmirs.fastq
	done
	echo "Running second and last alignment check"
	python $BASEDIR/scripts/hgenome_alignment_checker.py $BASEDIR $project
###########################################################################################################################
	mkdir $BASEDIR/$project/$project.non_human_extra
	mkdir $BASEDIR/$project/$project.human_extrachr
	cat $BASEDIR/data/projects/$project.txt | while read line
	do
		echo Filtering reads on human extra-chromosomal regions for $line in $project
		touch $BASEDIR/$project/$project.human_extrachr/$line.alignment_report.txt
		bowtie -x $BASEDIR/bt-index/hnonchr/hnonchr -n 0 -l 20 --norc --best --strata -m 1 -p $NUMTHREADS --un $BASEDIR/$project/$project.non_human_extra/$line.non_human_extra.fastq --al $BASEDIR/$project/$project.human_extrachr/$line.hextrachr.fastq -S $BASEDIR/$project/$project.non_human_chr/$line.nonhchr.fastq > $BASEDIR/$project/$project.human_extrachr/$line.sam 2>> $BASEDIR/$project/$project.human_extrachr/$line.alignment_report.txt
		echo "Removing SAM file and trimmed reads for ${line}"
		rm $BASEDIR/$project/$project.human_extrachr/$line.sam
		rm $BASEDIR/$project/$project.non_human_chr/$line.nonhchr.fastq
	done
##########################################################################################################################
	mkdir $BASEDIR/$project/$project.non_human
	mkdir $BASEDIR/$project/$project.human_cdna
	cat $BASEDIR/data/projects/$project.txt | while read line
	do
		echo "Filtering reads on human cDNA for $line in $project"
		touch $BASEDIR/$project/$project.human_cdna/$line.alignment_report.txt
		bowtie -x $BASEDIR/bt-index/hcdna/hcdna -n 0 -l 20 --norc --best --strata -m 1 -p $NUMTHREADS --un $BASEDIR/$project/$project.non_human/$line.non_human.fastq --al $BASEDIR/$project/$project.human_cdna/$line.hcdna.fastq -S $BASEDIR/$project/$project.non_human_extra/$line.non_human_extra.fastq > $BASEDIR/$project/$project.human_cdna/$line.sam 2>> $BASEDIR/$project/$project.human_cdna/$line.alignment_report.txt
		echo "Removing SAM file and trimmed reads for ${line}"
		rm $BASEDIR/$project/$project.human_cdna/$line.sam
		rm $BASEDIR/$project/$project.non_human_extra/$line.non_human_extra.fastq
	done
##########################################################################################################################
	echo "Performing QC on remaining non-human reads"
	mkdir $BASEDIR/$project/nonhumanQC/
	fastqc -t $NUMTHREADS -o $BASEDIR/$project/nonhumanQC/ $BASEDIR/$project/$project.non_human/*.fastq
	echo "Summarizing QC results for non-human reads"
	multiqc $BASEDIR/$project/nonhumanQC/ -o $BASEDIR/$project/nonhumanQC/
###########################################################################################################################
	echo "Starting alignment of non-human reads to plants genomes for $project"
	for plant in ${PLANTS[@]}
	do
		mkdir $BASEDIR/$project/$project.$plant
		cat $BASEDIR/data/projects/$project.txt | while read line
		do
			echo "Aligning non-human reads from $line from $project to $plant genome..."
			bowtie -x $BASEDIR/bt-index/$plant/$plant -n 0 -l 20 --norc --best --strata -m 1 -p 24 --un $BASEDIR/$project/$project.$plant/$line.unaligned.fastq -S $BASEDIR/$project/$project.non_human/$line.non_human.fastq > $BASEDIR/$project/$project.$plant/$line.$plant.sam 2>> $BASEDIR/$project/$project.$plant/$line.alignment_report.txt
		done
		echo "Extracting $plant miRNAs counts for $project"
		python $BASEDIR/scripts/mircounter.py $plant $project $BASEDIR 2> /dev/null
		echo "Summarizing results for $project"
		python $BASEDIR/scripts/fastparser.py $plant $project $BASEDIR
		echo "Collecting and summarizing results for $plant in $project"
		python $BASEDIR/scripts/sumproject.py $plant $project $BASEDIR
	done
done
###########################################################################################################################
for plant in ${PLANTS[@]}
do
	python $BASEDIR/scripts/sumdata.py $plant $BASEDIR
done
echo "Retrieving unique sequences"
python $BASEDIR/scripts/unique.py $BASEDIR
echo "Reporting sequences in fasta format"
python $BASEDIR/scripts/final_to_fasta.py $BASEDIR
echo "Building blast databases for sequences identification"
mkdir $BASEDIR/blast_db/mirbase
cp $BASEDIR/data/mirbase/mature.fa $BASEDIR/blast_db/mirbase
makeblastdb -in $BASEDIR/blast_db/mirbase/mature.fa -dbtype nucl -out $BASEDIR/blast_db/mirbase/mirbase
makeblastdb -in $BASEDIR/blast_db/hrefseq/hrefseq.fna -dbtype nucl -out $BASEDIR/blast_db/hrefseq/hrefseq
blastn -task blastn-short -num_threads $NUMTHREADS -outfmt 6 -db $BASEDIR/blast_db/mirbase/mirbase -query "${BASEDIR}/results/final_summary.fasta" -out $BASEDIR/results/final_blastres.tsv > $BASEDIR/blast_db/mirbase_blast_report.txt
blastn -task blastn-short -num_threads $NUMTHREADS -outfmt 6 -db $BASEDIR/blast_db/hrefseq/hrefseq -query "${BASEDIR}/results/final_summary.fasta" -out $BASEDIR/results/refseq_blast_final.tsv > $BASEDIR/blast_db/refseq_blast_report.txt
echo "Parsing BLAST results"
python $BASEDIR/scripts/final_blastparser.py $BASEDIR
python $BASEDIR/scripts/refseq_blastparser.py $BASEDIR
echo 'Retrieving full list of SRRids'
cat $BASEDIR/data/projects/*.txt > $BASEDIR/data/allids.txt
mkdir $BASEDIR/results/
mkdir $BASEDIR/results/individuals/
python $BASEDIR/scripts/individuals_seqs.py $BASEDIR
python $BASEDIR/scripts/individuals_count.py $BASEDIR
