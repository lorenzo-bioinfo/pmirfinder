from sys import argv
from io import StringIO
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Align.Applications import MuscleCommandline

infile = argv[1]
outfile = argv[2]
thresh = float(argv[3])
muscle = MuscleCommandline(input=infile)
muscleout = muscle()[0]
alignment = AlignIO.read(StringIO(muscleout), "fasta")
alignment_summary = AlignInfo.SummaryInfo(alignment)
consensus = alignment_summary.gap_consensus(threshold = thresh, ambiguous = 'N')
with open(outfile, 'w') as f:
	f.write(f'>Consensus sequence of length {len(consensus)} for {infile}\n{consensus}')
