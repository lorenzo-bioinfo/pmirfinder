from sys import argv
from os import listdir

basedir = argv[1]
proj = argv[2]

ls = listdir(f'{basedir}/{proj}/{proj}.human_mirs/')
files = [x for x in ls if x.endswith('txt')]
alignments = []
for file in files:
        with open(file) as f:
                for line in f:
                        if line.startswith('# reads with at least'):
                                alignments.append(float(line.strip().split(' ')[-1][1:-2]))
hgenome_meanscore = sum(alignments) / len(alignments)

if hgenome_meanscore >= 40:
        print("First alignment quality check passed: more than 40% (mean rate) of reads aligned to human miRBase")
else:
        print("##############################################################################################################\n")
        print("########## WARNING: Average alignment rate for human miRBase was < 40%. Consider removing {proj} from analysis\n")
        print("##############################################################################################################\n")
