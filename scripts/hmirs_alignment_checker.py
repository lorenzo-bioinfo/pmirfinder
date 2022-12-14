from sys import argv
from os import listdir

basedir = argv[1]
proj = argv[2]

ls = listdir(f'{basedir}/{proj}/{proj}.human_mirs/')
files = [x for x in ls if x.endswith('txt')]
alignments = []
for file in files:
        with open(f'{basedir}/{proj}/{proj}.human_mirs/{file}', 'r') as f:
                for line in f:
                        if line.startswith('# reads with at least'):
                                alignments.append(float(line.strip().split(' ')[-1][1:-2]))
hmirs_meanscore = sum(alignments) / len(alignments)

if hmirs_meanscore >= 40:
        print(f"First alignment quality check passed: more than 40% (mean rate) of reads aligned to human miRBase ({hmirs_meanscore:.2f}%)")
else:
        print("##############################################################################################################\n")
        print(f"########## WARNING: Average alignment rate for human miRBase was < 40% ({hmirs_meanscore:.2f}%). Consider removing {proj} from analysis\n")
        print("##############################################################################################################\n")
