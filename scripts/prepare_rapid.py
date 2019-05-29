
import pandas as pd
import os
import sys

metadata = sys.argv[1]


exists = os.path.isfile(metadata)
if not exists:
    print( metadata , "does not exist!" )
    exit()

print("reading", metadata)
meta = pd.read_csv(metadata, sep = '\t')

samples = meta['sample']
rapid_path = meta['rapid_path']


bam_dir = os.path.abspath("data/")
os.makedirs(bam_dir, exist_ok = True)

for i in range(len(samples)):
	print("sample:", samples[i])
	bam_source = os.path.abspath(rapid_path[i] + "/Processed/RAPiD/bams/" + samples[i] + ".bam")
	bam_target = bam_dir + "/" + samples[i] + ".bam"

	if not os.path.exists(bam_target):
		print(bam_target)
		os.symlink(bam_source, bam_target )

