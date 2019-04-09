
# snakefile for leafcutter pipeline

import pandas as pd
import os

# get variables out of config.yaml


# get stuff out of the config.yaml
leafcutterPath = config['leafcutterPath']
python2Path = config['python2Path']
python3Path = config['python3Path']

dataCode = config['dataCode']
inFolder = config['inFolder']
outFolder = config['outFolder']
stranded = config['stranded']

samples = config['samples']

refCondition = config['refCondition']
altCondition = config['altCondition']
# get sample information of support file
# using pandas

SAMPLES = pd.read_csv(samples, sep = '\t')['sample']

# extract sample IDs, BAM file paths

# will have to create the Leafcutter support file from this



# what does this do?
#localrules: all


# main rule - the final output of the pipeline should be the shiny RData file

# rule all:
# 	input:
# 		outFolder + '/shiny/{dataCode}_leafviz.RData'


rule all:
	#input: outFolder + "junctionList.txt"
	#input: outFolder + dataCode + "_perind_numers.counts.gz"
	input: outFolder + dataCode + "_ds_support.tsv"

# use regtools to extract junctions if not already completed
rule extractJunctions:
	input:
		inFolder + '{samples}.bam'
	output:
		outFolder + 'junctions/{samples}.junc'
	shell:
		"samtools index {input};"
		"regtools junctions extract -a 8 -m 50 -M 500000 -s {stranded} -o {output} {input}"


# Yang's script to cluster regtools junctions still uses python2
rule clusterJunctions:
	input: 
		expand(outFolder + 'junctions/{samples}.junc', samples = SAMPLES)
	output:
		junctionList=outFolder + "junctionList.txt",
		clusters=outFolder + dataCode + "_perind_numers.counts.gz",
		tempFiles = expand('{samples}.junc.{dataCode}.sorted.gz', samples = SAMPLES, dataCode = dataCode )
	shell:
		"touch {output.junctionList};"
		"for i in {input};"
		"do echo $i >> {output.junctionList};"
		"done;"
		"{python2Path} {leafcutterPath}/scripts/leafcutter_cluster_regtools.py -j {output.junctionList} -m 20 -o {outFolder}{dataCode} -l 500000;"

rule leafcutterDS:
	input:
		clusters=outFolder + dataCode + "_perind_numers.counts.gz",
		tempFiles=expand('{samples}.junc.{dataCode}.sorted.gz', samples = SAMPLES, dataCode = dataCode )
	output:
		support = outFolder + dataCode + "_ds_support.tsv"
	shell:
		"rm {input.tempFiles};"
		"Rscript ../sort_support.R "
		"	--samples {samples} "
		"	--dataCode {dataCode} "
		"	--refCondition {refCondition} "
		"	--altCondition {altCondition} "
		"	--outFolder {outFolder} "
		";"
		# "Rscript {leafcutterPath}/scripts/leafcutter_ds.R "
		# "	--output_prefix ${dataCode} "
		# "	--num_threads 4 "
	 #    "	--min_samples_per_intron 1 "
	 #    "	--min_samples_per_group 1 "
	 #    "	--min_coverage 5 "
		# "	{input.clusters} "
		# "	{input.support}"

# rule clusterIntrons:
# 	input:
# 		outFolder + "junctionList.txt"
# 	output:
# 		outFolder