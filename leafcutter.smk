# snakefile for leafcutter pipeline
# dependencies:
# samtools
# regtools
# python - pandas
# R - leafcutter
# R - a bunch of packages


import pandas as pd
import os

# get variables out of config.yaml


# get stuff out of the config.yaml
leafcutterPath = config['leafcutterPath']
#python2Path = config['python2Path']
python3Path = config['python3Path']

dataCode = config['dataCode']
inFolder = config['inFolder']
outFolder = config['outFolder']
stranded = config['stranded']

samples = config['samples']
bamSuffix = config['bamSuffix']

refCondition = config['refCondition']
altCondition = config['altCondition']

# annotation 
refFolder = config['refFolder']
refFile = config['refFile']
refCode = config['refCode']

# get sample information of support file
# using pandas
SAMPLES = pd.read_csv(samples, sep = '\t')['sample']

rule all:
	#input: outFolder + "junctionList.txt"
	#input: outFolder + dataCode + "_perind_numers.counts.gz"
	#input: outFolder + dataCode + "_ds_support.tsv"
	input: outFolder + dataCode + "_shiny.RData"

# use regtools to extract junctions if not already completed
rule extractJunctions:
	input:
		inFolder + '{samples}' + bamSuffix
	output:
		outFolder + 'junctions/{samples}.junc'
	shell:
		"samtools index {input};"
		#"regtools junctions extract -a 8 -m 50 -M 500000 -s {stranded} -o {output} {input}"
		# conda version of regtools uses i and I instead of m and M 
		"regtools junctions extract -a 8 -i 50 -I 500000 -s {stranded} -o {output} {input}"

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
		#"{python2Path} {leafcutterPath}/scripts/leafcutter_cluster_regtools.py -j {output.junctionList} -m 20 -o {outFolder}{dataCode} -l 500000;"
# python3 version - taken from a fork from https://github.com/mdshw5/leafcutter/blob/master/scripts/leafcutter_cluster_regtools.py
		"{python3Path} ../scripts/leafcutter_cluster_regtools.py -j {output.junctionList} -m 20 -o {outFolder}{dataCode} -l 500000;"


rule leafcutterDS:
	input:
		clusters=outFolder + dataCode + "_perind_numers.counts.gz",
		tempFiles=expand('{samples}.junc.{dataCode}.sorted.gz', samples = SAMPLES, dataCode = dataCode )
	output:
		support = outFolder + dataCode + "_ds_support.tsv",
		sigClusters = outFolder + dataCode + "_cluster_significance.txt",
		effectSizes = outFolder + dataCode + "_effect_sizes.txt"
	shell:
		"rm {input.tempFiles};"
		"Rscript ../scripts/sort_support.R "
		"	--samples {samples} "
		"	--dataCode {dataCode} "
		"	--refCondition {refCondition} "
		"	--altCondition {altCondition} "
		"	--outFolder {outFolder} "
		";"
		"Rscript {leafcutterPath}/scripts/leafcutter_ds.R "
		"	--output_prefix {outFolder}{dataCode} "
		"	--num_threads 1 "
		"	--min_samples_per_intron 1 "
		"	--min_samples_per_group 1 "
		"	--min_coverage 5 "
		"	{input.clusters} "
		"	{output.support} "

rule createRefs:
	input:
		refFolder + refFile
	output:
		refFolder + refCode + "_all_exons.txt.gz"
	shell:
		"perl {leafcutterPath}/leafviz/gtf2leafcutter.pl {input} -o {refFolder}{refCode}"

rule prepareShiny:
	input:
		clusterCounts = outFolder + dataCode + "_perind_numers.counts.gz",
		sigClusters = outFolder + dataCode + "_cluster_significance.txt",
		effectSizes = outFolder + dataCode + "_effect_sizes.txt",
		support = outFolder + dataCode + "_ds_support.tsv",
		exonFile = refFolder + refCode + "_all_exons.txt.gz"
	output:
		shinyData = outFolder + dataCode + "_shiny.RData"
	shell:
		"Rscript {leafcutterPath}/leafviz/prepare_results.R "
		"{input.clusterCounts} {input.sigClusters} {input.effectSizes} "
		"{refFolder}{refCode} "
		"-o {output.shinyData} "
		"-m {input.support} "
		"-c {dataCode}" 


