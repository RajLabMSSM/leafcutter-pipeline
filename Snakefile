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

leafcutterPath = config['leafcutterPath']
#python2Path = config['python2Path']
python3Path = config['python3Path']

dataCode = config['dataCode']

print(dataCode)

inFolder = config['inFolder']
# create outFolder path using dataCode
outFolder = config['outFolder'] + config['dataCode'] + '/'
stranded = config['stranded']

print(outFolder)

metadata = config['metadata']
bamSuffix = config['bamSuffix']

refCondition = config['refCondition']
altCondition = config['altCondition']

# annotation 
refFolder = config['refFolder']
refFile = config['refFile']
refCode = config['refCode']

# leafcutter options
leafcutterOpt = config['leafcutter']
# clustering options
minCluRatio = leafcutterOpt["minCluRatio"]
minCluReads = leafcutterOpt["minCluReads"]
intronMax = leafcutterOpt["intronMax"]
# ds options
samplesPerIntron = leafcutterOpt["samplesPerIntron"]
samplesPerGroup = leafcutterOpt["samplesPerGroup"]
minCoverage = leafcutterOpt["minCoverage"]

# get sample information of support file
# using pandas
samples = pd.read_csv(metadata, sep = '\t')['sample']

localrules: copyConfig

rule all:
	#input: outFolder + "junctionList.txt"
	#input: outFolder + dataCode + "_perind_numers.counts.gz"
	#input: outFolder + dataCode + "_ds_support.tsv"
	input: 
		outFolder + dataCode + "_shiny.RData",
		outFolder + "config.yaml"

# use regtools to extract junctions if not already completed
rule extractJunctions:
	input:
		inFolder + '{samples}' + bamSuffix
	output:
		'junctions/{samples}.junc'
	shell:
		'export PS1="";'
		"source activate leafcutterPipeline;"
		"samtools index {input};"
		#"regtools junctions extract -a 8 -m 50 -M 500000 -s {stranded} -o {output} {input}"
		# conda version of regtools uses i and I instead of m and M 
		"regtools junctions extract -a 8 -i 50 -I 500000 -s {stranded} -o {output} {input}"

# copy the config and samples files in to the outFolder for posterity
rule copyConfig:
	input: 
		config = "config.yaml",
		metadata = metadata

	output: 
		config = outFolder + "config.yaml",
		metadata = outFolder + "samples.tsv"
	shell:
		"cp {input.config} {output.config};"
		"cp {input.metadata} {output.metadata}"


# Yang's script to cluster regtools junctions still uses python2
# I took an updated version from a github fork and fixed the bugs
rule clusterJunctions:
	input: 
		expand('junctions/{samples}.junc', samples = samples)
	output:
		junctionList = outFolder + "junctionList.txt",
		clusters = outFolder + dataCode + "_perind_numers.counts.gz",
		#tempFiles = expand('{samples}.junc.{dataCode}.sorted.gz', samples = samples, dataCode = dataCode )
	shell:
                'export PS1="";'
                'source activate leafcutterPipeline;'
		'touch {output.junctionList};'
		'for i in {input};'
		'do echo $i >> {output.junctionList};'
		'done;'
		# from https://github.com/mdshw5/leafcutter/blob/master/scripts/leafcutter_cluster_regtools.py
		# now lives inside the leafcutter pipeline repo 
		'{python3Path} ../scripts/leafcutter_cluster_regtools.py '
		'-j {output.junctionList} --minclureads {minCluReads} '
		'--mincluratio {minCluRatio}  -o {outFolder}{dataCode} -l {intronMax};'
		#'rm {output.tempFiles}'

# run differential splicing
rule leafcutterDS:
	input:
		clusters=outFolder + dataCode + "_perind_numers.counts.gz",
		#tempFiles=expand('{samples}.junc.{dataCode}.sorted.gz', samples = samples, dataCode = dataCode )
	output:
		support = outFolder + dataCode + "_ds_support.tsv",
		sigClusters = outFolder + dataCode + "_cluster_significance.txt",
		effectSizes = outFolder + dataCode + "_effect_sizes.txt"
	shell:	
		"ml R;"
		#"rm {input.tempFiles};"
		"Rscript ../scripts/sort_support.R "
		"	--metadata {metadata} "
		"	--dataCode {dataCode} "
		"	--refCondition {refCondition} "
		"	--altCondition {altCondition} "
		"	--outFolder {outFolder} "
		";"
		'Rscript {leafcutterPath}/scripts/leafcutter_ds.R '
		'	--output_prefix {outFolder}{dataCode} '
		'	--num_threads 1 '
		'	--min_samples_per_intron {samplesPerIntron} '
		'	--min_samples_per_group {samplesPerGroup} '
		'	--min_coverage {minCoverage} '
		'	{input.clusters} '
		'	{output.support} '

# for the Shiny app 
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
		"ml R;"
		"Rscript {leafcutterPath}/leafviz/prepare_results.R "
		"{input.clusterCounts} {input.sigClusters} {input.effectSizes} "
		"{refFolder}{refCode} "
		"-o {output.shinyData} "
		"-m {input.support} "
		"-c {dataCode}" 


