# permcutter - permutation testing for leafcutter

# number of permutations  of the condition column to perform differential splicing on
nPerm = 1000

# snakefile for leafcutter pipeline
# dependencies:
# samtools
# regtools
# python - pandas
# R - leafcutter
# R - a bunch of packages


import pandas as pd
import os
import socket
# get variables out of config.yaml

leafcutterPath = config['leafcutterPath']
python2Path = config['python2Path']
python3Path = config['python3Path']
dataCode = config['dataCode']

print(dataCode)

inFolder = config['inFolder']
# create outFolder path using dataCode
outFolder = config['outFolder'] + config['dataCode'] + '/'
permFolder = config['outFolder'] + config['dataCode'] + '/permutation/'
stranded = config['stranded']
juncSuffix = config['juncSuffix']

if stranded != 0:
	strandParam = "--strand"
else:
	strandParam = ""

print(outFolder)

metadata = config['metadata']

refCondition = config['refCondition']
altCondition = config['altCondition']

# leafcutter options
leafcutterOpt = config['leafcutter']
# ds options
samplesPerIntron = leafcutterOpt["samplesPerIntron"]
samplesPerGroup = leafcutterOpt["samplesPerGroup"]
minCoverage = leafcutterOpt["minCoverage"]

# get sample information of support file
# using pandas
samples = pd.read_csv(metadata, sep = '\t')['sample']

# Chimera specific options
isChimera = "hpc.mssm.edu" in socket.getfqdn()

# not sure if this works when running in serial on interactive node
if isChimera:
	shell.prefix('export PS1="";source activate leafcutter-pipeline;ml R/3.6.0;')
#else:
#shell.prefix('conda activate leafcutterpipeline;')

# default R is now 3.6 - doesn't support leafcutter yet
#shell.prefix('ml R/3.6.0;')

clusterRegtools = config["clusterRegtools"]

if clusterRegtools == True:
        clusterScript = python3Path + " scripts/leafcutter_cluster_regtools.py"
	junctionMode = "regtools"
	strandParam = "" # strandParam only needed for normal clustering
else:        
        clusterScript = python2Path + " scripts/leafcutter_cluster.py" 
	junctionMode = "RAPiD"


localrules: copyConfig, writeJunctionList


rule all:
	#input: outFolder + "junctionList.txt"
	#input: outFolder + dataCode + "_perind_numers.counts.gz"
	#input: outFolder + dataCode + "_ds_support.tsv"
	input: 
		outFolder + "config.yaml",
                permFolder + dataCode + "_all_permutation_results.tsv"
		#expand(permFolder + dataCode + "_permutation_{permutation}_cluster_significance.txt", permutation = range(1,nPerm+1) )
		#outFolder + dataCode + "_classifications.tsv",
                #outFolder + dataCode + "_summary.tsv",
		#outFolder + dataCode + "_cassette_inclusion.tsv",
		#outFolder + dataCode + "_residual.counts.gz"

# prepare support table the way leafcutter likes it
# be wary of sample names - sometimes the juncSuffix will be appended
rule permuteMeta:
	input: metadata
	output: support = permFolder + dataCode + "_permutation_{permutation}_ds_support.tsv"
	params: js = juncSuffix,
		permDataCode = dataCode + "_permutation_{permutation}"
	shell:
		"Rscript scripts/sort_support.R "
		"	--metadata {metadata} "
		"	--dataCode {params.permDataCode} "
		"	--refCondition {refCondition} "
		"	--altCondition {altCondition} "
		"	--outFolder {permFolder} "
		"	--junctionMode {junctionMode} "
		"	--permutation "
		"       --juncSuffix \"{params.js}\" "
		";"


# run differential splicing
rule leafcutterDS:
	input:
		support = permFolder + dataCode + "_permutation_{permutation}_ds_support.tsv",
		clusters = outFolder + dataCode + "_filtered_perind_numers.counts.gz"
	output:
		sigClusters = permFolder + dataCode + "_permutation_{permutation}_cluster_significance.txt",
		effectSizes = permFolder + dataCode + "_permutation_{permutation}_effect_sizes.txt"
	params:
		perm = "{permutation}",
		n_threads = leafcutterOpt['n_threads']
	shell:	
		'Rscript {leafcutterPath}/scripts/leafcutter_ds.R '
		'	--output_prefix {permFolder}{dataCode}_permutation_{params.perm} '
		'	--num_threads {params.n_threads} '
		'	--min_samples_per_intron {samplesPerIntron} '
		'	--min_samples_per_group {samplesPerGroup} '
		'	--min_coverage {minCoverage} '
		'	{input.clusters} '
		'	{input.support} '

rule combinePermutationResults:
	input:
 		expand(permFolder + dataCode + "_permutation_{permutation}_cluster_significance.txt", permutation = range(1,nPerm+1) )
	output: 
		permFolder + dataCode + "_all_permutation_results.tsv"
	params:
		script = "scripts/combine_permutation_results.R"
	shell:
		"ml R/3.6.0;"
		"Rscript {params.script} "
		" --outFolder {permFolder} "
		" --dataCode {dataCode} "
