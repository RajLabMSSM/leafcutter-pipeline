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
stranded = config['stranded']

print(outFolder)

metadata = config['metadata']
bamSuffix = config['bamSuffix']
# default is '.junc'; for samples processed with RAPiD use filtered junctions : '.Aligned.Quality.Sorted.bam.junc' 
juncSuffix = config['juncSuffix']

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

# Chimera specific options
isChimera = "hpc.mssm.edu" in socket.getfqdn()

# not sure if this works when running in serial on interactive node
#if isChimera:
#	shell.prefix('export PS1="";source activate leafcutter-pipeline;ml R;')
#else:
#shell.prefix('conda activate leafcutterpipeline;')

# default R is now 3.6 - doesn't support leafcutter yet
shell.prefix('ml R/3.5.3;')

clusterRegtools = config["clusterRegtools"]

if clusterRegtools == True:
        clusterScript = python3Path + " scripts/leafcutter_cluster_regtools.py"
	junctionMode = "regtools"
else:        
        clusterScript = python2Path + " scripts/leafcutter_cluster.py"
	junctionMode = "RAPiD"


localrules: copyConfig

rule all:
	#input: outFolder + "junctionList.txt"
	#input: outFolder + dataCode + "_perind_numers.counts.gz"
	#input: outFolder + dataCode + "_ds_support.tsv"
	input: 
		outFolder + dataCode + "_shiny.RData",
		outFolder + "config.yaml",
		outFolder + dataCode + "_deltapsi_best.tsv",
                outFolder + dataCode + "_classifications.tsv",
                outFolder + dataCode + "_summary.tsv",
		outFolder + dataCode + "_cassette_inclusion.tsv"

# index bams if needed
rule indexBams:
	input:
		bam = inFolder + '{samples}' + bamSuffix
	output:
		inFolder + '{samples}' + bamSuffix + ".bai"
	shell:
		"samtools index {input.bam}"

# use regtools to extract junctions if not already completed
rule extractJunctions:
	input:
		bam = inFolder + '{samples}' + bamSuffix,
		bai = inFolder + '{samples}' + bamSuffix + ".bai"
	output:
		'junctions/{samples}' + juncSuffix
	shell:
		#"samtools index {input};"	redundant if indexes are present
		#"regtools junctions extract -a 8 -m 50 -M 500000 -s {stranded} -o {output} {input}"
		# conda version of regtools uses i and I instead of m and M 
		"regtools junctions extract -a 8 -i 50 -I 500000 -s {stranded} -o {output} {input.bam}"

# copy the config and samples files in to the outFolder for posterity
rule copyConfig:
	input: 
		config = workflow.overwrite_configfile,
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
		expand('junctions/{samples}{junc}', samples = samples, junc = juncSuffix)
	output:
		junctionList = outFolder + "junctionList.txt",
		clusters = outFolder + dataCode + "_perind_numers.counts.gz"
	params:
		tempFiles = expand('{samples}{junc}.{dataCode}.sorted.gz', samples = samples, junc = juncSuffix, dataCode = dataCode ),
		#script = "scripts/leafcutter_cluster_regtools.py"
		script = clusterScript
	shell:
                'touch {output.junctionList};'
		'for i in {input};'
		'do echo $i >> {output.junctionList};'
		'done;'
		# from https://github.com/mdshw5/leafcutter/blob/master/scripts/leafcutter_cluster_regtools.py
		# now lives inside the leafcutter pipeline repo 
		'{params.script} '
		'-j {output.junctionList} --minclureads {minCluReads} '
		'--mincluratio {minCluRatio}  -o {outFolder}{dataCode} -l {intronMax};'
		'rm {params.tempFiles}'

# prepare support table the way leafcutter likes it
# be wary of sample names - sometimes the juncSuffix will be appended
rule prepareMeta:
	input: metadata
	output: support = outFolder + dataCode + "_ds_support.tsv"
	params: js = juncSuffix
	shell:
		"Rscript scripts/sort_support.R "
		"	--metadata {metadata} "
		"	--dataCode {dataCode} "
		"	--refCondition {refCondition} "
		"	--altCondition {altCondition} "
		"	--outFolder {outFolder} "
		"	--junctionMode {junctionMode} "
		"       --juncSuffix \"{params.js}\" "
		";"


# run differential splicing
rule leafcutterDS:
	input:
		support = outFolder + dataCode + "_ds_support.tsv",
		clusters = outFolder + dataCode + "_perind_numers.counts.gz"
	output:
		sigClusters = outFolder + dataCode + "_cluster_significance.txt",
		effectSizes = outFolder + dataCode + "_effect_sizes.txt"
	params:
		n_threads = leafcutterOpt['n_threads']
	shell:	
		'Rscript {leafcutterPath}/scripts/leafcutter_ds.R '
		'	--output_prefix {outFolder}{dataCode} '
		'	--num_threads {params.n_threads} '
		'	--min_samples_per_intron {samplesPerIntron} '
		'	--min_samples_per_group {samplesPerGroup} '
		'	--min_coverage {minCoverage} '
		'	{input.clusters} '
		'	{input.support} '

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
		"Rscript {leafcutterPath}/leafviz/prepare_results.R "
		"{input.clusterCounts} {input.sigClusters} {input.effectSizes} "
		"{refFolder}{refCode} "
		"-o {output.shinyData} "
		"-m {input.support} "
		"-c {dataCode}" 


rule deltaPSI:
	input:
		app = outFolder + dataCode + "_shiny.RData"
	output:
		outFolder + dataCode + "_deltapsi_best.tsv",
		outFolder + dataCode + "_deltapsi_full.tsv"
	params:
		script = "scripts/create_dPSI_table.R"
	shell:
		"Rscript {params.script} "
		"       --app {input.app} "
                "       --dataCode {dataCode} "
                "       --refCondition {refCondition} "
                "       --altCondition {altCondition} "
                "       --outFolder {outFolder} "

rule classifyClusters:
	input: 
		app = outFolder + dataCode + "_shiny.RData",
		psi_results = outFolder + dataCode + "_deltapsi_full.tsv"
	output:
		outFolder + dataCode + "_classifications.tsv",
		outFolder + dataCode + "_summary.tsv",
		outFolder + dataCode + "_cassette_inclusion.tsv"
	params:
		script = "scripts/classify_clusters.R"
	shell:
		"Rscript {params.script} "
		" -o {outFolder}/{dataCode} "
                " {input.app} "
