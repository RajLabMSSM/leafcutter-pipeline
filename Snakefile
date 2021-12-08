# snakefile for leafcutter pipeline
# dependencies:
# samtools
# regtools
# python - pandas
# R - leafcutter
# R - a bunch of packages
from yaml import dump
import itertools
import pandas as pd
import os
import socket
import shutil 
# get variables out of config.yaml


leafcutterPath = config['leafcutterPath']
python2Path = config['python2Path']
python3Path = config['python3Path']
dataCode = config['dataCode']

print(" * Leafcutter Pipeline *")
print(" author: Jack Humphrey ")
print(" Dataset: %s " % dataCode)

inFolder = config['inFolder']
# create outFolder path using dataCode
outFolder = config['outFolder'] + config['dataCode'] + '/'
stranded = config['stranded']

if stranded != 0:
    strandParam = "--strand"
else:
    strandParam = ""

print( "output folder: ", outFolder)

metadataFile = config['metadata']
bamSuffix = config['bamSuffix']
# default is '.junc'; for samples processed with RAPiD use filtered junctions : '.Aligned.Quality.Sorted.bam.junc' 
juncSuffix = config['juncSuffix']

# get sample information of support file
# using pandas
metadata = pd.read_csv(metadataFile, sep = '\t')
samples = metadata['sample']

refCondition = config['refCondition']
altCondition = config['altCondition']
contrastSep = config['contrastSep']

# test whether all conditions are present in metadata
conditions = list(metadata['condition'])
overlaps = [i in conditions for i in refCondition + altCondition ]
if not all(overlaps):
    print(" * ERROR: some conditions in config.yaml are missing from the metadata ") 
    print(metadata)
    print("Expected conditions: ", conditions)
    exit()

# create all combinations of ref and alt.
# will be stored as two lists
# get product of two lists:
all_contrasts = list(itertools.product(refCondition,altCondition))
all_contrasts = [contrastSep.join(i) for i in all_contrasts]

print(" * %s contrasts detected: " % len(all_contrasts) )
print(" * %s " % all_contrasts ) 

#print(expand('junctions/{samples}{junc}', samples = samples, junc = juncSuffix))



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
# FDR threshold of shiny app
FDR_limit = leafcutterOpt["FDR_limit"]

# QC options
junction_qc = config["junctionQC"]
if junction_qc == True:
    print(" * extra junction QC selected")
    input_clusters = outFolder + dataCode + "_filtered_perind_numers.counts.gz"
else:
    input_clusters = outFolder + dataCode + "_perind_numers.counts.gz"

missingness = leafcutterOpt["missingness"]
maxsize = leafcutterOpt["maxsize"]
minratio = leafcutterOpt["minratio"]

# ds options
samplesPerIntron = leafcutterOpt["samplesPerIntron"]
samplesPerGroup = leafcutterOpt["samplesPerGroup"]
minCoverage = leafcutterOpt["minCoverage"]


# Chimera specific options
isChimera = "hpc.mssm.edu" in socket.getfqdn()

# not sure if this works when running in serial on interactive node
if isChimera:
    shell.prefix('export PS1="";source activate snakemake;ml R/3.6.0;')
#else:
#shell.prefix('conda activate leafcutterpipeline;')

# default R is now 3.6 - doesn't support leafcutter yet
#shell.prefix('ml R/3.6.0;')

clusterRegtools = config["clusterRegtools"]

if clusterRegtools == True:
    clusterScript = python3Path + " scripts/leafcutter_cluster_regtools.py"
    junctionMode = "regtools"
    #juncSuffix = ".junc" # enforce this
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
        refFolder + refCode + "/5UTRs.bed",  
        expand(outFolder + "{contrast}/" + dataCode + "_{contrast}_shiny.RData", contrast = all_contrasts),
        outFolder + "config.yaml",
        expand(outFolder + "{contrast}/" + dataCode + "_{contrast}_cassette_inclusion.tsv", contrast = all_contrasts),
        expand(outFolder + "{contrast}/" + dataCode + "_{contrast}_residual.counts.gz", contrast = all_contrasts)
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
        "samtools index {input};" # redundant if indexes are present
        "regtools junctions extract -a 8 -i 50 -I 500000 -s {stranded} -o {output} {input.bam}"
        # conda version of regtools uses i and I instead of m and M 
        # "ml regtools/0.5.1; "
        # "regtools junctions extract -a 8 -m 50 -M 500000 -s {stranded} -o {output} {input.bam}"


# remove weird contigs that cause add_chr() to break by adding "chr" to normal chr names
# now deprecated
rule stripContigs:
    input:
        'junctions/' + '{samples}' + juncSuffix
    output:
        'junctions/' + '{samples}' + juncSuffix + '.nocontigs'
    shell:
        "cat {input} |"
        "awk \'$1 ~/^chr/\' "
        " > {output} "

# copy the config and samples files in to the outFolder for posterity
rule copyConfig:
    input:
        metadata = metadataFile
    output: 
        config_out = outFolder + "config.yaml",
        metadata = outFolder + "samples.tsv"
    run:
        stream = open(output.config_out, 'w')
        dump(config, stream, default_flow_style = False)
        shutil.copy(input.metadata, output.metadata) 

# write junction list to file in python rather than bash
# bash has limit on length of shell command - 1000 samples in an array doesn't work.
rule writeJunctionList:
    input:
        juncFiles = expand('junctions/{samples}{junc}', samples = samples, junc = juncSuffix)
    output:
        junctionList = outFolder + "junctionList.txt",
        tempFileList = outFolder + "tempFileList.txt"
    params:
        tempFiles = expand('{samples}{junc}.{dataCode}.sorted.gz', samples = samples, junc = juncSuffix, dataCode = dataCode )
    run:
        # write juncFiles to junctionList, tempFiles to tempFiles list
        with open(output.junctionList, 'w') as f:
            for item in input.juncFiles:
                f.write("%s\n" % item)
        with open(output.tempFileList, 'w') as f:
            for item in params.tempFiles:
                f.write("%s\n" % item)
# Yang's script to cluster regtools junctions still uses python2
# I took an updated version from a github fork and fixed the bugs
# for samples processed with RAPiD just use the junctions from that and run the classic leafcutter_cluster.py
rule clusterJunctions:
    input: 
        junctionList = outFolder + "junctionList.txt",
        tempFileList = outFolder + "tempFileList.txt"
    output:
        clusters = outFolder + dataCode + "_perind_numers.counts.gz"
    params:
        #script = "scripts/leafcutter_cluster_regtools.py"
        script = clusterScript,
        strand = strandParam
    shell:
        #'touch {output.junctionList};'
        #'for i in {input};'
        #'do echo $i >> {output.junctionList};'
        #'done;'
        # from https://github.com/mdshw5/leafcutter/blob/master/scripts/leafcutter_cluster_regtools.py
        # now lives inside the leafcutter pipeline repo 
        '{params.script} '
        '-j {input.junctionList} --minclureads {minCluReads} '
        ' {params.strand} '
        ' --checkchrom '
        '--mincluratio {minCluRatio}  -o {outFolder}{dataCode} -l {intronMax};'
        # remove temporary files
        'for i in $(cat {input.tempFileList}); do rm $i; done'

## Perform junction and cluster level filtering before differential splicing analysis
## This is largely inspired by Balliu et al 2019
# maxsize - how many junctions can be in a cluster -default = 10
# missingness - what proportion of the samples can a junction be missing from? default = 0.1
# minratio - the minimum mean ratio a junction should contribute to a cluster. default = 0.05
rule junctionQC:
    input:
        clusters = outFolder + dataCode + "_perind_numers.counts.gz"
    output:
        outFolder + dataCode + "_filtered_perind_numers.counts.gz"
    params:
        script = "scripts/cluster_QC.R"
    shell:
        "ml R/3.6.0; "
        "Rscript {params.script} "
        "--outFolder {outFolder} "
        "--dataCode {dataCode} "
        "--missingness {missingness} "
        "--minratio {minratio} "
        "--maxsize {maxsize} "


# prepare support table the way leafcutter likes it
# be wary of sample names - sometimes the juncSuffix will be appended
rule prepareMeta:
    input: metadataFile
    output: support = outFolder + "{contrast}/" + dataCode + "_{contrast}_ds_support.tsv"
    params: js = juncSuffix
    shell:
        "Rscript scripts/sort_support.R "
        "   --metadata {input} "
        "   --dataCode {dataCode} "
        "   --contrast {wildcards.contrast} "
        "   --contrastSep {contrastSep} "
        "   --outFolder {outFolder} "
        "   --junctionMode {junctionMode} "
        "       --juncSuffix \"{params.js}\" "
        ";"


# run differential splicing
rule leafcutterDS:
    input:
        support = outFolder + "{contrast}/" + dataCode + "_{contrast}_ds_support.tsv",
        clusters = input_clusters #outFolder + dataCode + "_filtered_perind_numers.counts.gz"
    output:
        sigClusters = outFolder + "{contrast}/" + dataCode + "_{contrast}_cluster_significance.txt",
        effectSizes = outFolder + "{contrast}/" + dataCode + "_{contrast}_effect_sizes.txt"
    params:
        n_threads = leafcutterOpt['n_threads']
    shell:  
        'ml R/3.6.0; ' 
        'Rscript {leafcutterPath}/scripts/leafcutter_ds.R '
        '   --output_prefix {outFolder}{wildcards.contrast}/{dataCode}_{wildcards.contrast} '
        '   --num_threads {params.n_threads} '
        '   --min_samples_per_intron {samplesPerIntron} '
        '   --min_samples_per_group {samplesPerGroup} '
        '   --min_coverage {minCoverage} '
        '   {input.clusters} '
        '   {input.support} '

# for the Shiny app 
rule createRefs:
    input:
        refFolder + refFile
    output:
        refFolder + refCode + "_all_exons.txt.gz"
    shell:
        "perl {leafcutterPath}/leafviz/gtf2leafcutter.pl {input} -o {refFolder}{refCode}"

# for interpreting 2 junction clusters
# get first and last exons for each transcript
# first step converts GTF to genepred format
# second extracts first and last "full coding" exons from transcript
# this doesn't treat UTR sequence as exon which is silly
# alternative script produces bed files of 5 and 3 UTRs along with other files
rule getTerminalExons:
    input:
        refFolder + refFile
    params:
        gtftogenepred= "gtfToGenePred",
        genepredtobed = "scripts/genepred_to_bed.py",
        get_regions = "scripts/create_regions_from_gencode.R",
        outFolder = refFolder + refCode + "/"
    output:
        genepred = refFolder + refCode + ".genepred",
        starts = refFolder + refCode + "_first_exons.bed",
        ends = refFolder + refCode + "_last_exons.bed",
        utr5 = refFolder + refCode + "/5UTRs.bed" 
    shell:
        "{params.gtftogenepred} {input} {output.genepred};"
        "python {params.genepredtobed} --first_exon {output.genepred} > {output.starts} ; "
        "python {params.genepredtobed} --last_exon {output.genepred} > {output.ends}; "
        "mkdir {params.outFolder}; "
        "Rscript {params.get_regions} {input} {params.outFolder}"

# prepare results for shiny visualisation
rule prepareShiny:
    input:
        clusterCounts = input_clusters,
        sigClusters = outFolder + "{contrast}/" + dataCode + "_{contrast}_cluster_significance.txt",
        effectSizes = outFolder + "{contrast}/" + dataCode + "_{contrast}_effect_sizes.txt",
        support = outFolder + "{contrast}/" + dataCode + "_{contrast}_ds_support.tsv",
        exonFile = refFolder + refCode + "_all_exons.txt.gz"
    output:
        shinyData = outFolder + "{contrast}/" + dataCode + "_{contrast}_shiny.RData"
    shell:
        "Rscript {leafcutterPath}/leafviz/prepare_results.R "
        "{input.clusterCounts} {input.sigClusters} {input.effectSizes} "
        "{refFolder}{refCode} "
        "-o {output.shinyData} "
        "-m {input.support} "
        "-c {dataCode}_{wildcards.contrast} " 
        "--FDR {FDR_limit} "

# calculate delta PSI for each significant cluster
rule deltaPSI:
    input:
        app = outFolder + "{contrast}/" + dataCode + "_{contrast}_shiny.RData"
    output:
        outFolder + "{contrast}/" + dataCode + "_{contrast}_deltapsi_best.tsv",
        outFolder + "{contrast}/" + dataCode + "_{contrast}_deltapsi_full.tsv"
    params:
        script = "scripts/create_dPSI_table.R"
    shell:
        "Rscript {params.script} "
        "       --app {input.app} "
                "       --dataCode {dataCode} "
                "       --contrast {wildcards.contrast} "
                "       --contrastSep {contrastSep} "
                "       --outFolder {outFolder} "

rule classifyClusters:
    input: 
        app = outFolder + "{contrast}/" + dataCode + "_{contrast}_shiny.RData",
        psi_results = outFolder + "{contrast}/" + dataCode + "_{contrast}_deltapsi_full.tsv"
    output:
        outFolder + "{contrast}/" + dataCode + "_{contrast}_classifications.tsv",
        outFolder + "{contrast}/" + dataCode + "_{contrast}_summary.tsv",
        outFolder + "{contrast}/" + dataCode + "_{contrast}_cassette_inclusion.tsv"
    params:
        script = "scripts/classify_clusters.R"
    shell:
        "Rscript {params.script} "
        " -o {outFolder}{wildcards.contrast}/{dataCode}_{wildcards.contrast} "
                " {input.app} "

rule quantifyPSI:
    input:  
        support = outFolder + "{contrast}/" + dataCode + "_{contrast}_ds_support.tsv",
        counts = input_clusters
    output:
        outFolder + "{contrast}/" + dataCode + "_{contrast}_residual.counts.gz",
    params:
        script = "scripts/quantify_PSI.R",
        n_threads = leafcutterOpt['n_threads']
    shell:
        "Rscript {params.script} "
        " -c {input.support} "
        " -p {params.n_threads} "
        " -o {output} "
        " -j {juncSuffix} "
        " {input.counts} "
