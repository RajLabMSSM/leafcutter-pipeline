# Raj Lab Leafcutter Pipeline

Written by Jack Humphrey

## Installation

1. clone the repo


## Dependencies

Everything is managed through conda. The plan is to have a minerva-wide conda environment which will contain all the dependencies, both python and R packages.

Then you simply type:

`conda activate leafcutterPipeline`

Test this all works by 

```
cd example/
./snakejobs -n
```

-n specifies a dry run, so none of the scripts are executed.

## Pipeline

This pipeline takes a set of BAM files as input, along with some metadata.
For each bam, splice junction reads are extracted (extractJunctions).
All junctions are then clustered together (clusterJunctions).
Differential splicing is tested using Leafcutter (leafcutterDS).
An interactive shiny app is created (prepareShiny).
TODO: create an RMarkdown/plain text report file summarising the run and the results.

<p align="center">
  <img src="https://github.com/rajlabMSSM/LeafcutterPipeline/blob/master/dag.png">
</p>

## Input files

*metadata file* - this should be a tab-delimited text file. The first two columns, **sample** and **condition** are mandatory. Any other columns will be treated as covariates when fitting the differential splicing model.


*config.yaml* - this specifies the dataset-specific parameters.
For differential splicing analysis between two conditions, the reference condition **refCondition** and alternate condition **altCondition** must be set.

## Running on Chimera

`./snakejobs`

This is configured to run on Chimera using a script written by Brian FH from the Goate lab. Thanks Brian!
