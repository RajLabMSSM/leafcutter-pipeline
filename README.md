# Raj Lab Leafcutter Pipeline

Written by Jack Humphrey

## Installation

1. clone the repo

2. install leafcutter

```
ml R/3.6.0
R
install.packages("remotes")
remotes::install_github("stan-dev/rstantools")
remotes::install_github("davidaknowles/leafcutter/leafcutter", ref = "stanfixagain")
```

3. clone the leafcutter repo to the cluster. A path to the directory will be needed for the `config.yaml`.


## Dependencies

TODO: give full conda recipe

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


### Starting with BAM files

You need to then have a folder where all the BAMs are kept. If the BAMs are split across many folders, then create a single folder with symbolic links to each BAM. The pipeline requires that the names of each BAM are in the format:
	`<sampleID><bamSuffix>`
You specify the `bamSuffix` in the *config.yaml*
So a sample with ID 'sample1' and a bam file 'sample1_aligned.bam' would have a `bamSuffix` of '_aligned.bam'

### Starting with junction files

If you already have junction files prepared, then copy or symlink them to a folder called /junctions within the pipeline directory.
If the junctions were created using regtools, then set leafcutter clusterRegtools: to 'True'. 
If the junctions were created by leafcutter's bam2junc.sh, as is the case for RAPiD, then keep it set to 'False'.


## Running on Chimera

You can run the pipeline either in serial as a single process, or in parallel using Chimera's job submission.
In practice, only *extractJunctions* requires parallel execution. The following steps can all be run together.

1. Running each step in serial - assume that junctions have already been created

	`./snakelocal`

2. Running in parallel with LSF job submission

	`./snakejobs`


This is configured to run on Chimera using a script written by Brian FH from the Goate lab. Thanks Brian!

In *cluster.yaml* the resource requirements for each step are laid out.

