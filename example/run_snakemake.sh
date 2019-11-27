#conda activate default
# samtools and regtools are installed on this default environment
#ml R
# this version of R has the current version of leafcutter

snakemake  -s ../Snakefile --configfile config.yaml
