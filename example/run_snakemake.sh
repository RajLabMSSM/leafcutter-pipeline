conda activate default
ml R

snakemake -n -s ../leafcutter.smk --configfile config.yaml
