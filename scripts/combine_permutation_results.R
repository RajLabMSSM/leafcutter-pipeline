# aggregate permutation results together
# report number of significantly altered clusters at FDR < 0.05
# for each permutation

library(dplyr)
library(readr)
library(purrr)
library(optparse)

option_list <- list(
    make_option(c('--dataCode'), help='', default = "example"),
        make_option(c('--outFolder'), help='', default = "results/example/permutation/")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)


dataCode <- opt$dataCode
outFolder <- opt$outFolder

message("aggregating together permutation results")

#outFolder <- "results/example/permutation/"
outFile <- paste0(outFolder, dataCode, "_all_permutation_results.tsv")

# gather significance files
# for each - read in and count number of sig clusters at FDR < 0.05
res_files <- list.files(outFolder, pattern = "*significance.txt", full.names = TRUE )

results <- 
  map_df(1:length(res_files), ~{
	df <- read_tsv(res_files[.x], col_types = "cdddcd")
	n.sig = sum(df$p.adjust < 0.05)
	data.frame( i  = .x, n.sig = n.sig)	
})

write_tsv(x = results, path = outFile)
