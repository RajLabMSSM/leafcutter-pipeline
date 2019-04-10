
# arguments
# samples.tsv
# refCondition
# altCondition
# outFile

library(optparse)

option_list <- list(
    make_option(c('--samples'), help='', default = "samples.tsv"),
    make_option(c('--dataCode'), help='', default = "example"),
    make_option(c('--refCondition'), help='', default = "control"),
    make_option(c('--altCondition'), help='', default="case"),
	make_option(c('--outFolder'), help='', default = "leafcutter/")

)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)


samples <- opt$samples
dataCode <- opt$dataCode
refCondition <- opt$refCondition
altCondition <- opt$altCondition
outFolder <- opt$outFolder

outFile <- paste0(outFolder, "/", dataCode, "_ds_support.tsv")


df <- read.table(samples, header=TRUE, sep = "\t", stringsAsFactors = FALSE)

names(df)[2] <- "condition"

# remove any samples that have a condition not in ref or alt

df <- df[ df$condition %in% c(refCondition, altCondition),]

# leafcutter assumes first string in second column is the reference condition, so reorder table if necessary

df$condition <- factor(df$condition, levels = c(refCondition, altCondition))

df <- df[ order(df$condition),]

# write out

write.table(df, file = outFile, col.names=FALSE, sep = "\t", quote = FALSE, row.names = FALSE)