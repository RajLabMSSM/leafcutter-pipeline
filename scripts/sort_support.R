
# arguments
# samples.tsv
# refCondition
# altCondition
# outFile

library(optparse)

option_list <- list(
    make_option(c('--metadata'), help='', default = "samples.tsv"),
    make_option(c('--dataCode'), help='', default = "example"),
    make_option(c('--refCondition'), help='', default = "control"),
    make_option(c('--altCondition'), help='', default="case"),
	make_option(c('--outFolder'), help='', default = "leafcutter/")

)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)


metadata <- opt$metadata
dataCode <- opt$dataCode
refCondition <- opt$refCondition
altCondition <- opt$altCondition
outFolder <- opt$outFolder

outFile <- paste0(outFolder, "/", dataCode, "_ds_support.tsv")

message("Creating leafcutter support file from provided metadata")
message("-------------------------------------------------------")

if( ! file.exists(metadata ) ){ stop("metadata does not exist!" ) }


df <- read.table(metadata, header=TRUE, sep = "\t", stringsAsFactors = FALSE)

if( names(df)[2] != "condition" ){
	stop("the second column must be \"condition\" ")
}

# remove any samples that have a condition not in ref or alt

if( ! all( c(refCondition, altCondition) %in% df$condition) ){
	stop("neither condition is in metadata")
}


df <- df[ df$condition %in% c(refCondition, altCondition),]

# leafcutter assumes first string in second column is the reference condition, so reorder table if necessary


df$condition <- factor(df$condition, levels = c(refCondition, altCondition))

df <- df[ order(df$condition),]

if(nrow(df) == 0){
	stop("table is empty")
}

print(df)

# write out

message("creating support file")
message(paste0("saving as ", outFile) )
write.table(df, file = outFile, col.names=FALSE, sep = "\t", quote = FALSE, row.names = FALSE)
