
# arguments
# samples.tsv
# refCondition
# altCondition
# outFile

library(optparse)
library(stringr)
option_list <- list(
    make_option(c('--metadata'), help='', default = "samples.tsv"),
    make_option(c('--dataCode'), help='', default = "example"),
    make_option(c('--refCondition'), help='', default = "control"),
    make_option(c('--altCondition'), help='', default="case"),
	make_option(c('--contrast'), help='A contrast string, denoting two conditions separated by contrast_sep', default = "control_case"),
    make_option(c('--contrastSep'), help = "how the two conditions in the contrast string are separated", default = "_" ),
    make_option(c('--outFolder'), help='', default = "leafcutter/"),
	make_option(c('--junctionMode'), help = 'whether junctions are from regtools or from RAPiD', default = "RAPiD"),
	make_option(c('--juncSuffix'), help='suffix after sample name to name junctions', default = ''),
	make_option(c('--permutation'), action = "store_true", help = "whether to permute condition column", default = FALSE)
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)


metadata <- opt$metadata
dataCode <- opt$dataCode
contrast <- opt$contrast
contrast_sep <- opt$contrastSep
refCondition <- opt$refCondition
altCondition <- opt$altCondition
outFolder <- opt$outFolder
juncSuffix <- opt$juncSuffix
junctionMode <- opt$junctionMode
permutation <- opt$permutation

outFile <- paste0(outFolder, "/", contrast, "/", dataCode, "_", contrast, "_ds_support.tsv")

message("Creating leafcutter support file from provided metadata")
message("-------------------------------------------------------")

if( ! file.exists(metadata ) ){ stop("metadata does not exist!" ) }


df <- read.table(metadata, header=TRUE, sep = "\t", stringsAsFactors = FALSE)

# drop rapid_path if present

if( "rapid_path" %in% names(df) ){
df$rapid_path <- NULL
}

if( names(df)[2] != "condition" ){
	stop("the second column must be \"condition\" ")
}

# remove any samples that have a condition not in ref or alt
refCondition <- str_split_fixed(contrast, contrast_sep, 2)[,1]
altCondition <- str_split_fixed(contrast, contrast_sep, 2)[,2]



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


# if junctions are from regtools, then the sample IDs will be consistent in the junction counts matrix
# if they are from RAPiD, the sample IDs in the counts matrix will have the bamSuffix included.

if( junctionMode == "RAPiD"){
juncSuffix <- gsub('\\.junc', '', juncSuffix)

df$sample <- paste0(df$sample, juncSuffix)
}

# permutation

if(permutation == TRUE){
message("permutation mode selected")
message("condition column will be shuffled")

df$condition <- sample(df$condition)

}

# write out

print(df)

message("creating support file")
message(paste0("saving as ", outFile) )
write.table(df, file = outFile, col.names=FALSE, sep = "\t", quote = FALSE, row.names = FALSE)
