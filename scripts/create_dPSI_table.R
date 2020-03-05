#Jack Humphrey
#2019

#options(echo=TRUE)
library(tidyverse)
#library(leafcutter)
library(optparse)

option_list <- list(
    make_option(c('--app'), help='', default = "shiny.RData"),
    make_option(c('--dataCode'), help='', default = "example"),
   make_option(c('--contrast'), help='A contrast string, denoting two conditions separated by contrast_sep', default = "control_case"),
    make_option(c('--contrastSep'), help = "how the two conditions in the contrast string are separated", default = "_" ), 
    make_option(c('--refCondition'), help='', default = "control"),
    make_option(c('--altCondition'), help='', default="case"),
        make_option(c('--outFolder'), help='', default = "leafcutter/")

)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)


app <- opt$app
dataCode <- opt$dataCode
contrast <- opt$contrast
contrastSep <- opt$contrastSep
#refCondition <- opt$refCondition
#altCondition <- opt$altCondition
outFolder <- opt$outFolder


outFile <- paste0(outFolder, dataCode)

load(app)

# for each significant cluster, calculate delta PSI for each junction in cluster
# this requires:
  # counts - counts of junctions in each cluster in each subject
  # meta - which samples are from which group
  # introns_to_plot - gives you which junctions are in which cluster

# so for each cluster:
  # get junction counts for just junctions in cluster
  # normalise junction counts in each sample by total number of counts
  # for each junction take mean junction count in each group (PSI)
  # calculate PSI_PD - PSI_control (deltaPSI)
  # output table of junction deltaPSIs
  # then summarise by taking largest abs(dPSI) for each cluster
# how to treat all-zero cluster counts? for now ignore these samples

deltaPSI <- function(cluster, control_samples, case_samples){
  # testing
  #cluster <- "clu_1_-"
  cluster_counts <- counts[ introns_to_plot$clu == cluster, ]
  # normalise per sample
  norm_counts <- sweep(cluster_counts, MARGIN = 2, STATS = colSums(cluster_counts),FUN = "/")

  
  control_counts <- norm_counts[, names(norm_counts) %in% control_samples]
  case_counts <- norm_counts[, names(norm_counts) %in% case_samples]

  # calculate means
  control_means <- rowMeans(control_counts, na.rm=TRUE)
  case_means <- rowMeans(case_counts, na.rm=TRUE)
  
  # phaedra calculates PSI by taking counts over the median - look into
  
  # calculate deltaPSI
  dPSI <- case_means - control_means
  dPSI_table <- tibble(clusterID = cluster, intron = names(dPSI), control = control_means, case = case_means, dPSI = dPSI)
  
  return(dPSI_table)
}

# get ref and alt condition from contrast string
refCondition <- str_split_fixed(contrast, contrast_sep, 2)[,1]
altCondition <- str_split_fixed(contrast, contrast_sep, 2)[,2]

# get sample IDs - assumes only two groups
control_samples <- meta$sample[meta$group == refCondition]
case_samples <- meta$sample[meta$group == altCondition]

# cluster_ids is vector of clusters
dPSI_table <- purrr::map_df( clusters$clusterID, deltaPSI, control_samples = control_samples, case_samples = case_samples)
#save.image("test.RData")
# the deltaPSI values don't quite match those from leafviz, weirdly enough.
# it looks like the meanPSI values do and the leaviz deltaPSIs aren't matching those
dPSI_full <-
  introns %>%
  mutate(intron = paste(chr, start, end, clusterID, sep = ":") ) %>%
  left_join(dPSI_table, by = c("clusterID", "intron") )  %>%  
  rename(fitted_dPSI = "deltapsi")


# calculate per-cluster highest deltaPSI (abs)
dPSI_best <- 
  group_by(dPSI_table, clusterID) %>%
  summarise(dPSI = dPSI[which.max(abs(dPSI))] )  %>%
  arrange(dPSI) %>%
  left_join(clusters, by = "clusterID") %>%
  mutate(gene = gsub("<i>|</i>", "", gene) )


# write out tables
readr::write_tsv(dPSI_full, path = paste0(outFile, "_deltapsi_full.tsv" ) )
readr::write_tsv(dPSI_best, path = paste0(outFile, "_deltapsi_best.tsv" ) )

