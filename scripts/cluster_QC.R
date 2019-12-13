# junction clustering QC
library(data.table)
library(leafcutter)
library(dplyr)
library(purrr)
library(ggplot2)
library(optparse)
library(igraph)


readClustersFile <- function(file){
  clusters <- read.table(file, header=TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  clusters <- tibble::rownames_to_column(clusters, var = "junctionID")
  #clusters$junctionID <- row.names(clusters)
  #names(clusters)[1] <- "junctionID"
  
  meta <- leafcutter::get_intron_meta(clusters$junctionID)
  meta$ID <- clusters$junctionID
  clusters$clusterID <- meta$clu
  
  clusters <- as.data.table(clusters)
  
  return(clusters)
}


writeClusterFile <- function(clusters, outFile){
  message( paste("writing filtered clusters to", outFile))
  to_write <- clusters
  clusters$clusterID <- NULL
  rownames(clusters) <- clusters$junctionID
  clusters$junctionID <- NULL
  write.table(clusters, file = gzfile(outFile), quote = FALSE,sep = " ",row.names = TRUE, col.names = TRUE,  )
}


# use data.table to index junction count matrix by cluster ID


# first index by cluster
#setkey(clusters, "clusterID")

# count up numbers of junctions per cluster

#introns that are supported by at least 10% of the total number of reads supporting the clusters 
# they belong to in at least 25% of samples, considering each age separately. Clusters

junctionStats <- function(clusters, label = "."){
  n_junc <- length(unique(clusters$junctionID))
  n_cluster <- length(unique(clusters$clusterID))
  
  return(data.frame(label = label, junctions = n_junc, clusters = n_cluster, stringsAsFactors = FALSE))
}


# calculate per junction missingness - how many samples have 0 supporting reads for this junction?
calculateJunctionMissingness <- function(clusters){
  message("Calculating junction missingness...")
  coverage <-  clusters[, -"clusterID"]
  
  missingness_df <- 
    apply( coverage, MARGIN = 1, FUN = function(x){
      junctionID <- x[1]
      df <- as.numeric(x[-1]) # remove junctionID column

      data.frame(
          junctionID = junctionID, 
           missingness = ( sum(df == 0) / length(df)  ),
           stringsAsFactors = FALSE
           )
  }) %>% bind_rows()
    return(missingness_df)
  }

# calculate mean contribution for each junction to its cluster
calculateJunctionRatio <- function(clusters){
  message("Calculating average junction ratio...")
  # sum junctions across all samples
  sums <- data.table(junctionID = clusters$junctionID, clusterID = clusters$clusterID , total = rowSums(clusters[, -c("junctionID", "clusterID")] ) )
  setkey(sums, "clusterID")
  sums_list <- split(sums, by = "clusterID")
  
  ratios <-
    purrr::map_df( sums_list, ~{
      
      df <- .x[, -"clusterID"]
      df$ratio <- df$total / sum(df$total)
      return(df)

    } )
  
  return(ratios)
  
}

# for each cluster, assemble junctions into a graph
# flag disconnected junctions for removal
# if no junctions in cluster are connected then remove cluster




# make graphs from junction meta 
# using graph_from_edgelist
makeClusterGraph <- function(cluster){
  cluster %>%
    select(start,end) %>%
    as.matrix() %>%
    igraph::graph_from_edgelist()
}
  
# is_connected - logical test for connectivity
# if TRUE then do nothing
# if not fully connected
# then decompose() into subgraphs
# retain largest subgraph
# ecount - number of edges
testConnectivity <- function(graph){
  if( ecount(graph) == 1){
    return("remove")
  }
  if( is.connected(graph)){
    return("keep")
  }else{
    # decompose into subgraphs
    decomp <- decompose(graph)
    # count edges in subgraphs
    edge_counts <- map_dbl(decomp, ecount)
    # if all subgraphs are length 1 then remove cluster
    if( all(edge_counts == 1 )){
      return("remove")
    }else{
      return("modify")
    }
  }
}

# for graphs with disconnected junctions
# produce a list of disconnected junctions to remove
# apply on each cluster that has results == modify
junctions_to_prune <- function(graph){
  # decompose and pick the largest subgraph
  decomp <- decompose(graph)
  edge_numbers <- map_dbl(decomp, ecount)
  disconnected <- decomp[ which( edge_numbers != max(edge_numbers) )  ]
  # return to matrix
  map_df( disconnected, ~{ as.data.frame(as_edgelist(.x), stringsAsFactors = FALSE  )  }  )
}


# calculate total size of cluster
calculateClusterSize <- function(clusters){
  message("Calculating cluster sizes...")
  setkey(clusters, "clusterID")
  cluster_list <- split(clusters, by = "clusterID")
  
  cluster_lengths <- 
    purrr::map_df(cluster_list, ~{
    data.frame(n = nrow(.x), stringsAsFactors = FALSE)
    }, .id = "clusterID")
  
  return(cluster_lengths)
  
}

plotJunctionMissingness <- function( missingness_report){
  missingness_report %>%
    ggplot(aes( x = missingness)) +
    geom_histogram(binwidth = 0.1, colour = "white") +
    scale_x_continuous(labels=scales::percent) +
    theme_bw() +
    geom_vline(xintercept = missingness, linetype = 3, colour = "black") +
    stat_bin(binwidth=0.1, geom="text", colour="black", size=3.5,
             aes(label=..count..),
             vjust=-1.5 ) +
    labs(title = "Percent missingness per junction") 
}

plotJunctionRatio <- function(ratio_report){
  ratio_report %>%
    ggplot(aes( x = ratio)) +
    geom_histogram(binwidth = 0.1, colour = "white") +
    scale_x_continuous(labels=scales::percent) +
    theme_bw() +
    geom_vline(xintercept = minratio, linetype = 3, colour = "black") +
    labs(title = "Junction contribution to total cluster reads") +
    stat_bin(binwidth=0.1, geom="text", colour="black", size=3.5,
             aes(label=..count..),
             vjust=-1.5 ) 
}

plotClusterSize <- function(cluster_size_report, label){
  # hard limit at `15 junctions
  cluster_size_report$n[ cluster_size_report$n >= 15] <- 15
  
  cluster_size_report %>%
    ggplot(aes( x = n)) +
    geom_histogram(binwidth = 1) +
    #scale_x_continuous(labels=scales::percent) +
    theme_bw() +
    geom_vline(xintercept = 1 + 0.5, linetype = 3) +
    geom_vline(xintercept = maxsize + 0.5, linetype = 3) +
    labs(title = "Cluster size", subtitle = label) +
    #geom_text(stat='count', aes(label=..count..), position = position_stack(vjust = 1),size=4)
    #geom_histogram(aes(fill=cut), binwidth=1500, colour="grey20", lwd=0.2) +
    stat_bin(binwidth=1, geom="text", colour="black", size=3.5,
             aes(label=..count..),
             vjust=-1.1 ) 
}



## implement QC from Balliu et al
option_list <- list(
  make_option(c('--dataCode'), help='', default = "example"),
  make_option(c('--outFolder'), help='', default = "results/example/"),
  make_option(c('--missingness'), help = 'proportion of total sample size that are missing a junction for filtering: 0-1', default = 0.25),
  make_option(c('--minratio'), help = "the minimum contribution for a junction to its cluster: 0-1", default = 0.05),
  make_option(c('--maxsize'), help = "the maximum number of junctions allowed in a cluster", default = 10)
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)


dataCode <- opt$dataCode
outFolder <- opt$outFolder
missingness <- opt$missingness
minratio <- opt$minratio
maxsize <- opt$maxsize

inFile <- paste0(outFolder, dataCode, "_perind_numers.counts.gz")

stopifnot(file.exists(inFile))

clusters <- readClustersFile(inFile)

## calculate metrics

junction_missingness_baseline <- calculateJunctionMissingness(clusters)
junction_ratio_baseline <- calculateJunctionRatio(clusters)
cluster_size_baseline <- calculateClusterSize(clusters)


## apply junction filters
high_missingness <- junction_missingness_baseline$junctionID[ junction_missingness_baseline$missingness >= missingness]
low_ratio <- junction_ratio_baseline$junctionID[ junction_ratio_baseline$ratio < minratio]

junctions_to_remove <- clusters$junctionID[ clusters$junctionID %in% high_missingness | clusters$junctionID %in% low_ratio ]

message( paste("Removing ", length(junctions_to_remove), "junctions that fail QC") )

if(length(junctions_to_remove) == 0){
  junctions_post_qc <- clusters
}else{
  junctions_post_qc <- clusters[ ! clusters$junctionID %in% junctions_to_remove]
}


if(nrow(junctions_post_qc) == 0){
  stop("No junctions remain after junction QC filtering. consider relaxing your thresholds")
}

## Graph QC - check for connectivity

meta <- leafcutter::get_intron_meta(junctions_post_qc$junctionID)
meta$start <- as.character(meta$start)
meta$end <- as.character(meta$end)
meta_by_cluster <- split(meta, meta$clu)

cluster_graphs <- map(meta_by_cluster, makeClusterGraph)
results <- map(cluster_graphs, testConnectivity )

# when results == "modify" then prune the disconnected junctions from the cluster
to_prune <- cluster_graphs[ results == "modify" ] 


prune_df <- map_df(to_prune, junctions_to_prune )
names(prune_df) <- c("start", "end")
#prune_df$start <- as.numeric(prune_df$start)
#prune_df$end <- as.numeric(prune_df$end)

# match junctions in meta to prune list
meta$to_prune <- ifelse( meta$start %in% prune_df$start & meta$end %in% prune_df$end, TRUE, FALSE)
# remove these junctions from junction_post_qc
junctions_pruned <- junctions_post_qc[ meta$to_prune == FALSE,] 

meta_pruned <- leafcutter::get_intron_meta(junctions_pruned$junctionID)

# when results == remove, remove entire cluster
clusters_to_remove <- names(results)[results == "remove" ]

meta_pruned$to_remove <- meta_pruned$clu %in% clusters_to_remove

junctions_removed <- junctions_pruned[ meta_pruned$to_remove == FALSE,] 

# recheck connectivity
meta <- leafcutter::get_intron_meta(junctions_removed$junctionID)
meta$start <- as.character(meta$start)
meta$end <- as.character(meta$end)
meta_by_cluster <- split(meta, meta$clu)

cluster_graphs <- map(meta_by_cluster, makeClusterGraph)
results <- map_chr(cluster_graphs, testConnectivity )

# now left with clusters that have to be split into multiple clusters
# append letters to original clusterID
# if clu_001_+ becomes clu_001A_+  and clu_001B_+
splitClusters <- function(graph, clusterID){
  decomp <- decompose(graph)
  clusters_split <- map(decomp, as_edgelist)
  names(clusters_split) <- paste0(clusterID, "_", letters[1:length(clusters_split)] )
  clusters_split <- map2_df( clusters_split,  names(clusters_split), ~{
    df <- as.data.frame(.x, stringsAsFactors = FALSE)
    names(df) <- c("start", "end")
    df$clu <- .y
    df
  } )
  return(clusters_split)
}

# get the graphs that require modification
clusters_to_split <- cluster_graphs[names(results)[results == "modify" ] ]

# split these clusters only
split_clusters <- map2_df(clusters_to_split, names(clusters_to_split), splitClusters)

# add new split cluster ids to meta
# sneaky way - left join the split_clusters to the metadata
# this will create clu.x and clu.y 
# now coalesce the two values with clu.y going first - therefore if clu.y is a new cluster label that will be the label picked
meta_with_split <- left_join(meta, split_clusters, by = c("start", "end")) %>%
  mutate(clu = coalesce(clu.y, clu.x))

junctions_split <- junctions_removed
junctions_split$clusterID <- meta_with_split$clu




# Cluster QC: 
# after junction QC there will now be clusters that either have 1 or junction or are no longer connected together
# flag clusters that are no longer fully connected, or prune away non-connected junctions?

#exclude clusters with only a single remaining junction or greater than 10 junctions
cluster_size_post_qc <- calculateClusterSize(junctions_post_qc)

cluster_size_post_graph_qc <- calculateClusterSize(junctions_split)

#table(cluster_size_post_qc$n > 1 & cluster_size_post_qc$n < maxsize)

clusters_to_keep <- cluster_size_post_graph_qc$clusterID[cluster_size_post_graph_qc$n > 1 & cluster_size_post_graph_qc$n <= maxsize]

setkey(junctions_split, "clusterID")

message( paste("Removing ", nrow(cluster_size_post_graph_qc) - length(clusters_to_keep), "clusters that fail size QC") )

clusters_filtered <- junctions_split[ clusters_to_keep ]

message("Creating QC plots...")

## plots
plotFile <- paste0(outFolder, dataCode, "_junction_qc_plots.pdf")
pdf(plotFile, onefile = TRUE)
plotJunctionMissingness(junction_missingness_baseline)

plotJunctionRatio(junction_ratio_baseline)

plotClusterSize(cluster_size_baseline, label = "pre junction QC")

plotClusterSize(cluster_size_post_qc, label = "post junction QC")

plotClusterSize(cluster_size_post_graph_qc, label = "post graph QC")

dev.off()

# output counts of junctions and clusters
qc_stats <- list(
  junctionStats(clusters, label = "pre_QC"),
  junctionStats(junctions_post_qc, label = "post_junction_QC"),
  junctionStats(junctions_pruned, label = "pruned disconnected junctions"),
  junctionStats(junctions_removed, label = "removed totally disconnected clusters"),
  junctionStats(junctions_split, label = "splitting disconnected clusters"),
  junctionStats(clusters_filtered, label = "final_QC")
  
  
  
) %>%
  bind_rows()




statsFile <- paste0(outFolder, dataCode, "_qc_statistics.txt")

write.table(qc_stats, statsFile, quote = FALSE, row.names = FALSE, sep = "\t")

outFile <- paste0(outFolder, dataCode, "_filtered_perind_numers.counts.gz")

writeClusterFile(clusters_filtered, outFile)

print(qc_stats)



