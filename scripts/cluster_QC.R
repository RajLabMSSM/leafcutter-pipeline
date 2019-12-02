# junction clustering QC
library(data.table)
library(leafcutter)
library(dplyr)
library(purrr)
library(ggplot2)

## implement QC from Balliu et al

readClustersFile <- function(file){
  clusters <- fread(file)
  names(clusters)[1] <- "junctionID"
  
  meta <- leafcutter::get_intron_meta(clusters$junctionID)
  meta$ID <- clusters$junctionID
  clusters$clusterID <- meta$clu
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
  
  return(data.frame(label = label, junctions = n_junc, clusters = n_cluster))
}


# calculate per junction missingness - how many samples have 0 supporting reads for this junction?
calculateJunctionMissingness <- function(clusters){
  coverage <-  clusters[, -"clusterID"]
  
  missingness_df <- 
    apply( coverage, MARGIN = 1, FUN = function(x){
      junctionID <- x[1]
      df <- as.numeric(x[-1]) # remove junctionID column

      data.frame(
          junctionID = junctionID, 
           missingness = (1 - ( sum(df > 0) / length(df) ) ),
           stringsAsFactors = FALSE
           )
  }) %>% bind_rows()
    return(missingness_df)
  }

# calculate mean contribution for each junction to its cluster
calculateJunctionRatio <- function(clusters){
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

# calculate total size of cluster
calculateClusterSize <- function(clusters){
  
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
    #geom_vline(xintercept = 0.25, linetype = 1, colour = "red") +
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
    #geom_vline(xintercept = 0.05, linetype = 1, colour = "red") +
    labs(title = "Junction contribution to total cluster reads") +
    stat_bin(binwidth=0.1, geom="text", colour="black", size=3.5,
             aes(label=..count..),
             vjust=-1.5 ) 
}

plotClusterSize <- function(cluster_size_report){
  cluster_size_report %>%
    ggplot(aes( x = n)) +
    geom_histogram(binwidth = 5) +
    #scale_x_continuous(labels=scales::percent) +
    theme_bw() +
    #geom_vline(xintercept = 10, linetype = 3) +
    labs(title = "Cluster size") +
    #geom_text(stat='count', aes(label=..count..), position = position_stack(vjust = 1),size=4)
    #geom_histogram(aes(fill=cut), binwidth=1500, colour="grey20", lwd=0.2) +
    stat_bin(binwidth=5, geom="text", colour="black", size=3.5,
             aes(label=..count..),
             vjust=-1.1 )
  
}

test_file <- "data/FTD_FCX_mod2_perind_numers.counts.gz"

clusters <- readClustersFile(test_file)

# calculate metrics
junction_missingness_baseline <- calculateJunctionMissingness(clusters)
junction_ratio_baseline <- calculateJunctionRatio(clusters)
cluster_size_baseline <- calculateClusterSize(clusters)


# plots and tables

table(junction_missingness_baseline$missingness >= 0.25)
table(junction_ratio_baseline$ratio < 0.05)

## apply junction filters

high_missingness <- junction_missingness_baseline$junctionID[ junction_missingness_baseline$missingness >= 0.25]
low_ratio <- junction_ratio_baseline$junctionID[ junction_ratio_baseline$ratio < 0.05]

junctions_to_remove <- clusters$junctionID[ clusters$junctionID %in% high_missingness | clusters$junctionID %in% low_ratio ]

junctions_post_qc <- clusters[ ! clusters$junctionID %in% junctions_to_remove]

# Cluster QC: exclude clusters with only a single remaining junction or greater than 10 junctions

cluster_size_post_qc <- calculateClusterSize(junctions_post_qc)

table(cluster_size_cleaned$n > 1 & cluster_size_cleaned$n < 10)

clusters_to_keep <- cluster_size_cleaned$clusterID[cluster_size_cleaned$n > 1 & cluster_size_cleaned$n <= 10]

setkey(junctions_clean, "clusterID")

clusters_clean <- junctions_clean[ clusters_to_keep ]

## plots

plotJunctionMissingness(junction_missingness_baseline)

plotJunctionRatio(junction_ratio_baseline)

plotClusterSize(cluster_size_baseline)

plotClusterSize(cluster_size_post_qc)


# output counts of junctions and clusters
junctionStats(clusters, label = "pre-QC")

junctionStats(junctions_clean, label = "post-junction QC")

junctionStats(clusters_clean, label = "post-cluster QC")





