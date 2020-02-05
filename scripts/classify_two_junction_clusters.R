# Classify 2-junction clusters

# Jack Humphrey
# 2020


load(shiny_results.RData)

meta <- leafcutter::get_intron_meta(clusters$clusterID)

clusters <- clusters[ clusters$N == 2, ]

exons <- ""
UTR5 <- ""
UTR3 <- ""

# match all junctions to exons, 3'UTR and 5'UTR
# split by strand first
by_strand <- split(junctions, junctions$strand)

pos <- by_strand$`+`
neg <- by_strand$`-`

pos$exon <- paste(pos$chr, pos$end) %in% paste(exon$chr, exon$start) | paste(pos$chr, post$start) %in% paste(exon$chr, exon$end) 
pos$utr3 <- paste(pos$chr, pos$end) %in% paste(utr3$chr, utr3$start) 
pos$utr5 <- paste(pos$chr, pos$start) %in% paste(utr5$chr, utr5$end)

neg$exon <- paste(neg$chr, neg$end) %in% paste(exon$chr, exon$start) | paste(neg$chr, negt$start) %in% paste(exon$chr, exon$end)
neg$utr3 <- paste(neg$chr, neg$start) %in% paste(utr3$chr, utr3$end)  
neg$utr5 <- paste(neg$chr, neg$end) %in% paste(utr5$chr, utr5$start)

annotated_junctions <- rbind(pos,neg)

by_cluster <- split(annotated_junctions, annotated_junctions$clusterID)

classify_2_junction_cluster <- function(clu){
    # figure out parent and child
    widths <- clu$end - clu$start
    parent <- clu[ which(widths == max(widths) ),]
    child <-  clu[ which(widths == min(widths) ),]
    # find alignment - left or right
    if(length(unique(clu$start))  == 1 ){
        alignment <- "left"
    }else{
        alignment <- "right"
    }
    strand <- clu$strand

    verdict <- NULL
    #   if cluster$orientation == "left_anchored" & cluster$strand == "+" | cluster$orientation == "right_anchored" & cluster$strand == "-":
#       if parent is 3'UTR and child is 3'UTR and overlap = FALSE; verdict = alternate 3'end
#       if parent is exon and child is exon and overlap = TRUE; verdict = alternate 3'SS
#   if cluster$orientation == "right_anchored" & cluster$strand == "+" | cluster$orientation == "left_anchored" & cluster$strand == "-":
#       if parent is 5'UTR exon and child is 5'UTR exon and the two don't overlap: verdict = alternate start site
#       if parent is exon and child is exon and there is overlap  = alternate 5' splice site    

    if( strand == "+" & aligment == "left" | strand == "-" & alignment == "right" ){
       if( parent$utr3 == TRUE & child$utr3 == TRUE & overlap == FALSE ){
        verdict <- "alt_3_end"    
    }
   #etc 


    return(verdict)
}




# extends classify_clusters.R to look at 2 junction clusters
# these can be 1 of 5 tractable classes:
# alternate start site
# alternate end site
# alternate 5' splice site
# alternate 3' splice site

# this requires lists of known exons and known 3'UTRs and 5'UTRs

# for each 2 junction cluster:

# longer junction is parent, shorter is child



# for each junction 

# get starts and ends and strand

# match starts and ends to exons, 5'UTRs and 3'UTRs

# test if start of child and start of parent overlap in a 5'UTR or exon
# test if end of child and end of parent overlap in a 3'UTR or exon

#   if cluster$orientation == "left_anchored" & cluster$strand == "+" | cluster$orientation == "right_anchored" & cluster$strand == "-":
#       if parent is 3'UTR and child is 3'UTR and overlap = FALSE; verdict = alternate 3'end
#       if parent is exon and child is exon and overlap = TRUE; verdict = alternate 3'SS
#   if cluster$orientation == "right_anchored" & cluster$strand == "+" | cluster$orientation == "left_anchored" & cluster$strand == "-":
#       if parent is 5'UTR exon and child is 5'UTR exon and the two don't overlap: verdict = alternate start site
#       if parent is exon and child is exon and there is overlap  = alternate 5' splice site


