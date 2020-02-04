# Classify 2-junction clusters

# Jack Humphrey
# 2020

# extends classify_clusters.R to look at 2 junction clusters
# these can be 1 of 5 tractable classes:
# alternate start site
# alternate end site
# alternate 5' splice site
# alternate 3' splice site

# this requires lists of known exons and known 3'UTRs and 5'UTRs

# for each 2 junction cluster:

# get starts and ends and strand

# match starts and ends to exons, 5'UTRs and 3'UTRs