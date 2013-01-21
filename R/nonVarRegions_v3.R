# nonVarRegions_v3.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


nonVarRegions_v3 <- function(
  regionsV3bedFilename        = "/data/malariagen2/plasmodium/pf-crosses/meta/regions_v3_rdp.bed"
) {
  regionsBedV3 <- read.table(regionsV3bedFilename, header=FALSE, sep="\t", as.is=TRUE)
  regions_v3GR <- GRanges(
    seqnames=regionsBedV3[[1]],
    ranges=IRanges(
      start=regionsBedV3[[2]],
      end=regionsBedV3[[3]]
    ),
    region=regionsBedV3[[4]]
  )
  nonVarRegions <- regions_v3GR[values(regions_v3GR)[["region"]] == "O"]
  return(nonVarRegions)
}

