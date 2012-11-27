# analyseRecombinations.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


analyseRecombinations <- function(
  recombinations,
  plotFilestem                = "7G8xGB4",
  geneConversionLengthCutoff  = 25000,
  geneConversionLengthComeron = 15000
) {
#  subjectRecombinationsList <- sapply(recombinations, function(x) x[["recombinationsGRL"]], USE.NAMES=TRUE, simplify=FALSE)
#  subjectRecombinationsUncertainties <- sapply(
#    recombinations,
#    function(x) {
#      unlist(x[["widthsList"]])
#    }
#  )
#  
#  subjectRecombinationsBySample <- sapply(
#    names(subjectRecombinationsList[[1]]),
#    function(sampleID) {
#        lapply(
#          subjectRecombinationsList,
#          function(x) {
#            x[[sampleID]]
#          }
#        )
#    },
#    simplify=FALSE,
#    USE.NAMES=TRUE
#  )
#  subjectRecombinationsBySampleGRlist <- GRangesList(
#    sapply(
#      names(subjectRecombinationsList[[1]]),
#      function(sampleID) {
#        GenomicRanges::unlist(do.call(GRangesList, subjectRecombinationsBySample[[sampleID]][sapply(subjectRecombinationsBySample[[sampleID]], length) > 0]))
#      },
#      simplify=FALSE,
#      USE.NAMES=TRUE
#    )
#  )

  allRecombinations <- GenomicRanges::unlist(recombinations[["recombinationsBySampleGRlist"]])
  values(allRecombinations)[["haplotypeLengthToLeft"]] <- ifelse(
    is.na(values(allRecombinations)[["nearestMidpointLeft"]]),
    values(allRecombinations)[["midpoint"]],
    values(allRecombinations)[["midpoint"]] - values(allRecombinations)[["nearestMidpointLeft"]]
  )
  values(allRecombinations)[["haplotypeLengthToRight"]] <- ifelse(
    is.na(values(allRecombinations)[["nearestMidpointRight"]]),
    seqlengths(allRecombinations)[as.character(seqnames(allRecombinations))] - values(allRecombinations)[["midpoint"]],
    values(allRecombinations)[["nearestMidpointRight"]] - values(allRecombinations)[["midpoint"]]
  )
    
  haplotypeLengths <- c(
    values(allRecombinations)[["haplotypeLengthToLeft"]],
    values(allRecombinations)[["haplotypeLengthToRight"]][is.na(values(allRecombinations)[["nearestMidpointRight"]])]
  ) 

  
  haplotypeLengthsGaps <- BiocGenerics::sapply(recombinations[["recombinationsBySampleGRlist"]], function(x) width(gaps(x)))
  
  require(ggplot2)
  pdf(paste(plotFilestem, "haplotypeLengthsGaps.pdf", sep="."), height=4, width=8)
  print(
    qplot(unlist(haplotypeLengthsGaps)+2, log="x", xlab="Length of haplotypes (bp)", ylab="Frequency (number of haplotypes)")
  + geom_vline(xintercept = geneConversionLengthCutoff, colour="red")
  + theme_bw()
  )
  dev.off()
  
  pdf(paste(plotFilestem, "haplotypeLengths.pdf", sep="."), height=4, width=8)
  print(
    qplot(haplotypeLengths, log="x", xlab="Length of haplotypes (bp)", ylab="Frequency (number of haplotypes)")
  + geom_vline(xintercept = geneConversionLengthComeron, colour="red")
  + geom_vline(xintercept = geneConversionLengthCutoff, colour="blue")
  + theme_bw()
  )
  dev.off()
  
  values(allRecombinations)[["isCrossoverEvent"]] <- values(allRecombinations)[["haplotypeLengthToLeft"]] >= geneConversionLengthCutoff & values(allRecombinations)[["haplotypeLengthToRight"]] >= geneConversionLengthCutoff
  fusionGenes <- allRecombinations[which(values(allRecombinations)[["startForwardGene"]]==values(allRecombinations)[["endForwardGene"]] & values(allRecombinations)[["startForwardGene"]]!="")]
  fusionGenes <- fusionGenes[GenomicRanges::order(fusionGenes)]
  
  
  cat("There are", length(allRecombinations), "recombination points.")
  cat("Of these", length(fusionGenes), "are within", length(unique(values(fusionGenes)[["startForwardGene"]])), "unique genes.")
  cat("Of these", length(which(values(fusionGenes)[["isCrossoverEvent"]])), "are crossover events")
  
  

}
