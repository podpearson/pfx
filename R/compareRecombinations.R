# compareRecombinations.R
# 
# Package: pfx
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


compareRecombinations <- function(
  subjectRecombinations,
  comparisonRecombinations,
  plotFilestem                = "7G8xGB4vsJiangComparison"
) {
#  subjectRecombinationsList <- sapply(subjectRecombinations, function(x) x[["recombinationsGRL"]], USE.NAMES=TRUE, simplify=FALSE)
#  subjectRecombinationsUncertainties <- sapply(
#    subjectRecombinations,
#    function(x) {
#      unlist(x[["widthsList"]])
#    }
#  )
#  stem(log10(unlist(subjectRecombinationsUncertainties)))
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
#  allSubjectRecombinations <- GenomicRanges::unlist(subjectRecombinationsBySampleGRlist)
  allSubjectRecombinations <- GenomicRanges::unlist(subjectRecombinations[["recombinationsBySampleGRlist"]])

#  comparisonRecombinationsList <- sapply(comparisonRecombinations, function(x) x[["recombinationsGRL"]], USE.NAMES=TRUE, simplify=FALSE)
#  comparisonRecombinationsBySample <- sapply(
#    names(comparisonRecombinationsList[[1]]),
#    function(sampleID) {
#        lapply(
#          comparisonRecombinationsList,
#          function(x) {
#            x[[sampleID]]
#          }
#        )
#    },
#    simplify=FALSE,
#    USE.NAMES=TRUE
#  )
#  comparisonRecombinationsBySampleGRlist <- GRangesList(
#    sapply(
#      names(comparisonRecombinationsList[[1]]),
#      function(sampleID) {
#        GenomicRanges::unlist(do.call(GRangesList, comparisonRecombinationsBySample[[sampleID]][sapply(comparisonRecombinationsBySample[[sampleID]], length) > 0]))
#      },
#      simplify=FALSE,
#      USE.NAMES=TRUE
#    )
#  )
#  allComparisonRecombinations <- GenomicRanges::unlist(comparisonRecombinationsBySampleGRlist)
  allComparisonRecombinations <- GenomicRanges::unlist(comparisonRecombinations[["recombinationsBySampleGRlist"]])
  
  values(allSubjectRecombinations)[["numberOfOverlappingCrossovers"]] <- countOverlaps(allSubjectRecombinations, allSubjectRecombinations, type="equal")
  values(allComparisonRecombinations)[["numberOfOverlappingCrossovers"]] <- countOverlaps(allComparisonRecombinations, allComparisonRecombinations, type="equal")

  plotDF <- rbind(
    data.frame(dataset="Malariagen", numberOfOverlappingCrossovers=values(allSubjectRecombinations)[["numberOfOverlappingCrossovers"]]),
    data.frame(dataset="Jiang", numberOfOverlappingCrossovers=values(allComparisonRecombinations)[["numberOfOverlappingCrossovers"]])
  )
  require(ggplot2)
  pdf(paste(plotFilestem, "overlappingRecombinationPointsByCallset.pdf", sep="."), height=6, width=8)
  print(
    qplot(numberOfOverlappingCrossovers, facets=dataset~., data=plotDF, binwidth=1, xlab="Number of samples with recombination", ylab="Frequency (number of recombinations)")
    + theme_bw()
  )
  dev.off()
  
  plotDF4 <- rbind(
    data.frame(dataset="Malariagen", uncertaintyInCrossoverPosition=width(allSubjectRecombinations)),
    data.frame(dataset="Jiang", uncertaintyInCrossoverPosition=width(allComparisonRecombinations))
  )
  pdf(paste(plotFilestem, "uncertaintyInRecombinationPoints.pdf", sep="."), height=6, width=8)
  print(
    qplot(uncertaintyInCrossoverPosition, fill=dataset, facets=dataset~., data=plotDF4, binwidth=0.2, log="x", xlab="Uncertainty in recombination position (bp)", ylab="Frequency (number of recombinations)")
  + theme_bw()
  )
  dev.off()
  
  subjectAll <- GenomicRanges::unlist(subjectRecombinations[["recombinationsBySampleGRlist"]])
  subjectFusionGenes <- subjectAll[which(values(subjectAll)[["startForwardGene"]]==values(subjectAll)[["endForwardGene"]] & values(subjectAll)[["startForwardGene"]]!="")]
  subjectFusionGenes <- subjectFusionGenes[GenomicRanges::order(subjectFusionGenes)]
  comparisonAll <- GenomicRanges::unlist(comparisonRecombinations[["recombinationsBySampleGRlist"]])
  comparisonFusionGenes <- comparisonAll[which(values(comparisonAll)[["startForwardGene"]]==values(comparisonAll)[["endForwardGene"]] & values(comparisonAll)[["startForwardGene"]]!="")]
  comparisonFusionGenes <- comparisonFusionGenes[GenomicRanges::order(comparisonFusionGenes)]
  print(length(subjectFusionGenes))
  print(length(unique(values(subjectFusionGenes)[["startForwardGene"]])))
  print(unique(values(subjectFusionGenes)[["startForwardGene"]]))
  print(length(comparisonFusionGenes))
  print(length(unique(values(comparisonFusionGenes)[["startForwardGene"]])))
  print(unique(values(comparisonFusionGenes)[["startForwardGene"]]))
  comparisonFusionGenes[GenomicRanges::match(subjectFusionGenes, comparisonFusionGenes)[!is.na(GenomicRanges::match(subjectFusionGenes, comparisonFusionGenes))]]
  subjectFusionGenes[GenomicRanges::match(comparisonFusionGenes, subjectFusionGenes)[!is.na(GenomicRanges::match(comparisonFusionGenes, subjectFusionGenes))]]
  
#  browser()
  values(allSubjectRecombinations)[["minimumHapLength"]] <- pmin(
    values(allSubjectRecombinations)[["midpoint"]]-values(allSubjectRecombinations)[["nearestMidpointLeft"]],
    values(allSubjectRecombinations)[["nearestMidpointRight"]]-values(allSubjectRecombinations)[["midpoint"]],
    na.rm = TRUE
  )
  plotDFTemp <- GenomicRanges::as.data.frame(values(allSubjectRecombinations))
  plotDFTemp[is.na(plotDFTemp[["minimumHapLength"]]), "minimumHapLength"] <- 1
  pdf(paste(plotFilestem, "minHapLengthByNumberOfCrossovers.pdf", sep="."), height=14, width=6)
  print(
    qplot(minimumHapLength, fill=numberOfOverlappingCrossovers, facets=numberOfOverlappingCrossovers~., data=plotDFTemp, log="x", xlab="Length of shortest haplotype created by recombination", ylab="Frequency (number of recombinations)")
  + geom_vline(xintercept = 15000, colour="red")
  + theme_bw()
  )
  dev.off()
  
  allSubjectRecombinations[values(allSubjectRecombinations)[["numberOfOverlappingCrossovers"]]==22]
  values(allSubjectRecombinations)[["suspicious"]] <- with(GenomicRanges::as.data.frame(values(allSubjectRecombinations)), numberOfOverlappingCrossovers > 1 & (minimumHapLength <= 30000 | is.na(nearestMidpointLeft) | is.na(nearestMidpointRight)))
  plotDFTemp2 <- GenomicRanges::as.data.frame(values(allSubjectRecombinations[!values(allSubjectRecombinations)[["suspicious"]]]))
  pdf(paste(plotFilestem, "minHapLengthByNumberOfCrossoversAfterSuspiciousRemoval.pdf", sep="."), height=14, width=6)
  print(
    qplot(minimumHapLength, fill=numberOfOverlappingCrossovers, facets=numberOfOverlappingCrossovers~., data=plotDFTemp2, log="x", xlab="Length of shortest haplotype created by recombination", ylab="Frequency (number of recombinations)")
  + geom_vline(xintercept = 15000, colour="red")
  + theme_bw()
  )
  dev.off()
  
#  allSubjectRecombinations[!values(allSubjectRecombinations)[["suspicious"]] & values(allSubjectRecombinations)[["numberOfOverlappingCrossovers"]]==2]
#  table(as.character(seqnames(allSubjectRecombinations)), useNA="ifany")
#  seqlevels(allSubjectRecombinations) <- paste("MAL", 1:14, sep="")
#  table(as.character(seqnames(allSubjectRecombinations)), useNA="ifany")
#  suspiciousRecombinationPoints <- BiocGenerics::unique(allSubjectRecombinations[values(allSubjectRecombinations)[["suspicious"]]])[BiocGenerics::order(BiocGenerics::unique(allSubjectRecombinations[values(allSubjectRecombinations)[["suspicious"]]]))]
#  suspiciousRecombinationPoints[1:19]
#  suspiciousRecombinationPoints[20:38]
#  suspiciousRecombinationPoints[39:47]
#  suspiciousRecombinationPoints[values(suspiciousRecombinationPoints)[["numberOfOverlappingCrossovers"]] > 3]
#  GenomicRanges::reduce(suspiciousRecombinationPoints[values(suspiciousRecombinationPoints)[["numberOfOverlappingCrossovers"]] >= 3])
#  
#  suspiciousRegions <- c(
#    GRanges("MAL6",  IRanges(1,       74373),   seqinfo=seqinfo(allSubjectRecombinations)),
#    GRanges("MAL9",  IRanges(1214920, 1241458), seqinfo=seqinfo(allSubjectRecombinations)),
#    GRanges("MAL11", IRanges(1902245, 1935028), seqinfo=seqinfo(allSubjectRecombinations)),
#    GRanges("MAL12", IRanges(763134,  776156),  seqinfo=seqinfo(allSubjectRecombinations)),
#    GRanges("MAL13", IRanges(1578814, 1591642), seqinfo=seqinfo(allSubjectRecombinations)),
#    GRanges("MAL14", IRanges(1,       35424),   seqinfo=seqinfo(allSubjectRecombinations)),
#    GRanges("MAL14", IRanges(3193783, 3199050), seqinfo=seqinfo(allSubjectRecombinations))
#  )
#  additionalSuspiciousRegion <- GRanges("MAL14", IRanges(237968, 240153), seqinfo=seqinfo(allSubjectRecombinations))
  
  
  return(
    c(
      medianWidthInSubject = median(width(allSubjectRecombinations)),
      medianWidthInComparison = median(width(allComparisonRecombinations))
    )
  )
}
